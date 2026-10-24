// SPDX-License-Identifier: GPL-2.0-or-later
#include <nanobind/nanobind.h>

#include <array>
#include <charconv>
#include <cstdlib>
#include <cstring>

#include "ducc0/infra/simd.h"
#include "ducc0/sht/sht.h"
#include "internal.hpp"

namespace SPHARM_CPP_NAMESPACE
{
namespace
{

template <typename T>
double coefficient_sign(std::size_t m)
{
  return (m % 2 != 0) ? -1.0 : 1.0;
}

template <typename T>
std::complex<T> scale_coefficient(const std::complex<T>& value, double factor)
{
  return {static_cast<T>(static_cast<double>(value.real()) * factor),
          static_cast<T>(static_cast<double>(value.imag()) * factor)};
}

std::vector<std::size_t> ducc_mstart(std::size_t lmax)
{
  std::vector<std::size_t> result(lmax + 1);
  std::size_t start = 0;
  for (std::size_t m = 0; m <= lmax; ++m)
  {
    result[m] = start - m;
    start += lmax - m + 1;
  }
  return result;
}

std::size_t transform_nthreads()
{
  static const std::size_t configured = []
  {
    const char* value = std::getenv("DUCC0_NUM_THREADS");
    if (!value)
      value = std::getenv("OMP_NUM_THREADS");
    if (!value)
      return std::size_t(1);
    const char* end = value + std::strlen(value);
    std::size_t parsed = 0;
    const auto result = std::from_chars(value, end, parsed);
    if (result.ec != std::errc{} || result.ptr != end)
      throw std::runtime_error(
          "DUCC0_NUM_THREADS/OMP_NUM_THREADS must be a nonnegative integer");
    return parsed;
  }();
  return configured;
}

template <typename T>
DUCC0_NAMESPACE::cmav<T, 3> as_ducc_map(const CNpArr& input, std::size_t nlat,
                                        std::size_t nlon, std::size_t field)
{
  const auto st = DUCC0_NAMESPACE::detail_pybind::copy_strides<T, false>(input);
  std::array<std::size_t, 3> shape{1, nlat, nlon};
  std::array<std::ptrdiff_t, 3> stride{0, st[0], st[1]};
  const auto* ptr = reinterpret_cast<const T*>(input.data());
  if (input.ndim() == 3)
    ptr += static_cast<std::ptrdiff_t>(field) * st[2];
  return DUCC0_NAMESPACE::cmav<T, 3>(ptr, shape, stride);
}

template <typename T>
DUCC0_NAMESPACE::vmav<T, 3> as_ducc_map(NpArr& output, std::size_t nlat,
                                        std::size_t nlon, std::size_t field,
                                        bool batched)
{
  auto storage = DUCC0_NAMESPACE::detail_pybind::to_vfmav<T>(output);
  std::array<std::size_t, 3> shape{1, nlat, nlon};
  std::array<std::ptrdiff_t, 3> stride{0, storage.stride(0), storage.stride(1)};
  auto* ptr = storage.data();
  if (batched)
    ptr += field * storage.stride(2);
  return DUCC0_NAMESPACE::vmav<T, 3>(ptr, shape, stride);
}

template <typename T>
DUCC0_NAMESPACE::cmav<std::complex<T>, 2> as_ducc_alms(const CNpArr& input,
                                                       std::size_t nalm,
                                                       std::size_t field)
{
  const auto st =
      DUCC0_NAMESPACE::detail_pybind::copy_strides<std::complex<T>, false>(input);
  std::array<std::size_t, 2> shape{1, nalm};
  std::array<std::ptrdiff_t, 2> stride{0, st[0]};
  auto* ptr = reinterpret_cast<const std::complex<T>*>(input.data());
  if (input.ndim() == 2)
    ptr += static_cast<std::ptrdiff_t>(field) * st[1];
  return DUCC0_NAMESPACE::cmav<std::complex<T>, 2>(ptr, shape, stride);
}

template <typename T>
DUCC0_NAMESPACE::vmav<std::complex<T>, 2> as_ducc_alms(NpArr& output, std::size_t nalm,
                                                       std::size_t field, bool batched)
{
  auto storage = DUCC0_NAMESPACE::detail_pybind::to_vfmav<std::complex<T>>(output);
  std::array<std::size_t, 2> shape{1, nalm};
  std::array<std::ptrdiff_t, 2> stride{0, storage.stride(0)};
  auto* ptr = storage.data();
  if (batched)
    ptr += field * storage.stride(1);
  return DUCC0_NAMESPACE::vmav<std::complex<T>, 2>(ptr, shape, stride);
}

template <typename T>
void convert_ducc_to_py(NpArr& array, std::size_t lmax, std::size_t nfield,
                        std::size_t nlon)
{
  auto view = DUCC0_NAMESPACE::detail_pybind::to_vfmav<std::complex<T>>(array);
  for (std::size_t f = 0; f < nfield; ++f)
    for (std::size_t m = 0; m <= lmax; ++m)
      for (std::size_t l = m; l <= lmax; ++l)
      {
        const std::size_t idx = m * (2 * lmax - m + 3) / 2 + (l - m);
        const auto at = array.ndim() == 1
                            ? static_cast<std::ptrdiff_t>(idx)
                            : static_cast<std::ptrdiff_t>(idx) * view.stride(0) +
                                  static_cast<std::ptrdiff_t>(f) * view.stride(1);
        const T nyquist_factor = (nlon % 2 == 0 && m == nlon / 2) ? T(0.5) : T(1);
        view.data()[at] *=
            static_cast<T>(coefficient_sign<T>(m) * nyquist_factor / sqrt_2pi);
      }
}

template <typename T>
NpArr grid_to_spectral(const CNpArr& grid, std::size_t nlat, std::size_t nlon,
                       const std::string& gridtype, std::size_t lmax)
{
  const std::size_t nfield = grid.ndim() == 2 ? 1 : grid.shape(2);
  const std::size_t nalm = coefficient_count(lmax);
  NpArr result = grid.ndim() == 2 ? new_array<std::complex<T>>({nalm})
                                  : new_array<std::complex<T>>({nalm, nfield});
  std::vector<std::size_t> starts = ducc_mstart(lmax);
  std::vector<double> ringfactor(nlat, 1.0);
  DUCC0_NAMESPACE::cmav<std::size_t, 1> mstart(
      starts.data(), std::array<std::size_t, 1>{starts.size()});
  DUCC0_NAMESPACE::cmav<double, 1> rfactor(
      ringfactor.data(), std::array<std::size_t, 1>{ringfactor.size()});
  for (std::size_t f = 0; f < nfield; ++f)
  {
    auto alms = as_ducc_alms<T>(result, nalm, f, grid.ndim() == 3);
    auto maps = as_ducc_map<T>(grid, nlat, nlon, f);
    {
      py::gil_scoped_release release;
      DUCC0_NAMESPACE::analysis_2d(alms, maps, 0, lmax, mstart, 1,
                                   gridtype == "regular" ? "CC" : "GL", 0.0, rfactor,
                                   transform_nthreads());
    }
  }
  convert_ducc_to_py<T>(result, lmax, nfield, nlon);
  return result;
}

template <typename T>
NpArr spectral_to_grid(const CNpArr& spec, std::size_t nlat, std::size_t nlon,
                       const std::string& gridtype, std::size_t lmax)
{
  const bool batched = spec.ndim() == 2;
  const std::size_t nfield = batched ? spec.shape(1) : 1;
  const std::size_t nalm = coefficient_count(lmax);
  NpArr result =
      batched ? new_array<T>({nlat, nlon, nfield}) : new_array<T>({nlat, nlon});
  NpArr normalized = new_array<std::complex<T>>({nalm, nfield});
  for (std::size_t f = 0; f < nfield; ++f)
  {
    auto source = as_ducc_alms<T>(spec, nalm, f);
    auto dest = as_ducc_alms<T>(normalized, nalm, f, true);
    for (std::size_t m = 0; m <= lmax; ++m)
      for (std::size_t l = m; l <= lmax; ++l)
      {
        const auto idx = m * (2 * lmax - m + 3) / 2 + (l - m);
        dest(0, idx) =
            source(0, idx) * static_cast<T>(sqrt_2pi * coefficient_sign<T>(m));
      }
  }
  std::vector<std::size_t> starts = ducc_mstart(lmax);
  std::vector<double> ringfactor(nlat, 1.0);
  DUCC0_NAMESPACE::cmav<std::size_t, 1> mstart(
      starts.data(), std::array<std::size_t, 1>{starts.size()});
  DUCC0_NAMESPACE::cmav<double, 1> rfactor(
      ringfactor.data(), std::array<std::size_t, 1>{ringfactor.size()});
  for (std::size_t f = 0; f < nfield; ++f)
  {
    auto dest = as_ducc_alms<T>(normalized, nalm, f, true);
    auto maps = as_ducc_map<T>(result, nlat, nlon, f, batched);
    {
      py::gil_scoped_release release;
      DUCC0_NAMESPACE::synthesis_2d(dest, maps, 0, lmax, mstart, 1,
                                    gridtype == "regular" ? "CC" : "GL", 0.0, rfactor,
                                    transform_nthreads(), DUCC0_NAMESPACE::STANDARD);
    }
  }
  return result;
}

template <typename T>
std::pair<NpArr, NpArr> spectral_gradient(const CNpArr& spec, std::size_t nlat,
                                          std::size_t nlon, const std::string& gridtype,
                                          double rsphere)
{
  const bool batched = spec.ndim() == 2;
  const std::size_t nfield = batched ? spec.shape(1) : 1;
  const std::size_t nalm = spec.shape(0);
  const std::size_t lmax = truncation_from_count(nalm);
  NpArr east =
      batched ? new_array<T>({nlat, nlon, nfield}) : new_array<T>({nlat, nlon});
  NpArr north =
      batched ? new_array<T>({nlat, nlon, nfield}) : new_array<T>({nlat, nlon});
  auto eout = DUCC0_NAMESPACE::detail_pybind::to_vfmav<T>(east);
  auto nout = DUCC0_NAMESPACE::detail_pybind::to_vfmav<T>(north);
  const std::size_t npix = nlat * nlon;
  std::vector<T> map_storage(2 * npix);
  std::vector<std::complex<T>> alm_storage(nalm);
  std::vector<std::size_t> starts = ducc_mstart(lmax);
  std::vector<double> ringfactor(nlat, 1.0);
  DUCC0_NAMESPACE::cmav<std::size_t, 1> mstart(
      starts.data(), std::array<std::size_t, 1>{starts.size()});
  DUCC0_NAMESPACE::cmav<double, 1> rfactor(
      ringfactor.data(), std::array<std::size_t, 1>{ringfactor.size()});
  const std::array<std::size_t, 2> alm_shape{1, nalm};
  const std::array<std::size_t, 3> map_shape{2, nlat, nlon};
  const std::array<std::ptrdiff_t, 3> map_stride{static_cast<std::ptrdiff_t>(npix),
                                                 static_cast<std::ptrdiff_t>(nlon), 1};
  for (std::size_t f = 0; f < nfield; ++f)
  {
    auto source = as_ducc_alms<T>(spec, nalm, f);
    for (std::size_t m = 0; m <= lmax; ++m)
      for (std::size_t l = m; l <= lmax; ++l)
      {
        const auto idx = m * (2 * lmax - m + 3) / 2 + (l - m);
        alm_storage[idx] =
            source(0, idx) * static_cast<T>(sqrt_2pi * coefficient_sign<T>(m));
      }
    std::fill(map_storage.begin(), map_storage.end(), T(0));
    if (lmax > 0)
    {
      DUCC0_NAMESPACE::cmav<std::complex<T>, 2> alm(alm_storage.data(), alm_shape);
      DUCC0_NAMESPACE::vmav<T, 3> map(map_storage.data(), map_shape, map_stride);
      py::gil_scoped_release release;
      DUCC0_NAMESPACE::synthesis_2d(alm, map, 1, lmax, mstart, 1,
                                    gridtype == "regular" ? "CC" : "GL", 0.0, rfactor,
                                    transform_nthreads(), DUCC0_NAMESPACE::DERIV1);
    }
    for (std::size_t j = 0; j < nlat; ++j)
      for (std::size_t i = 0; i < nlon; ++i)
      {
        const auto pyidx = batched ? j * nlon * nfield + i * nfield + f : j * nlon + i;
        const auto pix = j * nlon + i;
        eout.data()[pyidx] =
            static_cast<T>(static_cast<double>(map_storage[npix + pix]) / rsphere);
        nout.data()[pyidx] =
            static_cast<T>(-static_cast<double>(map_storage[pix]) / rsphere);
      }
  }
  return {east, north};
}

template <typename T>
std::pair<NpArr, NpArr> analyze_winds(const CNpArr& ugrid, const CNpArr& vgrid,
                                      std::size_t nlat, std::size_t nlon,
                                      const std::string& gridtype, double rsphere,
                                      std::size_t lmax)
{
  const bool batched = ugrid.ndim() == 3;
  const std::size_t nfield = batched ? ugrid.shape(2) : 1;
  const std::size_t nalm = coefficient_count(lmax);
  NpArr vort = batched ? new_array<std::complex<T>>({nalm, nfield})
                       : new_array<std::complex<T>>({nalm});
  NpArr div = batched ? new_array<std::complex<T>>({nalm, nfield})
                      : new_array<std::complex<T>>({nalm});
  auto vout = DUCC0_NAMESPACE::detail_pybind::to_vfmav<std::complex<T>>(vort);
  auto dout = DUCC0_NAMESPACE::detail_pybind::to_vfmav<std::complex<T>>(div);
  const std::size_t npix = nlat * nlon;
  std::vector<T> map_storage(2 * npix);
  std::vector<std::complex<T>> alm_storage(2 * nalm);
  std::vector<std::size_t> starts = ducc_mstart(lmax);
  std::vector<double> ringfactor(nlat, 1.0);
  DUCC0_NAMESPACE::cmav<std::size_t, 1> mstart(
      starts.data(), std::array<std::size_t, 1>{starts.size()});
  DUCC0_NAMESPACE::cmav<double, 1> rfactor(
      ringfactor.data(), std::array<std::size_t, 1>{ringfactor.size()});
  const std::array<std::size_t, 2> alm_shape{2, nalm};
  const std::array<std::size_t, 3> map_shape{2, nlat, nlon};
  const std::array<std::ptrdiff_t, 3> map_stride{static_cast<std::ptrdiff_t>(npix),
                                                 static_cast<std::ptrdiff_t>(nlon), 1};
  for (std::size_t f = 0; f < nfield; ++f)
  {
    const auto us = DUCC0_NAMESPACE::detail_pybind::copy_strides<T, false>(ugrid);
    const auto vs = DUCC0_NAMESPACE::detail_pybind::copy_strides<T, false>(vgrid);
    const auto* uptr = reinterpret_cast<const T*>(ugrid.data());
    const auto* vptr = reinterpret_cast<const T*>(vgrid.data());
    for (std::size_t j = 0; j < nlat; ++j)
      for (std::size_t i = 0; i < nlon; ++i)
      {
        const auto uidx = static_cast<std::ptrdiff_t>(j) * us[0] +
                          static_cast<std::ptrdiff_t>(i) * us[1] +
                          (batched ? static_cast<std::ptrdiff_t>(f) * us[2] : 0);
        const auto vidx = static_cast<std::ptrdiff_t>(j) * vs[0] +
                          static_cast<std::ptrdiff_t>(i) * vs[1] +
                          (batched ? static_cast<std::ptrdiff_t>(f) * vs[2] : 0);
        map_storage[j * nlon + i] = -vptr[vidx];
        map_storage[npix + j * nlon + i] = uptr[uidx];
      }
    if (lmax > 0)
    {
      DUCC0_NAMESPACE::vmav<std::complex<T>, 2> alms(alm_storage.data(), alm_shape);
      DUCC0_NAMESPACE::cmav<T, 3> maps(map_storage.data(), map_shape, map_stride);
      py::gil_scoped_release release;
      DUCC0_NAMESPACE::analysis_2d(alms, maps, 1, lmax, mstart, 1,
                                   gridtype == "regular" ? "CC" : "GL", 0.0, rfactor,
                                   transform_nthreads());
    }
    for (std::size_t m = 0; m <= lmax; ++m)
      for (std::size_t l = m; l <= lmax; ++l)
      {
        const auto idx = m * (2 * lmax - m + 3) / 2 + (l - m);
        const double degree = static_cast<double>(l);
        const double lambda = degree * (degree + 1.0);
        const double coefficient_scale =
            -std::sqrt(lambda) / (rsphere * sqrt_2pi) * coefficient_sign<T>(m);
        const double scale = l == 0 ? 0.0 : coefficient_scale;
        const auto outidx = batched ? idx * nfield + f : idx;
        vout.data()[outidx] = scale_coefficient(alm_storage[nalm + idx], scale);
        dout.data()[outidx] = scale_coefficient(alm_storage[idx], scale);
      }
  }
  return {vort, div};
}

template <typename T>
std::pair<NpArr, NpArr> synthesize_winds(const CNpArr& vrtspec, const CNpArr& divspec,
                                         std::size_t nlat, std::size_t nlon,
                                         const std::string& gridtype, double rsphere)
{
  const bool batched = vrtspec.ndim() == 2;
  const std::size_t nfield = batched ? vrtspec.shape(1) : 1;
  const std::size_t nalm = vrtspec.shape(0);
  const std::size_t lmax = truncation_from_count(nalm);
  NpArr u = batched ? new_array<T>({nlat, nlon, nfield}) : new_array<T>({nlat, nlon});
  NpArr v = batched ? new_array<T>({nlat, nlon, nfield}) : new_array<T>({nlat, nlon});
  auto uout = DUCC0_NAMESPACE::detail_pybind::to_vfmav<T>(u);
  auto vout = DUCC0_NAMESPACE::detail_pybind::to_vfmav<T>(v);
  const std::size_t npix = nlat * nlon;
  std::vector<T> map_storage(2 * npix);
  std::vector<std::complex<T>> alm_storage(2 * nalm);
  std::vector<std::size_t> starts = ducc_mstart(lmax);
  std::vector<double> ringfactor(nlat, 1.0);
  DUCC0_NAMESPACE::cmav<std::size_t, 1> mstart(
      starts.data(), std::array<std::size_t, 1>{starts.size()});
  DUCC0_NAMESPACE::cmav<double, 1> rfactor(
      ringfactor.data(), std::array<std::size_t, 1>{ringfactor.size()});
  const std::array<std::size_t, 2> alm_shape{2, nalm};
  const std::array<std::size_t, 3> map_shape{2, nlat, nlon};
  const std::array<std::ptrdiff_t, 3> map_stride{static_cast<std::ptrdiff_t>(npix),
                                                 static_cast<std::ptrdiff_t>(nlon), 1};
  for (std::size_t f = 0; f < nfield; ++f)
  {
    auto vort = as_ducc_alms<T>(vrtspec, nalm, f);
    auto div = as_ducc_alms<T>(divspec, nalm, f);
    std::fill(alm_storage.begin(), alm_storage.end(), std::complex<T>(0, 0));
    for (std::size_t m = 0; m <= lmax; ++m)
      for (std::size_t l = m; l <= lmax; ++l)
      {
        const auto idx = m * (2 * lmax - m + 3) / 2 + (l - m);
        if (l == 0)
          continue;
        const double degree = static_cast<double>(l);
        const double scale = -rsphere * sqrt_2pi * coefficient_sign<T>(m) /
                             std::sqrt(degree * (degree + 1.0));
        alm_storage[idx] = scale_coefficient(div(0, idx), scale);
        alm_storage[nalm + idx] = scale_coefficient(vort(0, idx), scale);
      }
    std::fill(map_storage.begin(), map_storage.end(), T(0));
    if (lmax > 0)
    {
      DUCC0_NAMESPACE::cmav<std::complex<T>, 2> alms(alm_storage.data(), alm_shape);
      DUCC0_NAMESPACE::vmav<T, 3> maps(map_storage.data(), map_shape, map_stride);
      py::gil_scoped_release release;
      DUCC0_NAMESPACE::synthesis_2d(alms, maps, 1, lmax, mstart, 1,
                                    gridtype == "regular" ? "CC" : "GL", 0.0, rfactor,
                                    transform_nthreads(), DUCC0_NAMESPACE::STANDARD);
    }
    for (std::size_t j = 0; j < nlat; ++j)
      for (std::size_t i = 0; i < nlon; ++i)
      {
        const auto pyidx = batched ? j * nlon * nfield + i * nfield + f : j * nlon + i;
        const auto pix = j * nlon + i;
        uout.data()[pyidx] = map_storage[npix + pix];
        vout.data()[pyidx] = -map_storage[pix];
      }
  }
  return {u, v};
}

template <typename T>
NpArr apply_inverse_laplacian(const CNpArr& spec, double rsphere)
{
  const auto nalm = spec.shape(0);
  const auto lmax = truncation_from_count(nalm);
  NpArr result = spec.ndim() == 1 ? new_array<std::complex<T>>({nalm})
                                  : new_array<std::complex<T>>({nalm, spec.shape(1)});
  auto input = DUCC0_NAMESPACE::detail_pybind::to_cfmav<std::complex<T>>(spec);
  auto output = DUCC0_NAMESPACE::detail_pybind::to_vfmav<std::complex<T>>(result);
  const auto nfield = spec.ndim() == 1 ? 1 : spec.shape(1);
  for (std::size_t f = 0; f < nfield; ++f)
    for (std::size_t m = 0; m <= lmax; ++m)
      for (std::size_t l = m; l <= lmax; ++l)
      {
        const auto idx = m * (2 * lmax - m + 3) / 2 + (l - m);
        const auto from = spec.ndim() == 1
                              ? static_cast<std::ptrdiff_t>(idx)
                              : static_cast<std::ptrdiff_t>(idx) * input.stride(0) +
                                    static_cast<std::ptrdiff_t>(f) * input.stride(1);
        const auto to = spec.ndim() == 1
                            ? static_cast<std::ptrdiff_t>(idx)
                            : static_cast<std::ptrdiff_t>(idx) * output.stride(0) +
                                  static_cast<std::ptrdiff_t>(f) * output.stride(1);
        const double degree = static_cast<double>(l);
        const double factor =
            l == 0 ? 0.0 : -rsphere * rsphere / (degree * (degree + 1.0));
        output.data()[to] = scale_coefficient(input.data()[from], factor);
      }
  return result;
}

template <typename T>
NpArr apply_smoothing(const CNpArr& spec, const CNpArr& smooth)
{
  const auto nalm = spec.shape(0);
  const auto lmax = truncation_from_count(nalm);
  if (smooth.ndim() != 1 || smooth.shape(0) < lmax + 1)
    throw py::value_error("smooth must be rank one and cover every coefficient degree");
  validate_array_dtype<T>(smooth,
                          std::is_same_v<T, float> ? "numpy.float32" : "numpy.float64");
  const bool batched = spec.ndim() == 2;
  const auto nfield = batched ? spec.shape(1) : 1;
  NpArr result = batched ? new_array<std::complex<T>>({nalm, nfield})
                         : new_array<std::complex<T>>({nalm});
  auto input = DUCC0_NAMESPACE::detail_pybind::to_cfmav<std::complex<T>>(spec);
  auto factors = DUCC0_NAMESPACE::detail_pybind::to_cfmav<T>(smooth);
  auto output = DUCC0_NAMESPACE::detail_pybind::to_vfmav<std::complex<T>>(result);
  for (std::size_t f = 0; f < nfield; ++f)
    for (std::size_t m = 0; m <= lmax; ++m)
      for (std::size_t l = m; l <= lmax; ++l)
      {
        const auto idx = m * (2 * lmax - m + 3) / 2 + (l - m);
        const auto from = batched ? static_cast<std::ptrdiff_t>(idx) * input.stride(0) +
                                        static_cast<std::ptrdiff_t>(f) * input.stride(1)
                                  : static_cast<std::ptrdiff_t>(idx);
        const auto to = batched ? static_cast<std::ptrdiff_t>(idx) * output.stride(0) +
                                      static_cast<std::ptrdiff_t>(f) * output.stride(1)
                                : static_cast<std::ptrdiff_t>(idx);
        output.data()[to] = input.data()[from] * factors(l);
      }
  return result;
}

template <typename To, typename From>
NpArr convert_spectral_precision(const CNpArr& input)
{
  const std::size_t nalm = input.shape(0);
  const bool batched = input.ndim() == 2;
  const std::size_t nfield = batched ? input.shape(1) : 1;
  NpArr result = batched ? new_array<std::complex<To>>({nalm, nfield})
                         : new_array<std::complex<To>>({nalm});
  const auto in_stride =
      DUCC0_NAMESPACE::detail_pybind::copy_strides<std::complex<From>, false>(input);
  auto output = DUCC0_NAMESPACE::detail_pybind::to_vfmav<std::complex<To>>(result);
  const auto* source = reinterpret_cast<const std::complex<From>*>(input.data());
  for (std::size_t l = 0; l < nalm; ++l)
    for (std::size_t f = 0; f < nfield; ++f)
    {
      const auto from = batched ? static_cast<std::ptrdiff_t>(l) * in_stride[0] +
                                      static_cast<std::ptrdiff_t>(f) * in_stride[1]
                                : static_cast<std::ptrdiff_t>(l) * in_stride[0];
      const auto to = batched ? static_cast<std::ptrdiff_t>(l) * output.stride(0) +
                                    static_cast<std::ptrdiff_t>(f) * output.stride(1)
                              : static_cast<std::ptrdiff_t>(l);
      output.data()[to] = static_cast<std::complex<To>>(source[from]);
    }
  return result;
}

}  // namespace

NpArr grid_to_spectral_f32(const CNpArr& grid, std::size_t nlat, std::size_t nlon,
                           const std::string& gridtype, std::size_t lmax)
{
  return grid_to_spectral<float>(grid, nlat, nlon, gridtype, lmax);
}

NpArr grid_to_spectral_f64(const CNpArr& grid, std::size_t nlat, std::size_t nlon,
                           const std::string& gridtype, std::size_t lmax)
{
  return grid_to_spectral<double>(grid, nlat, nlon, gridtype, lmax);
}

NpArr spectral_to_grid_f32(const CNpArr& spec, std::size_t nlat, std::size_t nlon,
                           const std::string& gridtype, std::size_t lmax)
{
  return spectral_to_grid<float>(spec, nlat, nlon, gridtype, lmax);
}

NpArr spectral_to_grid_f64(const CNpArr& spec, std::size_t nlat, std::size_t nlon,
                           const std::string& gridtype, std::size_t lmax)
{
  return spectral_to_grid<double>(spec, nlat, nlon, gridtype, lmax);
}

std::pair<NpArr, NpArr> spectral_gradient_f32(const CNpArr& spec, std::size_t nlat,
                                              std::size_t nlon,
                                              const std::string& gridtype,
                                              double rsphere)
{
  return spectral_gradient<float>(spec, nlat, nlon, gridtype, rsphere);
}

std::pair<NpArr, NpArr> spectral_gradient_f64(const CNpArr& spec, std::size_t nlat,
                                              std::size_t nlon,
                                              const std::string& gridtype,
                                              double rsphere)
{
  return spectral_gradient<double>(spec, nlat, nlon, gridtype, rsphere);
}

std::pair<NpArr, NpArr> analyze_winds_f32(const CNpArr& ugrid, const CNpArr& vgrid,
                                          std::size_t nlat, std::size_t nlon,
                                          const std::string& gridtype, double rsphere,
                                          std::size_t lmax)
{
  return analyze_winds<float>(ugrid, vgrid, nlat, nlon, gridtype, rsphere, lmax);
}

std::pair<NpArr, NpArr> analyze_winds_f64(const CNpArr& ugrid, const CNpArr& vgrid,
                                          std::size_t nlat, std::size_t nlon,
                                          const std::string& gridtype, double rsphere,
                                          std::size_t lmax)
{
  return analyze_winds<double>(ugrid, vgrid, nlat, nlon, gridtype, rsphere, lmax);
}

std::pair<NpArr, NpArr> synthesize_winds_f32(const CNpArr& vrtspec,
                                             const CNpArr& divspec, std::size_t nlat,
                                             std::size_t nlon,
                                             const std::string& gridtype,
                                             double rsphere)
{
  return synthesize_winds<float>(vrtspec, divspec, nlat, nlon, gridtype, rsphere);
}

std::pair<NpArr, NpArr> synthesize_winds_f64(const CNpArr& vrtspec,
                                             const CNpArr& divspec, std::size_t nlat,
                                             std::size_t nlon,
                                             const std::string& gridtype,
                                             double rsphere)
{
  return synthesize_winds<double>(vrtspec, divspec, nlat, nlon, gridtype, rsphere);
}

NpArr apply_inverse_laplacian_f32(const CNpArr& spec, double rsphere)
{
  return apply_inverse_laplacian<float>(spec, rsphere);
}

NpArr apply_inverse_laplacian_f64(const CNpArr& spec, double rsphere)
{
  return apply_inverse_laplacian<double>(spec, rsphere);
}

NpArr apply_smoothing_f32(const CNpArr& spec, const CNpArr& smooth)
{
  return apply_smoothing<float>(spec, smooth);
}

NpArr apply_smoothing_f64(const CNpArr& spec, const CNpArr& smooth)
{
  return apply_smoothing<double>(spec, smooth);
}

NpArr convert_spectral_precision_f32_from_f64(const CNpArr& input)
{
  return convert_spectral_precision<float, double>(input);
}

NpArr convert_spectral_precision_f64_from_f32(const CNpArr& input)
{
  return convert_spectral_precision<double, float>(input);
}

}  // namespace SPHARM_CPP_NAMESPACE
