// SPDX-FileCopyrightText: 2026 Albert M. W. Yau
// SPDX-License-Identifier: GPL-2.0-or-later
#include <nanobind/nanobind.h>

#include <algorithm>
#include <array>
#include <cstdint>

#include "ducc0/sht/sht.h"
#include "internal.hpp"

namespace SPHARM_CPP_NAMESPACE
{
namespace
{

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

struct TransformSetup
{
  std::vector<std::size_t> starts;
  std::vector<double> ringfactor;

  TransformSetup(std::size_t lmax, std::size_t nlat)
      : starts(ducc_mstart(lmax)), ringfactor(nlat, 1.0)
  {
  }

  DUCC0_NAMESPACE::cmav<std::size_t, 1> mstart() const
  {
    return {starts.data(), std::array<std::size_t, 1>{starts.size()}};
  }

  DUCC0_NAMESPACE::cmav<double, 1> rfactor() const
  {
    return {ringfactor.data(), std::array<std::size_t, 1>{ringfactor.size()}};
  }
};

template <typename T, std::size_t N>
NpArr make_component_view(T* data, const std::array<std::size_t, N>& shape,
                          const std::array<std::int64_t, N>& strides, py::handle owner,
                          std::size_t ndim)
{
  return NpArr(DUCC0_NAMESPACE::detail_pybind::NpArrT<T>(data, ndim, shape.data(),
                                                         owner, strides.data()));
}

struct FieldExecutionPlan
{
  std::size_t outer_threads;
  std::size_t inner_threads;
};

FieldExecutionPlan field_execution_plan(std::size_t nfield, std::size_t nlat,
                                        std::size_t nlon)
{
  const auto capacity = configured_nthreads();
  constexpr std::size_t min_parallel_fields = 4;
  constexpr std::size_t min_pixels = 32768;  // Medium Gaussian benchmark grid.
  const bool enough_pixels =
      nlon != 0 && nlat >= min_pixels / nlon + (min_pixels % nlon != 0);
  if (nfield >= min_parallel_fields && capacity > 1 && enough_pixels)
    return {std::min(nfield, capacity), 1};
  return {1, capacity};
}

template <typename Func>
void for_each_field_range(std::size_t nfield, const FieldExecutionPlan& plan,
                          Func&& func)
{
  if (plan.outer_threads == 1)
  {
    func(0, 0, nfield);
    return;
  }
  DUCC0_NAMESPACE::execParallel(0, nfield, plan.outer_threads,
                                std::forward<Func>(func));
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
void convert_ducc_to_py(DUCC0_NAMESPACE::vmav<std::complex<T>, 2>& alms,
                        std::size_t lmax, std::size_t nlon)
{
  for (std::size_t m = 0; m <= lmax; ++m)
    for (std::size_t l = m; l <= lmax; ++l)
    {
      const auto idx = packed_coefficient_index(lmax, m, l);
      const T nyquist_factor = (nlon % 2 == 0 && m == nlon / 2) ? T(0.5) : T(1);
      alms(0, idx) *= static_cast<T>(coefficient_sign(m) * nyquist_factor / sqrt_2pi);
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
  const TransformSetup setup(lmax, nlat);
  const auto mstart = setup.mstart();
  const auto rfactor = setup.rfactor();
  std::vector<DUCC0_NAMESPACE::vmav<std::complex<T>, 2>> alms;
  std::vector<DUCC0_NAMESPACE::cmav<T, 3>> maps;
  alms.reserve(nfield);
  maps.reserve(nfield);
  for (std::size_t f = 0; f < nfield; ++f)
  {
    alms.push_back(as_ducc_alms<T>(result, nalm, f, grid.ndim() == 3));
    maps.push_back(as_ducc_map<T>(grid, nlat, nlon, f));
  }
  const auto plan = field_execution_plan(nfield, nlat, nlon);
  {
    py::gil_scoped_release release;
    for_each_field_range(nfield, plan,
                         [&](std::size_t, std::size_t lo, std::size_t hi)
                         {
                           for (std::size_t f = lo; f < hi; ++f)
                           {
                             DUCC0_NAMESPACE::analysis_2d(
                                 alms[f], maps[f], 0, lmax, mstart, 1,
                                 gridtype == "regular" ? "CC" : "GL", 0.0, rfactor,
                                 plan.inner_threads);
                             convert_ducc_to_py<T>(alms[f], lmax, nlon);
                           }
                         });
  }
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
  std::vector<DUCC0_NAMESPACE::cmav<std::complex<T>, 2>> sources;
  std::vector<DUCC0_NAMESPACE::vmav<T, 3>> maps;
  sources.reserve(nfield);
  maps.reserve(nfield);
  for (std::size_t f = 0; f < nfield; ++f)
  {
    sources.push_back(as_ducc_alms<T>(spec, nalm, f));
    maps.push_back(as_ducc_map<T>(result, nlat, nlon, f, batched));
  }
  const TransformSetup setup(lmax, nlat);
  const auto mstart = setup.mstart();
  const auto rfactor = setup.rfactor();
  const auto plan = field_execution_plan(nfield, nlat, nlon);
  const std::array<std::size_t, 2> alm_shape{1, nalm};
  std::vector<std::vector<std::complex<T>>> workspace(
      plan.outer_threads, std::vector<std::complex<T>>(nalm));
  {
    py::gil_scoped_release release;
    for_each_field_range(
        nfield, plan,
        [&](std::size_t thread, std::size_t lo, std::size_t hi)
        {
          auto& storage = workspace[thread];
          DUCC0_NAMESPACE::cmav<std::complex<T>, 2> alms(storage.data(), alm_shape);
          for (std::size_t f = lo; f < hi; ++f)
          {
            for (std::size_t m = 0; m <= lmax; ++m)
              for (std::size_t l = m; l <= lmax; ++l)
              {
                const auto idx = packed_coefficient_index(lmax, m, l);
                storage[idx] =
                    sources[f](0, idx) * static_cast<T>(sqrt_2pi * coefficient_sign(m));
              }
            DUCC0_NAMESPACE::synthesis_2d(
                alms, maps[f], 0, lmax, mstart, 1, gridtype == "regular" ? "CC" : "GL",
                0.0, rfactor, plan.inner_threads, DUCC0_NAMESPACE::STANDARD);
          }
        });
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
  const std::size_t npix = nlat * nlon;
  const std::size_t output_plane = npix * nfield;
  NpArr output_storage = new_array<T>({2, nlat, nlon, nfield});
  auto* output_data = reinterpret_cast<T*>(output_storage.data());
  const auto owner = output_storage.cast();
  const std::array<std::size_t, 3> output_shape{nlat, nlon, nfield};
  const std::array<std::int64_t, 3> output_strides{
      static_cast<std::int64_t>(nlon * nfield), static_cast<std::int64_t>(nfield), 1};
  NpArr north = make_component_view(output_data, output_shape, output_strides, owner,
                                    batched ? 3 : 2);
  NpArr east = make_component_view(output_data + output_plane, output_shape,
                                   output_strides, owner, batched ? 3 : 2);
  if (lmax == 0)
  {
    std::fill_n(output_data, 2 * output_plane, T(0));
    return {east, north};
  }
  std::vector<DUCC0_NAMESPACE::cmav<std::complex<T>, 2>> sources;
  sources.reserve(nfield);
  for (std::size_t f = 0; f < nfield; ++f)
    sources.push_back(as_ducc_alms<T>(spec, nalm, f));
  const TransformSetup setup(lmax, nlat);
  const auto mstart = setup.mstart();
  const auto rfactor = setup.rfactor();
  const auto plan = field_execution_plan(nfield, nlat, nlon);
  const std::array<std::size_t, 2> alm_shape{1, nalm};
  const std::array<std::size_t, 3> map_shape{2, nlat, nlon};
  const std::array<std::ptrdiff_t, 3> map_stride{
      static_cast<std::ptrdiff_t>(output_plane),
      static_cast<std::ptrdiff_t>(nlon * nfield), static_cast<std::ptrdiff_t>(nfield)};
  std::vector<std::vector<std::complex<T>>> workspace(
      plan.outer_threads, std::vector<std::complex<T>>(nalm));
  {
    py::gil_scoped_release release;
    for_each_field_range(
        nfield, plan,
        [&](std::size_t thread, std::size_t lo, std::size_t hi)
        {
          auto& scratch = workspace[thread];
          DUCC0_NAMESPACE::cmav<std::complex<T>, 2> alm(scratch.data(), alm_shape);
          for (std::size_t f = lo; f < hi; ++f)
          {
            for (std::size_t m = 0; m <= lmax; ++m)
              for (std::size_t l = m; l <= lmax; ++l)
              {
                const auto idx = packed_coefficient_index(lmax, m, l);
                scratch[idx] =
                    sources[f](0, idx) * static_cast<T>(sqrt_2pi * coefficient_sign(m));
              }
            DUCC0_NAMESPACE::vmav<T, 3> map(output_data + f, map_shape, map_stride);
            DUCC0_NAMESPACE::synthesis_2d(
                alm, map, 1, lmax, mstart, 1, gridtype == "regular" ? "CC" : "GL", 0.0,
                rfactor, plan.inner_threads, DUCC0_NAMESPACE::DERIV1);
            for (std::size_t j = 0; j < nlat; ++j)
              for (std::size_t i = 0; i < nlon; ++i)
              {
                const auto output_index =
                    batched ? j * nlon * nfield + i * nfield + f : j * nlon + i;
                output_data[output_plane + output_index] = static_cast<T>(
                    static_cast<double>(output_data[output_plane + output_index]) /
                    rsphere);
                output_data[output_index] = static_cast<T>(
                    -static_cast<double>(output_data[output_index]) / rsphere);
              }
          }
        });
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
  const std::size_t output_plane = nalm * nfield;
  NpArr output_storage = new_array<std::complex<T>>({2, nalm, nfield});
  auto* output_data = reinterpret_cast<std::complex<T>*>(output_storage.data());
  const auto owner = output_storage.cast();
  const std::array<std::size_t, 2> output_shape{nalm, nfield};
  const std::array<std::int64_t, 2> output_strides{static_cast<std::int64_t>(nfield),
                                                   1};
  NpArr div = make_component_view(output_data, output_shape, output_strides, owner,
                                  batched ? 2 : 1);
  NpArr vort = make_component_view(output_data + output_plane, output_shape,
                                   output_strides, owner, batched ? 2 : 1);
  if (lmax == 0)
  {
    std::fill_n(output_data, 2 * output_plane, std::complex<T>(0, 0));
    return {vort, div};
  }
  const std::size_t npix = nlat * nlon;
  std::vector<DUCC0_NAMESPACE::cmav<T, 3>> umaps;
  std::vector<DUCC0_NAMESPACE::cmav<T, 3>> vmaps;
  umaps.reserve(nfield);
  vmaps.reserve(nfield);
  for (std::size_t f = 0; f < nfield; ++f)
  {
    umaps.push_back(as_ducc_map<T>(ugrid, nlat, nlon, f));
    vmaps.push_back(as_ducc_map<T>(vgrid, nlat, nlon, f));
  }
  const TransformSetup setup(lmax, nlat);
  const auto mstart = setup.mstart();
  const auto rfactor = setup.rfactor();
  const auto plan = field_execution_plan(nfield, nlat, nlon);
  const std::array<std::size_t, 2> alm_shape{2, nalm};
  const std::array<std::ptrdiff_t, 2> alm_stride{
      static_cast<std::ptrdiff_t>(output_plane), static_cast<std::ptrdiff_t>(nfield)};
  const std::array<std::size_t, 3> map_shape{2, nlat, nlon};
  const std::array<std::ptrdiff_t, 3> map_stride{static_cast<std::ptrdiff_t>(npix),
                                                 static_cast<std::ptrdiff_t>(nlon), 1};
  std::vector<std::vector<T>> workspace(plan.outer_threads, std::vector<T>(2 * npix));
  {
    py::gil_scoped_release release;
    for_each_field_range(
        nfield, plan,
        [&](std::size_t thread, std::size_t lo, std::size_t hi)
        {
          auto& scratch = workspace[thread];
          DUCC0_NAMESPACE::cmav<T, 3> maps(scratch.data(), map_shape, map_stride);
          for (std::size_t f = lo; f < hi; ++f)
          {
            for (std::size_t j = 0; j < nlat; ++j)
              for (std::size_t i = 0; i < nlon; ++i)
              {
                scratch[j * nlon + i] = -vmaps[f](0, j, i);
                scratch[npix + j * nlon + i] = umaps[f](0, j, i);
              }
            DUCC0_NAMESPACE::vmav<std::complex<T>, 2> alms(output_data + f, alm_shape,
                                                           alm_stride);
            DUCC0_NAMESPACE::analysis_2d(alms, maps, 1, lmax, mstart, 1,
                                         gridtype == "regular" ? "CC" : "GL", 0.0,
                                         rfactor, plan.inner_threads);
            for (std::size_t m = 0; m <= lmax; ++m)
              for (std::size_t l = m; l <= lmax; ++l)
              {
                const auto idx = packed_coefficient_index(lmax, m, l);
                if (l == 0)
                {
                  alms(1, idx) = std::complex<T>(0, 0);
                  alms(0, idx) = std::complex<T>(0, 0);
                  continue;
                }
                const double degree = static_cast<double>(l);
                const double lambda = degree * (degree + 1.0);
                const double scale =
                    -std::sqrt(lambda) / (rsphere * sqrt_2pi) * coefficient_sign(m);
                alms(1, idx) = scale_coefficient(alms(1, idx), scale);
                alms(0, idx) = scale_coefficient(alms(0, idx), scale);
              }
          }
        });
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
  const std::size_t output_plane = nlat * nlon * nfield;
  NpArr output_storage = new_array<T>({2, nlat, nlon, nfield});
  auto* output_data = reinterpret_cast<T*>(output_storage.data());
  const auto owner = output_storage.cast();
  const std::array<std::size_t, 3> output_shape{nlat, nlon, nfield};
  const std::array<std::int64_t, 3> output_strides{
      static_cast<std::int64_t>(nlon * nfield), static_cast<std::int64_t>(nfield), 1};
  NpArr v = make_component_view(output_data, output_shape, output_strides, owner,
                                batched ? 3 : 2);
  NpArr u = make_component_view(output_data + output_plane, output_shape,
                                output_strides, owner, batched ? 3 : 2);

  if (lmax == 0)
  {
    std::fill_n(output_data, 2 * output_plane, T(0));
    return {u, v};
  }

  std::vector<std::complex<T>> alm_storage(2 * nalm);
  const TransformSetup setup(lmax, nlat);
  const auto mstart = setup.mstart();
  const auto rfactor = setup.rfactor();
  const std::array<std::size_t, 2> alm_shape{2, nalm};
  const std::array<std::size_t, 3> map_shape{2, nlat, nlon};
  const std::array<std::ptrdiff_t, 3> map_stride{
      static_cast<std::ptrdiff_t>(output_plane),
      static_cast<std::ptrdiff_t>(nlon * nfield), static_cast<std::ptrdiff_t>(nfield)};
  for (std::size_t f = 0; f < nfield; ++f)
  {
    auto vort = as_ducc_alms<T>(vrtspec, nalm, f);
    auto div = as_ducc_alms<T>(divspec, nalm, f);
    for (std::size_t m = 0; m <= lmax; ++m)
      for (std::size_t l = m; l <= lmax; ++l)
      {
        const auto idx = packed_coefficient_index(lmax, m, l);
        if (l == 0)
          continue;
        const double degree = static_cast<double>(l);
        const double scale = -rsphere * sqrt_2pi * coefficient_sign(m) /
                             std::sqrt(degree * (degree + 1.0));
        alm_storage[idx] = scale_coefficient(div(0, idx), scale);
        alm_storage[nalm + idx] = scale_coefficient(vort(0, idx), scale);
      }
    DUCC0_NAMESPACE::cmav<std::complex<T>, 2> alms(alm_storage.data(), alm_shape);
    DUCC0_NAMESPACE::vmav<T, 3> maps(output_data + f, map_shape, map_stride);
    py::gil_scoped_release release;
    DUCC0_NAMESPACE::synthesis_2d(alms, maps, 1, lmax, mstart, 1,
                                  gridtype == "regular" ? "CC" : "GL", 0.0, rfactor,
                                  configured_nthreads(), DUCC0_NAMESPACE::STANDARD);
  }
  for (std::size_t idx = 0; idx < output_plane; ++idx)
    output_data[idx] = -output_data[idx];
  return {u, v};
}

template <typename T>
NpArr apply_inverse_laplacian(const CNpArr& spec, double rsphere)
{
  const std::size_t nalm = spec.shape(0);
  const std::size_t lmax = truncation_from_count(nalm);
  const bool batched = spec.ndim() == 2;
  const std::size_t nfield = batched ? spec.shape(1) : 1;
  NpArr result = batched ? new_array<std::complex<T>>({nalm, nfield})
                         : new_array<std::complex<T>>({nalm});
  for (std::size_t f = 0; f < nfield; ++f)
  {
    auto source = as_ducc_alms<T>(spec, nalm, f);
    auto dest = as_ducc_alms<T>(result, nalm, f, batched);
    for (std::size_t m = 0; m <= lmax; ++m)
      for (std::size_t l = m; l <= lmax; ++l)
      {
        const auto idx = packed_coefficient_index(lmax, m, l);
        const double degree = static_cast<double>(l);
        const double factor =
            l == 0 ? 0.0 : -rsphere * rsphere / (degree * (degree + 1.0));
        dest(0, idx) = scale_coefficient(source(0, idx), factor);
      }
  }
  return result;
}

template <typename T>
NpArr apply_smoothing(const CNpArr& spec, const CNpArr& smooth)
{
  const std::size_t nalm = spec.shape(0);
  const std::size_t lmax = truncation_from_count(nalm);
  if (smooth.ndim() != 1 || smooth.shape(0) < lmax + 1)
    throw py::value_error("smooth must be rank one and cover every coefficient degree");
  validate_array_dtype<T>(smooth,
                          std::is_same_v<T, float> ? "numpy.float32" : "numpy.float64");
  const bool batched = spec.ndim() == 2;
  const std::size_t nfield = batched ? spec.shape(1) : 1;
  NpArr result = batched ? new_array<std::complex<T>>({nalm, nfield})
                         : new_array<std::complex<T>>({nalm});
  const auto factors = DUCC0_NAMESPACE::detail_pybind::to_cfmav<T>(smooth);
  for (std::size_t f = 0; f < nfield; ++f)
  {
    auto source = as_ducc_alms<T>(spec, nalm, f);
    auto dest = as_ducc_alms<T>(result, nalm, f, batched);
    for (std::size_t m = 0; m <= lmax; ++m)
      for (std::size_t l = m; l <= lmax; ++l)
      {
        const auto idx = packed_coefficient_index(lmax, m, l);
        dest(0, idx) = source(0, idx) * factors(l);
      }
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
  for (std::size_t f = 0; f < nfield; ++f)
  {
    auto source = as_ducc_alms<From>(input, nalm, f);
    auto dest = as_ducc_alms<To>(result, nalm, f, batched);
    for (std::size_t l = 0; l < nalm; ++l)
      dest(0, l) = static_cast<std::complex<To>>(source(0, l));
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
