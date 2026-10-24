// SPDX-License-Identifier: GPL-2.0-or-later
#include <nanobind/nanobind.h>
#include <nanobind/stl/complex.h>
#include <nanobind/stl/vector.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <utility>
#include <vector>

#include "ducc0/infra/simd.h"
#include "ducc0/math/gl_integrator.h"
#include "internal.hpp"

namespace SPHARM_CPP_NAMESPACE
{
using namespace nanobind::literals;

namespace
{

constexpr double legendre_scale_high = 0x1p400;
constexpr double legendre_scale_low = 0x1p-400;
constexpr int legendre_scale_shift = 800;

// Keep tiny Pbar(m,m) seeds representable and carry their binary exponent
// separately, following DUCC's power-of-two recurrence scaling strategy.
void normalize_legendre_scale(double& value, int& exponent)
{
  while (std::abs(value) > legendre_scale_high)
  {
    value = std::scalbn(value, -legendre_scale_shift);
    exponent += legendre_scale_shift;
  }
  while (value != 0. && std::abs(value) < legendre_scale_low)
  {
    value = std::scalbn(value, legendre_scale_shift);
    exponent -= legendre_scale_shift;
  }
}

void normalize_legendre_scale(double& previous, double& current, int& exponent)
{
  double largest = std::max(std::abs(previous), std::abs(current));
  while (largest > legendre_scale_high)
  {
    previous = std::scalbn(previous, -legendre_scale_shift);
    current = std::scalbn(current, -legendre_scale_shift);
    exponent += legendre_scale_shift;
    largest = std::max(std::abs(previous), std::abs(current));
  }
  while (largest != 0. && largest < legendre_scale_low)
  {
    previous = std::scalbn(previous, legendre_scale_shift);
    current = std::scalbn(current, legendre_scale_shift);
    exponent -= legendre_scale_shift;
    largest = std::max(std::abs(previous), std::abs(current));
  }
}

double restore_legendre_scale(double value, int exponent)
{
  return exponent == 0 ? value : std::scalbn(value, exponent);
}

template <typename Tv>
void normalize_legendre_scale(Tv& value, std::array<int, Tv::size()>& exponent)
{
  Tv factor(1.);
  for (std::size_t lane = 0; lane < Tv::size(); ++lane)
  {
    auto magnitude = std::abs(value[lane]);
    int shift = 0;
    while (magnitude > legendre_scale_high)
    {
      magnitude = std::scalbn(magnitude, -legendre_scale_shift);
      shift -= legendre_scale_shift;
    }
    while (magnitude != 0. && magnitude < legendre_scale_low)
    {
      magnitude = std::scalbn(magnitude, legendre_scale_shift);
      shift += legendre_scale_shift;
    }
    factor[lane] = std::scalbn(1., shift);
    exponent[lane] -= shift;
  }
  value *= factor;
}

template <typename Tv>
void normalize_legendre_scale(Tv& previous, Tv& current,
                              std::array<int, Tv::size()>& exponent)
{
  bool needs_rescale = false;
  for (std::size_t lane = 0; lane < Tv::size(); ++lane)
  {
    const auto magnitude = std::max(std::abs(previous[lane]), std::abs(current[lane]));
    if (magnitude > legendre_scale_high ||
        (magnitude != 0. && magnitude < legendre_scale_low))
    {
      needs_rescale = true;
      break;
    }
  }
  if (!needs_rescale)
    return;

  Tv factor(1.);
  for (std::size_t lane = 0; lane < Tv::size(); ++lane)
  {
    auto magnitude = std::max(std::abs(previous[lane]), std::abs(current[lane]));
    if (magnitude <= legendre_scale_high &&
        (magnitude == 0. || magnitude >= legendre_scale_low))
      continue;
    int shift = 0;
    while (magnitude > legendre_scale_high)
    {
      magnitude = std::scalbn(magnitude, -legendre_scale_shift);
      shift -= legendre_scale_shift;
    }
    while (magnitude != 0. && magnitude < legendre_scale_low)
    {
      magnitude = std::scalbn(magnitude, legendre_scale_shift);
      shift += legendre_scale_shift;
    }
    if (shift != 0)
    {
      factor[lane] = std::scalbn(1., shift);
      exponent[lane] -= shift;
    }
  }
  previous *= factor;
  current *= factor;
}

template <typename T>
NpArr legendre_values(double latitude, std::size_t lmax)
{
  NpArr result = new_array<T>({coefficient_count(lmax)});
  auto values = DUCC0_NAMESPACE::detail_pybind::to_vfmav<T>(result);
  if (latitude == 90.0 || latitude == -90.0)
  {
    // At a pole only m=0 survives, and P_l(+-1) is exactly (+-)^l.
    // Avoid the O(l) recurrence's accumulated roundoff in this known case.
    std::fill_n(values.data(), coefficient_count(lmax), T(0));
    const double sign = latitude < 0 ? -1.0 : 1.0;
    for (std::size_t l = 0; l <= lmax; ++l)
      values.data()[l] =
          static_cast<T>((sign < 0 && l % 2 != 0 ? -1.0 : 1.0) *
                         std::sqrt((2.0 * static_cast<double>(l) + 1.0) / 2.0));
    return result;
  }
  const double phi = latitude * pi / 180.0;
  const double x = std::sin(phi);
  const double s = std::cos(phi);
  double diagonal = 1.0 / std::sqrt(2.0);
  int diagonal_exponent = 0;
  for (std::size_t m = 0; m <= lmax; ++m)
  {
    const auto base = m * (2 * lmax - m + 3) / 2;
    values.data()[base] =
        static_cast<T>(restore_legendre_scale(diagonal, diagonal_exponent));
    if (m < lmax)
    {
      double previous = diagonal;
      double current = x * std::sqrt(2.0 * static_cast<double>(m) + 3) * diagonal;
      int exponent = diagonal_exponent;
      normalize_legendre_scale(previous, current, exponent);
      values.data()[base + 1] =
          static_cast<T>(restore_legendre_scale(current, exponent));
      for (std::size_t l = m + 2; l <= lmax; ++l)
      {
        const double ld = static_cast<double>(l), md = static_cast<double>(m);
        const double lp = ld + md, lm = ld - md;
        const double a = std::sqrt(((2 * ld - 1) * (2 * ld + 1)) / (lp * lm));
        const double b =
            std::sqrt(((2 * ld + 1) * (lp - 1) * (lm - 1)) / ((2 * ld - 3) * lp * lm));
        const double next = a * x * current - b * previous;
        previous = current;
        current = next;
        normalize_legendre_scale(previous, current, exponent);
        values.data()[base + (l - m)] =
            static_cast<T>(restore_legendre_scale(current, exponent));
      }
      const double md = static_cast<double>(m + 1);
      diagonal *= s * std::sqrt((2.0 * md + 1) / (2.0 * md));
      normalize_legendre_scale(diagonal, diagonal_exponent);
    }
  }
  return result;
}

template <typename T>
NpArr legendre_many_values(const std::vector<double>& latitudes, std::size_t lmax)
{
  const auto nalm = coefficient_count(lmax);
  if (nalm && latitudes.size() > std::numeric_limits<std::size_t>::max() / nalm)
    throw py::value_error("latitude batch is too large");
  NpArr result = new_array<T>({latitudes.size(), nalm});
  auto output = DUCC0_NAMESPACE::detail_pybind::to_vfmav<T>(result);
  std::vector<double> x(latitudes.size()), s(latitudes.size());
  for (std::size_t i = 0; i < latitudes.size(); ++i)
  {
    const auto lat = latitudes[i];
    if (!std::isfinite(lat) || lat < -90.0 || lat > 90.0)
      throw py::value_error("latitudes must be finite and between -90 and 90 degrees");
    const double phi = lat * pi / 180.0;
    x[i] = (lat == 90.0) ? 1.0 : ((lat == -90.0) ? -1.0 : std::sin(phi));
    s[i] = (lat == 90.0 || lat == -90.0) ? 0.0 : std::cos(phi);
  }

  using Tv = DUCC0_NAMESPACE::native_simd<double>;
  constexpr std::size_t vlen = Tv::size();
  std::size_t start = 0;
  {
    py::gil_scoped_release release;
    for (; start + vlen <= latitudes.size(); start += vlen)
    {
      const Tv xv = DUCC0_NAMESPACE::loadu<Tv>(x.data() + start);
      const Tv sv = DUCC0_NAMESPACE::loadu<Tv>(s.data() + start);
      Tv diagonal(1.0 / std::sqrt(2.0));
      std::array<int, vlen> diagonal_exponent{};
      for (std::size_t m = 0; m <= lmax; ++m)
      {
        const auto base = m * (2 * lmax - m + 3) / 2;
        const auto store =
            [&](std::size_t coeff, Tv value, const std::array<int, vlen>& exponent)
        {
          for (std::size_t lane = 0; lane < vlen; ++lane)
            output
                .data()[(start + lane) * output.stride(0) + coeff * output.stride(1)] =
                static_cast<T>(restore_legendre_scale(value[lane], exponent[lane]));
        };
        const Tv pmm = diagonal;
        store(base, pmm, diagonal_exponent);
        if (m < lmax)
        {
          Tv pm1 = xv * std::sqrt(2.0 * static_cast<double>(m) + 3.0) * pmm;
          Tv pm2 = pmm;
          auto exponent = diagonal_exponent;
          normalize_legendre_scale(pm2, pm1, exponent);
          store(base + 1, pm1, exponent);
          for (std::size_t l = m + 2; l <= lmax; ++l)
          {
            const double ld = static_cast<double>(l), md = static_cast<double>(m);
            const double lp = ld + md, lm = ld - md;
            const double a = std::sqrt(((2 * ld - 1) * (2 * ld + 1)) / (lp * lm));
            const double b = std::sqrt(((2 * ld + 1) * (lp - 1) * (lm - 1)) /
                                       ((2 * ld - 3) * lp * lm));
            const Tv pl = a * xv * pm1 - b * pm2;
            pm2 = pm1;
            pm1 = pl;
            normalize_legendre_scale(pm2, pm1, exponent);
            store(base + (l - m), pm1, exponent);
          }
          const double md = static_cast<double>(m + 1);
          diagonal = sv * std::sqrt((2.0 * md + 1) / (2.0 * md)) * diagonal;
          normalize_legendre_scale(diagonal, diagonal_exponent);
        }
      }
    }
  }
  for (; start < latitudes.size(); ++start)
  {
    auto row = legendre_values<T>(latitudes[start], lmax);
    auto source = DUCC0_NAMESPACE::detail_pybind::to_cfmav<T>(row);
    for (std::size_t coeff = 0; coeff < nalm; ++coeff)
      output.data()[start * output.stride(0) + coeff * output.stride(1)] =
          source(coeff);
  }
  for (std::size_t row = 0; row < latitudes.size(); ++row)
    if (latitudes[row] == 90.0 || latitudes[row] == -90.0)
    {
      for (std::size_t coeff = 0; coeff < nalm; ++coeff)
        output.data()[row * output.stride(0) + coeff * output.stride(1)] = T(0);
      const double sign = latitudes[row] < 0 ? -1.0 : 1.0;
      for (std::size_t l = 0; l <= lmax; ++l)
        output.data()[row * output.stride(0) + l * output.stride(1)] =
            static_cast<T>((sign < 0 && l % 2 != 0 ? -1.0 : 1.0) *
                           std::sqrt((2.0 * static_cast<double>(l) + 1.0) / 2.0));
    }
  return result;
}

template <typename T>
struct Vec3
{
  T x, y, z;
};

template <typename T>
Vec3<T> spherical_roundtrip(Vec3<T> p)
{
  const T r1 = p.x * p.x + p.y * p.y;
  T theta = T(0), phi = T(0);
  if (r1 == T(0))
  {
    if (p.z < T(0))
      theta = T(4) * std::atan(T(1));
  }
  else
  {
    const T radial = std::sqrt(r1);
    theta = std::atan2(radial, p.z);
    phi = std::atan2(p.y, p.x);
  }
  return {std::sin(theta) * std::cos(phi), std::sin(theta) * std::sin(phi),
          std::cos(theta)};
}

template <typename T>
std::pair<std::vector<T>, std::vector<T>> geodesic_points(std::size_t edge_points)
{
  // Construct all points in one native call. Preserve the ordered scalar
  // arithmetic of each selected precision; vectorizing these conversions
  // would add lane fallbacks without useful speed.
  if (edge_points == 0)
    throw py::value_error("m must be at least 1");
  if (edge_points == 1)
    return {{T(90), T(-90)}, {T(0), T(0)}};
  if (edge_points > (std::numeric_limits<std::size_t>::max() - 1) / 2)
    throw py::value_error("m is too large");
  const std::size_t n = edge_points - 1;
  if (n > 0 && n > (std::numeric_limits<std::size_t>::max() - 2) / 10 / n)
    throw py::value_error("m is too large");
  const std::size_t expected = 10 * n * n + 2;
  std::vector<T> latitudes(expected), longitudes(expected);
  const std::size_t idp = 2 * edge_points - 1, jdp = edge_points;
  if (idp > std::numeric_limits<std::size_t>::max() / jdp / 5)
    throw py::value_error("m is too large");
  std::vector<Vec3<T>> xyz(idp * jdp * 5);
  auto at = [&](std::size_t i, std::size_t j, std::size_t k) -> Vec3<T>&
  { return xyz[(i * jdp + j) * 5 + k]; };
  auto stoc = [](T theta, T phi)
  {
    return Vec3<T>{std::sin(theta) * std::cos(phi), std::sin(theta) * std::sin(phi),
                   std::cos(theta)};
  };
  const T pif = T(4) * std::atan(T(1));
  const T dphi = T(0.4) * pif, hdphi = T(0.5) * dphi, tdphi = T(3) * hdphi;
  const T beta = std::cos(dphi);
  const T theta1 = std::acos(beta / (T(1) - beta)), theta2 = pif - theta1;
  for (std::size_t k = 0; k < 5; ++k)
  {
    const T phi = static_cast<T>(k) * dphi;
    const Vec3<T> x1 = stoc(theta2, phi), x2 = stoc(pif, phi + hdphi),
                  x3 = stoc(theta2, phi + dphi), x4 = stoc(theta1, phi + hdphi),
                  x5 = stoc(theta1, phi + tdphi), x6 = stoc(T(0), phi + dphi);
    auto delta = [](Vec3<T> a, Vec3<T> b, T divisor)
    {
      return Vec3<T>{(b.x - a.x) / divisor, (b.y - a.y) / divisor,
                     (b.z - a.z) / divisor};
    };
    auto add = [](Vec3<T> a, Vec3<T> b, T scale)
    { return Vec3<T>{a.x + scale * b.x, a.y + scale * b.y, a.z + scale * b.z}; };
    auto dxi = delta(x1, x2, static_cast<T>(n));
    auto dyj = delta(x2, x3, static_cast<T>(n));
    for (std::size_t i = 0; i < edge_points; ++i)
    {
      const auto xs = add(x1, dxi, static_cast<T>(i));
      for (std::size_t j = 0; j <= i; ++j)
        at(j, i, k) = add(xs, dyj, static_cast<T>(j));
    }
    dxi = delta(x4, x3, static_cast<T>(n));
    auto dxj = delta(x1, x4, static_cast<T>(n));
    for (std::size_t j = 0; j < edge_points; ++j)
    {
      const auto xs = add(x1, dxj, static_cast<T>(j));
      for (std::size_t i = 0; i <= j; ++i)
        at(j, i, k) = add(xs, dxi, static_cast<T>(i));
    }
    dxi = delta(x4, x3, static_cast<T>(n));
    dxj = delta(x3, x5, static_cast<T>(n));
    for (std::size_t i = 0; i < edge_points; ++i)
    {
      const auto xs = add(x4, dxi, static_cast<T>(i));
      for (std::size_t j = 0; j <= i; ++j)
        at(j + edge_points - 1, i, k) = add(xs, dxj, static_cast<T>(j));
    }
    dxi = delta(x6, x5, static_cast<T>(n));
    dxj = delta(x4, x6, static_cast<T>(n));
    for (std::size_t j = 0; j < edge_points; ++j)
    {
      const auto xs = add(x4, dxj, static_cast<T>(j));
      for (std::size_t i = 0; i <= j; ++i)
        at(j + edge_points - 1, i, k) = add(xs, dxi, static_cast<T>(i));
    }
  }
  std::size_t out = 2;
  for (std::size_t i = 0; i < 2 * (edge_points - 1); ++i)
    for (std::size_t j = 0; j < edge_points - 1; ++j)
      for (std::size_t k = 0; k < 5; ++k)
      {
        const auto p = spherical_roundtrip(at(i, j, k));
        const T radial = std::sqrt(p.x * p.x + p.y * p.y);
        const T rad2dg = T(180) / pif;
        const T xtmp = (p.x != T(0) || p.y != T(0)) ? p.x : T(1);
        latitudes[out] = std::atan2(radial, p.z) * rad2dg - T(90);
        longitudes[out] = std::atan2(p.y, xtmp) * rad2dg + T(180);
        ++out;
      }
  latitudes[0] = T(90);
  latitudes[1] = T(-90);
  longitudes[0] = T(0);
  longitudes[1] = T(0);
  if (out != expected)
    throw std::runtime_error("icosahedral subdivision produced an invalid point count");
  return {latitudes, longitudes};
}

template <typename T>
double spectral_interpolate(double longitude, const CNpArr& spec,
                            const CNpArr& legfuncs)
{
  if (spec.ndim() != 1 || legfuncs.ndim() != 1)
    throw py::value_error("dataspec and legfuncs must be rank-one arrays");
  if (spec.shape(0) == 0 || legfuncs.shape(0) == 0)
    throw py::value_error("dataspec and legfuncs must not be empty");
  const auto ntrunc = truncation_from_count(spec.shape(0));
  if (truncation_from_count(legfuncs.shape(0)) != ntrunc)
    throw py::value_error(
        "dataspec and legfuncs imply inconsistent spectral truncations");
  const auto ss =
      DUCC0_NAMESPACE::detail_pybind::copy_strides<std::complex<T>, false>(spec);
  const auto ps = DUCC0_NAMESPACE::detail_pybind::copy_strides<T, false>(legfuncs);
  const auto* coeff = reinterpret_cast<const std::complex<T>*>(spec.data());
  const auto* pbar = reinterpret_cast<const T*>(legfuncs.data());
  const double reduced_degrees = std::remainder(longitude, 360.0);
  const double lon = reduced_degrees * (pi / 180.0);
  double result = 0;
  for (std::size_t m = 0; m <= ntrunc; ++m)
  {
    const double c = std::cos(static_cast<double>(m) * lon);
    const double s = std::sin(static_cast<double>(m) * lon);
    for (std::size_t l = m; l <= ntrunc; ++l)
    {
      const auto idx = m * (2 * ntrunc - m + 3) / 2 + (l - m);
      const auto co = coeff[static_cast<std::ptrdiff_t>(idx) * ss[0]];
      const double pl = pbar[static_cast<std::ptrdiff_t>(idx) * ps[0]];
      result += pl * (m == 0 ? static_cast<double>(co.real())
                             : 2.0 * (static_cast<double>(co.real()) * c -
                                      static_cast<double>(co.imag()) * s));
    }
  }
  return result;
}

template <typename T>
NpArr spectral_interpolate_many(const std::vector<double>& longitudes,
                                const CNpArr& spec, const CNpArr& legfuncs)
{
  const auto ntrunc = truncation_from_count(spec.shape(0));
  auto coefficients = DUCC0_NAMESPACE::detail_pybind::to_cfmav<std::complex<T>>(spec);
  auto legendre = DUCC0_NAMESPACE::detail_pybind::to_cfmav<T>(legfuncs);
  NpArr result = new_array<T>({longitudes.size()});
  auto output = DUCC0_NAMESPACE::detail_pybind::to_vfmav<T>(result);
  using Tv = DUCC0_NAMESPACE::native_simd<double>;
  constexpr std::size_t vlen = Tv::size();
  std::size_t point = 0;
  {
    py::gil_scoped_release release;
    for (; point + vlen <= longitudes.size(); point += vlen)
    {
      std::array<double, vlen> reduced_degrees{};
      for (std::size_t lane = 0; lane < vlen; ++lane)
        reduced_degrees[lane] = std::remainder(longitudes[point + lane], 360.0);
      const Tv lon = DUCC0_NAMESPACE::loadu<Tv>(reduced_degrees.data()) * (pi / 180.0);
      Tv cosine(1.0), sine(0.0);
      Tv cosine_step(1.0), sine_step(0.0);
      if (ntrunc > 0)
      {
        cosine_step = DUCC0_NAMESPACE::detail_simd::apply(
            lon, [](double value) { return std::cos(value); });
        sine_step = DUCC0_NAMESPACE::detail_simd::apply(
            lon, [](double value) { return std::sin(value); });
      }
      Tv sum(0.0);
      for (std::size_t m = 0; m <= ntrunc; ++m)
      {
        if (m > 0)
        {
          if (m % 32 == 0)
          {
            const Tv phase = lon * static_cast<double>(m);
            cosine = DUCC0_NAMESPACE::detail_simd::apply(
                phase, [](double value) { return std::cos(value); });
            sine = DUCC0_NAMESPACE::detail_simd::apply(
                phase, [](double value) { return std::sin(value); });
          }
          else
          {
            const Tv next_cosine = cosine * cosine_step - sine * sine_step;
            const Tv next_sine = sine * cosine_step + cosine * sine_step;
            cosine = next_cosine;
            sine = next_sine;
          }
        }
        Tv sum_re(0.0), sum_im(0.0);
        for (std::size_t l = m; l <= ntrunc; ++l)
        {
          const auto idx = m * (2 * ntrunc - m + 3) / 2 + (l - m);
          const auto coeff = coefficients(idx);
          std::array<double, vlen> pvalues{};
          for (std::size_t lane = 0; lane < vlen; ++lane)
            pvalues[lane] = static_cast<double>(legendre(point + lane, idx));
          const Tv p = DUCC0_NAMESPACE::loadu<Tv>(pvalues.data());
          sum_re += p * static_cast<double>(coeff.real());
          // GCC 13.3 v4 lowered the former scalar m>0 guard to k1=0x01 on this
          // vector multiply, so only lane 0 accumulated the imaginary term.
          // Use a zero factor at m=0 and update all lanes unconditionally.
          const double imag_coefficient =
              (m == 0) ? 0.0 : static_cast<double>(coeff.imag());
          sum_im += p * imag_coefficient;
        }
        if (m == 0)
          sum += sum_re;
        else
          sum += 2.0 * (sum_re * cosine - sum_im * sine);
      }
      for (std::size_t lane = 0; lane < vlen; ++lane)
        output(point + lane) = static_cast<T>(sum[lane]);
    }
    for (; point < longitudes.size(); ++point)
    {
      double sum = 0.0;
      const double lon = std::remainder(longitudes[point], 360.0) * pi / 180.0;
      for (std::size_t m = 0; m <= ntrunc; ++m)
      {
        const double cosine = std::cos(static_cast<double>(m) * lon);
        const double sine = std::sin(static_cast<double>(m) * lon);
        double sum_re = 0.0, sum_im = 0.0;
        for (std::size_t l = m; l <= ntrunc; ++l)
        {
          const auto idx = m * (2 * ntrunc - m + 3) / 2 + (l - m);
          const auto coeff = coefficients(idx);
          const double p = static_cast<double>(legendre(point, idx));
          sum_re += p * static_cast<double>(coeff.real());
          if (m > 0)
            sum_im += p * static_cast<double>(coeff.imag());
        }
        if (m == 0)
          sum += sum_re;
        else
          sum += 2.0 * (sum_re * cosine - sum_im * sine);
      }
      output(point) = static_cast<T>(sum);
    }
  }
  return result;
}

}  // namespace

void bind_numerics(py::module_& m)
{
  m.def(
      "gaussian_lats_wts",
      [](std::size_t nlat)
      {
        if (nlat == 0)
          throw py::value_error("nlat must be positive");
        const DUCC0_NAMESPACE::GL_Integrator integrator(nlat);
        const auto theta = integrator.thetas();
        const auto weights = integrator.weights();
        NpArr lats = new_array<double>({nlat});
        NpArr wts = new_array<double>({nlat});
        auto la = DUCC0_NAMESPACE::detail_pybind::to_vfmav<double>(lats);
        auto wa = DUCC0_NAMESPACE::detail_pybind::to_vfmav<double>(wts);
        // DUCC stores theta south-to-north; the public API uses north-to-south.
        for (std::size_t i = 0; i < (nlat + 1) / 2; ++i)
        {
          const auto opposite = nlat - 1 - i;
          const double latitude = 90.0 - theta[opposite] * (180.0 / pi);
          la.data()[i] = latitude;
          la.data()[opposite] = -latitude;
          wa.data()[i] = weights[i];
          wa.data()[opposite] = weights[i];
        }
        return py::make_tuple(lats, wts);
      },
      "nlat"_a,
      py::sig(
          R"(def gaussian_lats_wts(nlat: int) -> tuple[numpy.typing.NDArray[numpy.float64], numpy.typing.NDArray[numpy.float64]])"),
      R"doc(Return Gauss-Legendre latitudes and quadrature weights.

Parameters
----------
nlat : int
    Positive number of Gaussian latitude rings.

Returns
-------
(latitudes, weights) : tuple of ndarray
    Float64 arrays of length nlat. Latitude is in degrees and
    decreases from north to south. Weights integrate over sine
    latitude and sum to 2 up to floating-point roundoff.)doc");

  m.def(
      "getspecindx",
      [](std::int64_t ntrunc)
      {
        if (ntrunc < 0)
          throw py::value_error("ntrunc must be nonnegative");
        const auto lmax = static_cast<std::size_t>(ntrunc);
        const auto n = coefficient_count(lmax);
        static_assert(sizeof(std::intptr_t) == sizeof(void*),
                      "getspecindx requires pointer-width NumPy indices");
        NpArr indxm = new_array<std::intptr_t>({n});
        NpArr indxn = new_array<std::intptr_t>({n});
        auto mm = DUCC0_NAMESPACE::detail_pybind::to_vfmav<std::intptr_t>(indxm);
        auto nn = DUCC0_NAMESPACE::detail_pybind::to_vfmav<std::intptr_t>(indxn);
        std::size_t k = 0;
        for (std::size_t m = 0; m <= lmax; ++m)
          for (std::size_t l = m; l <= lmax; ++l, ++k)
          {
            mm.data()[k] = static_cast<std::intptr_t>(m);
            nn.data()[k] = static_cast<std::intptr_t>(l);
          }
        return py::make_tuple(indxm, indxn);
      },
      "ntrunc"_a,
      py::sig(
          R"(def getspecindx(ntrunc: int) -> tuple[numpy.typing.NDArray[numpy.intp], numpy.typing.NDArray[numpy.intp]])"),
      R"doc(Return longitudinal orders and degrees of packed spectral coefficients.

Parameters
----------
ntrunc : int
    Nonnegative maximum spherical harmonic degree.

Returns
-------
(order, degree) : tuple of ndarray
    Integer arrays of length (ntrunc + 1) * (ntrunc + 2) / 2.
    Order m is the outer index and degree l runs from m to ntrunc.)doc");

  m.def(
      "legendre",
      [](double latitude, std::int64_t ntrunc, py::object dtype)
      {
        if (ntrunc < 0)
          throw py::value_error("ntrunc must be nonnegative");
        if (!std::isfinite(latitude) || latitude < -90.0 || latitude > 90.0)
          throw py::value_error("lat must be finite and between -90 and 90 degrees");
        if (dtype.is_none())
          dtype = py::module_::import_("numpy").attr("float64");
        auto name = py::cast<std::string>(
            py::module_::import_("numpy").attr("dtype")(dtype).attr("name"));
        if (name == "float32")
          return legendre_values<float>(latitude, static_cast<std::size_t>(ntrunc));
        if (name == "float64")
          return legendre_values<double>(latitude, static_cast<std::size_t>(ntrunc));
        throw py::type_error("dtype must be numpy.float32 or numpy.float64");
      },
      "lat"_a, "ntrunc"_a, py::kw_only(), "dtype"_a = py::none(),
      py::sig(
          R"(def legendre(lat: float, ntrunc: int, *, dtype: object | None = None) -> numpy.typing.NDArray[numpy.float32] | numpy.typing.NDArray[numpy.float64])"),
      R"doc(Evaluate normalized associated Legendre functions at one latitude.

Parameters
----------
lat : float
    Latitude in degrees, between -90 and 90.
ntrunc : int
    Nonnegative maximum spherical harmonic degree.
dtype : numpy dtype, optional
    Output float32 or float64 (default).

Returns
-------
ndarray
    Real values in triangular coefficient order, of length
    (ntrunc + 1) * (ntrunc + 2) / 2.

Notes
-----
The normalization excludes the Condon-Shortley phase. In particular,
the degree-zero value is 1 / sqrt(2).)doc");

  m.def(
      "legendre_many",
      [](const CNpArr& latitude_array, std::int64_t ntrunc, py::object dtype)
      {
        if (ntrunc < 0)
          throw py::value_error("ntrunc must be nonnegative");
        auto latitudes = real_vector(latitude_array, "latitudes");
        if (dtype.is_none())
          dtype = py::module_::import_("numpy").attr("float64");
        auto name = py::cast<std::string>(
            py::module_::import_("numpy").attr("dtype")(dtype).attr("name"));
        if (name == "float32")
          return legendre_many_values<float>(latitudes,
                                             static_cast<std::size_t>(ntrunc));
        if (name == "float64")
          return legendre_many_values<double>(latitudes,
                                              static_cast<std::size_t>(ntrunc));
        throw py::type_error("dtype must be numpy.float32 or numpy.float64");
      },
      "latitudes"_a, "ntrunc"_a, py::kw_only(), "dtype"_a = py::none(),
      py::sig(
          R"(def legendre_many(latitudes: numpy.typing.NDArray[numpy.float32] | numpy.typing.NDArray[numpy.float64], ntrunc: int, *, dtype: object | None = None) -> numpy.typing.NDArray[numpy.float32] | numpy.typing.NDArray[numpy.float64])"),
      R"doc(Evaluate normalized associated Legendre functions at several latitudes.

Parameters
----------
latitudes : ndarray
    One-dimensional real latitude array in degrees, between -90 and 90.
ntrunc : int
    Nonnegative maximum spherical harmonic degree.
dtype : numpy dtype, optional
    Output float32 or float64 (default).

Returns
-------
ndarray
    Real array of shape (npoints, ncoeff), using the same
    normalization and packing as legendre.)doc");

  m.def(
      "getgeodesicpts",
      [](std::int64_t m, py::object dtype)
      {
        if (m < 1)
          throw py::value_error("m must be at least 1");
        if (dtype.is_none())
          dtype = py::module_::import_("numpy").attr("float64");
        const auto name = py::cast<std::string>(
            py::module_::import_("numpy").attr("dtype")(dtype).attr("name"));
        if (name == "float32")
        {
          auto result = geodesic_points<double>(static_cast<std::size_t>(m));
          NpArr latitudes = new_array<float>({result.first.size()});
          NpArr longitudes = new_array<float>({result.second.size()});
          auto lat = DUCC0_NAMESPACE::detail_pybind::to_vfmav<float>(latitudes);
          auto lon = DUCC0_NAMESPACE::detail_pybind::to_vfmav<float>(longitudes);
          for (std::size_t i = 0; i < result.first.size(); ++i)
          {
            lat.data()[i] = static_cast<float>(result.first[i]);
            lon.data()[i] = static_cast<float>(result.second[i]);
          }
          return py::make_tuple(latitudes, longitudes);
        }
        if (name == "float64")
        {
          auto result = geodesic_points<double>(static_cast<std::size_t>(m));
          NpArr latitudes = new_array<double>({result.first.size()});
          NpArr longitudes = new_array<double>({result.second.size()});
          auto lat = DUCC0_NAMESPACE::detail_pybind::to_vfmav<double>(latitudes);
          auto lon = DUCC0_NAMESPACE::detail_pybind::to_vfmav<double>(longitudes);
          std::copy(result.first.begin(), result.first.end(), lat.data());
          std::copy(result.second.begin(), result.second.end(), lon.data());
          return py::make_tuple(latitudes, longitudes);
        }
        throw py::type_error("dtype must be numpy.float32 or numpy.float64");
      },
      "m"_a, py::kw_only(), "dtype"_a = py::none(),
      py::sig(
          R"(def getgeodesicpts(m: int, *, dtype: object | None = None) -> tuple[numpy.typing.NDArray[numpy.float32] | numpy.typing.NDArray[numpy.float64], numpy.typing.NDArray[numpy.float32] | numpy.typing.NDArray[numpy.float64]])"),
      R"doc(Return the vertices of a subdivided icosahedral geodesic grid.

Parameters
----------
m : int
    Number of points per icosahedral edge, including endpoints;
    must be at least 1.
dtype : numpy dtype, optional
    float64 (default) or float32 output coordinates.

Returns
-------
(latitudes, longitudes) : tuple of ndarray
    Coordinate arrays in degrees with 10 * (m - 1)**2 + 2 points.
    Coordinates are computed in float64 before optional conversion.)doc");

  m.def(
      "specintrp",
      [](double longitude, const CNpArr& spec, const CNpArr& legfuncs)
      {
        if (!std::isfinite(longitude))
          throw py::value_error("lon must be finite");
        if (spec.ndim() != 1 || legfuncs.ndim() != 1)
          throw py::value_error("dataspec and legfuncs must be rank-one arrays");
        if (spec.shape(0) == 0 || legfuncs.shape(0) == 0)
          throw py::value_error("dataspec and legfuncs must not be empty");
        if (DUCC0_NAMESPACE::detail_pybind::isPyarr<std::complex<float>>(spec))
        {
          validate_array_dtype<float>(legfuncs, "numpy.float32");
          return py::float_(spectral_interpolate<float>(longitude, spec, legfuncs));
        }
        if (DUCC0_NAMESPACE::detail_pybind::isPyarr<std::complex<double>>(spec))
        {
          validate_array_dtype<double>(legfuncs, "numpy.float64");
          return py::float_(spectral_interpolate<double>(longitude, spec, legfuncs));
        }
        throw py::type_error(
            "dataspec must have dtype numpy.complex64 or numpy.complex128");
      },
      "lon"_a, "dataspec"_a, "legfuncs"_a,
      py::sig(
          R"(def specintrp(lon: float, dataspec: numpy.typing.NDArray[numpy.complex64] | numpy.typing.NDArray[numpy.complex128], legfuncs: numpy.typing.NDArray[numpy.float32] | numpy.typing.NDArray[numpy.float64]) -> float)"),
      R"doc(Evaluate a scalar spherical harmonic expansion at one point.

Parameters
----------
lon : float
    Longitude in degrees.
dataspec : ndarray
    One-dimensional complex coefficients in triangular packing.
legfuncs : ndarray
    One-dimensional associated Legendre values at the point's
    latitude, as returned by legendre. Precision must match
    dataspec (float32/complex64 or float64/complex128).

Returns
-------
float
    Real interpolated value at the specified longitude and
    latitude represented by legfuncs.)doc");

  m.def(
      "specintrp_many",
      [](const CNpArr& longitude_array, const CNpArr& spec, const CNpArr& legfuncs)
      {
        auto longitudes = real_vector(longitude_array, "longitudes");
        if (spec.ndim() != 1 || legfuncs.ndim() != 2)
          throw py::value_error("dataspec must be rank one and legfuncs rank two");
        if (spec.shape(0) == 0)
          throw py::value_error("dataspec must not be empty");
        truncation_from_count(spec.shape(0));
        if (legfuncs.shape(0) != longitudes.size() ||
            legfuncs.shape(1) != spec.shape(0))
          throw py::value_error("legfuncs must have shape (npoints, ncoefficients)");
        for (const auto longitude : longitudes)
          if (!std::isfinite(longitude))
            throw py::value_error("longitudes must be finite");
        if (DUCC0_NAMESPACE::detail_pybind::isPyarr<std::complex<float>>(spec))
        {
          validate_array_dtype<float>(legfuncs, "numpy.float32");
          return spectral_interpolate_many<float>(longitudes, spec, legfuncs);
        }
        if (DUCC0_NAMESPACE::detail_pybind::isPyarr<std::complex<double>>(spec))
        {
          validate_array_dtype<double>(legfuncs, "numpy.float64");
          return spectral_interpolate_many<double>(longitudes, spec, legfuncs);
        }
        throw py::type_error(
            "dataspec must have dtype numpy.complex64 or numpy.complex128");
      },
      "longitudes"_a, "dataspec"_a, "legfuncs"_a,
      py::sig(
          R"(def specintrp_many(longitudes: numpy.typing.NDArray[numpy.float32] | numpy.typing.NDArray[numpy.float64], dataspec: numpy.typing.NDArray[numpy.complex64] | numpy.typing.NDArray[numpy.complex128], legfuncs: numpy.typing.NDArray[numpy.float32] | numpy.typing.NDArray[numpy.float64]) -> numpy.typing.NDArray[numpy.float32] | numpy.typing.NDArray[numpy.float64])"),
      R"doc(Evaluate one scalar spherical harmonic expansion at several points.

Parameters
----------
longitudes : ndarray
    One-dimensional longitude array in degrees, length npoints.
dataspec : ndarray
    One-dimensional complex coefficients in triangular packing.
legfuncs : ndarray
    Associated Legendre values of shape (npoints, ncoeff), e.g.,
    from legendre_many. Precision must match dataspec.

Returns
-------
ndarray
    Interpolated values of shape (npoints,), with float32 or float64
    precision matching dataspec and legfuncs.)doc");
}

}  // namespace SPHARM_CPP_NAMESPACE
