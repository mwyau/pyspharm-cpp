// SPDX-FileCopyrightText: 2026 Albert M. W. Yau
// SPDX-License-Identifier: GPL-2.0-or-later
#pragma once

#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>

#include <charconv>
#include <complex>
#include <cstddef>
#include <cstdlib>
#include <cstring>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "ducc0/bindings/pybind_utils.h"
#include "ducc0/infra/threading.h"

namespace SPHARM_CPP_NAMESPACE
{

namespace py = nanobind;
using NpArr = DUCC0_NAMESPACE::detail_pybind::NpArr;
using CNpArr = DUCC0_NAMESPACE::detail_pybind::CNpArr;

inline constexpr double sqrt_2pi = 2.506628274631000502415765284811045253;

std::size_t checked_truncation(const py::object&, std::size_t, std::size_t,
                               const char*);
std::size_t truncation_from_count(std::size_t);
std::size_t coefficient_count(std::size_t);
std::vector<double> real_vector(const CNpArr&, const char*);

inline std::size_t packed_coefficient_index(std::size_t lmax, std::size_t m,
                                            std::size_t l)
{
  return m * (2 * lmax - m + 3) / 2 + (l - m);
}

inline std::size_t configured_nthreads()
{
  static const std::size_t requested = []
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
  return DUCC0_NAMESPACE::adjust_nthreads(requested);
}

template <typename T>
inline void validate_array_dtype(const CNpArr& arr, const char* expected)
{
  if (!DUCC0_NAMESPACE::detail_pybind::isPyarr<T>(arr))
    throw py::type_error((std::string("array dtype must be ") + expected).c_str());
}

template <typename T>
inline NpArr new_array(const std::vector<std::size_t>& shape)
{
  return DUCC0_NAMESPACE::detail_pybind::make_Pyarr<T>(shape);
}

void bind_numerics(py::module_& m);

NpArr grid_to_spectral_f32(const CNpArr&, std::size_t, std::size_t, const std::string&,
                           std::size_t);
NpArr grid_to_spectral_f64(const CNpArr&, std::size_t, std::size_t, const std::string&,
                           std::size_t);
NpArr spectral_to_grid_f32(const CNpArr&, std::size_t, std::size_t, const std::string&,
                           std::size_t);
NpArr spectral_to_grid_f64(const CNpArr&, std::size_t, std::size_t, const std::string&,
                           std::size_t);
std::pair<NpArr, NpArr> spectral_gradient_f32(const CNpArr&, std::size_t, std::size_t,
                                              const std::string&, double);
std::pair<NpArr, NpArr> spectral_gradient_f64(const CNpArr&, std::size_t, std::size_t,
                                              const std::string&, double);
std::pair<NpArr, NpArr> analyze_winds_f32(const CNpArr&, const CNpArr&, std::size_t,
                                          std::size_t, const std::string&, double,
                                          std::size_t);
std::pair<NpArr, NpArr> analyze_winds_f64(const CNpArr&, const CNpArr&, std::size_t,
                                          std::size_t, const std::string&, double,
                                          std::size_t);
std::pair<NpArr, NpArr> synthesize_winds_f32(const CNpArr&, const CNpArr&, std::size_t,
                                             std::size_t, const std::string&, double);
std::pair<NpArr, NpArr> synthesize_winds_f64(const CNpArr&, const CNpArr&, std::size_t,
                                             std::size_t, const std::string&, double);
NpArr apply_inverse_laplacian_f32(const CNpArr&, double);
NpArr apply_inverse_laplacian_f64(const CNpArr&, double);
NpArr apply_smoothing_f32(const CNpArr&, const CNpArr&);
NpArr apply_smoothing_f64(const CNpArr&, const CNpArr&);
NpArr convert_spectral_precision_f32_from_f64(const CNpArr&);
NpArr convert_spectral_precision_f64_from_f32(const CNpArr&);

}  // namespace SPHARM_CPP_NAMESPACE
