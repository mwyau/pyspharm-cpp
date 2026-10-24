// SPDX-License-Identifier: GPL-2.0-or-later
#include <nanobind/nanobind.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/vector.h>

#include <stdexcept>
#include <string>
#include <vector>

#include "multiarch.h"

namespace py = nanobind;

#ifdef DUCC0_MULTIARCH
namespace spharm_cpp_v1
{
void add_api(py::module_& m);
}
namespace spharm_cpp_v3
{
void add_api(py::module_& m);
}
namespace spharm_cpp_v4
{
void add_api(py::module_& m);
}
#else
namespace spharm_cpp_single
{
void add_api(py::module_& m);
}
#endif

namespace
{

const char* single_profile_name()
{
#ifdef SPHARM_CPP_IS_NATIVE
  return "native";
#else
  return "portable";
#endif
}

void add_cpu_diagnostics(py::module_& m, const ducc0_multiarch::cpu_capabilities& cpu,
                         bool multiarch, int active_profile, int configured_limit)
{
  const std::string active_name =
      multiarch ? ducc0_multiarch::profile_name(active_profile) : single_profile_name();
  std::vector<std::string> compiled;
  std::vector<std::string> available;
  if (multiarch)
  {
    for (int p : ducc0_multiarch::compiled_profiles())
      compiled.emplace_back(ducc0_multiarch::profile_name(p));
    for (int p : ducc0_multiarch::available_profiles(
             ducc0_multiarch::ducc_compiled_profiles_mask, cpu.x86_psabi_level))
      available.emplace_back(ducc0_multiarch::profile_name(p));
  }
  else
  {
    compiled.emplace_back(active_name);
    available.emplace_back(active_name);
  }
  m.def(
      "cpu_info",
      [=]()
      {
        py::dict result;
        result["architecture"] = ducc0_multiarch::architecture_name();
        result["features"] = cpu.features;
        result["multiarch"] = multiarch;
        if (multiarch)
        {
          result["compiled_profiles"] = compiled;
          result["available_profiles"] = available;
          result["configured_limit"] = ducc0_multiarch::profile_name(configured_limit);
          result["active_profile"] = active_name;
        }
        return result;
      },
      py::sig("def cpu_info() -> dict[str, object]"),
      "Return information about the current CPU and DUCC dispatch configuration.\n\n"
      "The ``features`` entry lists SIMD/vector features recognized by DUCC that are "
      "available to the current process.\n\n"
      "For multiarch builds, the result also reports the compiled and available "
      "x86-64 profiles, the configured profile limit, and the profile selected for "
      "the current process.");
}

}  // namespace

// The binding macro uses signed bit flags internally; the application code does not.
NB_MODULE(_spharm_cpp, m)  // NOLINT(bugprone-signed-bitwise)
{
  const auto cpu = ducc0_multiarch::detect_cpu_capabilities();
#ifdef DUCC0_MULTIARCH
  const auto state = ducc0_multiarch::current_profile_state(cpu.x86_psabi_level);
  switch (state.active_profile)
  {
    case 4:
      spharm_cpp_v4::add_api(m);
      break;
    case 3:
      spharm_cpp_v3::add_api(m);
      break;
    case 1:
      spharm_cpp_v1::add_api(m);
      break;
    default:
      throw std::runtime_error("no compatible DUCC multiarch profile");
  }
  add_cpu_diagnostics(m, cpu, true, state.active_profile, state.configured_limit);
#else
  spharm_cpp_single::add_api(m);
  add_cpu_diagnostics(m, cpu, false, 0, 0);
#endif
}
