# SPDX-License-Identifier: GPL-2.0-or-later
cmake_minimum_required(VERSION 3.30.3)
project(pyspharm_cpp LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)
set(CMAKE_CXX_VISIBILITY_PRESET hidden)
set(CMAKE_INTERPROCEDURAL_OPTIMIZATION OFF)

if(NOT DEFINED SPHARM_CPP_BUILD_PROFILE)
  set(SPHARM_CPP_BUILD_PROFILE "auto" CACHE STRING "Build profile: auto, portable, native, or multiarch")
endif()
set_property(CACHE SPHARM_CPP_BUILD_PROFILE PROPERTY STRINGS auto portable native multiarch)

option(SPHARM_CPP_ENABLE_SANITIZERS
  "Enable AddressSanitizer and UndefinedBehaviorSanitizer (GCC/Clang only)"
  OFF)
option(SPHARM_CPP_ENABLE_COVERAGE
  "Instrument first-party C++ sources for gcov coverage (GCC only)"
  OFF)
option(SPHARM_CPP_STABLE_ABI
  "Build standard CPython extensions with the linked abi3 ABI (CPython 3.12+)"
  OFF)
if(SPHARM_CPP_ENABLE_SANITIZERS AND NOT CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
  message(FATAL_ERROR "SPHARM_CPP_ENABLE_SANITIZERS requires GCC or Clang")
endif()
if(SPHARM_CPP_ENABLE_COVERAGE AND NOT CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
  message(FATAL_ERROR "SPHARM_CPP_ENABLE_COVERAGE requires GCC/gcov")
endif()
if(SPHARM_CPP_ENABLE_COVERAGE AND SPHARM_CPP_ENABLE_SANITIZERS)
  message(FATAL_ERROR "Coverage and sanitizer instrumentation cannot be enabled together")
endif()

if(SPHARM_CPP_BUILD_PROFILE STREQUAL "auto")
  if(CMAKE_SYSTEM_NAME STREQUAL "Linux" AND CMAKE_SYSTEM_PROCESSOR MATCHES "^(x86_64|amd64|AMD64)$")
    set(_profile "multiarch")
  else()
    set(_profile "portable")
  endif()
else()
  set(_profile "${SPHARM_CPP_BUILD_PROFILE}")
endif()

if(NOT _profile MATCHES "^(portable|native|multiarch)$")
  message(FATAL_ERROR "Unknown SPHARM_CPP_BUILD_PROFILE='${SPHARM_CPP_BUILD_PROFILE}'")
endif()
if(SPHARM_CPP_ENABLE_COVERAGE AND NOT _profile STREQUAL "portable")
  message(FATAL_ERROR "SPHARM_CPP_ENABLE_COVERAGE requires the portable build profile")
endif()
if(_profile STREQUAL "native" AND NOT CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
  message(FATAL_ERROR "The native profile requires GCC or Clang support for -march=native; use portable with this compiler")
endif()
if(_profile STREQUAL "multiarch" AND NOT (CMAKE_SYSTEM_NAME STREQUAL "Linux" AND CMAKE_SYSTEM_PROCESSOR MATCHES "^(x86_64|amd64|AMD64)$"))
  message(FATAL_ERROR "The multiarch profile currently requires Linux x86-64; use portable or native on this platform")
endif()

message(STATUS "pyspharm-cpp build profile: ${_profile}")
message(STATUS "Pinned DUCC source: ${CMAKE_CURRENT_SOURCE_DIR}/third_party/ducc")

find_package(Python 3.11 REQUIRED COMPONENTS Interpreter Development.Module
  OPTIONAL_COMPONENTS Development.SABIModule)
find_package(nanobind CONFIG REQUIRED)
find_package(Threads REQUIRED)
execute_process(
  COMMAND "${Python_EXECUTABLE}" -c
    "import sysconfig; print(int(bool(sysconfig.get_config_var('Py_GIL_DISABLED'))))"
  OUTPUT_STRIP_TRAILING_WHITESPACE
  OUTPUT_VARIABLE SPHARM_CPP_PYTHON_FREE_THREADED
  COMMAND_ERROR_IS_FATAL ANY)
if(SPHARM_CPP_PYTHON_FREE_THREADED STREQUAL "1")
  message(STATUS "Building for free-threaded CPython")
  set(SPHARM_CPP_IS_FREE_THREADED ON)
else()
  message(STATUS "Building for standard CPython")
  set(SPHARM_CPP_IS_FREE_THREADED OFF)
endif()
if(SPHARM_CPP_STABLE_ABI)
  if(Python_VERSION VERSION_LESS 3.12)
    message(FATAL_ERROR "SPHARM_CPP_STABLE_ABI requires CPython 3.12 or newer")
  endif()
  if(SPHARM_CPP_IS_FREE_THREADED)
    message(FATAL_ERROR "Linked abi3 is unavailable for free-threaded Python")
  endif()
  if(NOT Python_INTERPRETER_ID STREQUAL "Python" OR
      NOT TARGET Python::SABIModule)
    message(FATAL_ERROR
      "SPHARM_CPP_STABLE_ABI requires CPython and Python::SABIModule")
  endif()
  if(NOT DEFINED SKBUILD_SABI_VERSION OR
      SKBUILD_SABI_VERSION VERSION_LESS 3.12)
    message(FATAL_ERROR
      "SPHARM_CPP_STABLE_ABI requires scikit-build-core wheel.py-api=cp312 or newer")
  endif()
  message(STATUS "Building for the linked CPython stable ABI (abi3)")
elseif(DEFINED SKBUILD_SABI_VERSION AND
    NOT SKBUILD_SABI_VERSION STREQUAL "")
  message(FATAL_ERROR
    "wheel.py-api requests an abi3 wheel, but SPHARM_CPP_STABLE_ABI is OFF")
endif()
if(WIN32 AND SPHARM_CPP_IS_FREE_THREADED
    AND NOT "Py_GIL_DISABLED=1" IN_LIST Python_DEFINITIONS)
  message(FATAL_ERROR
    "CMake did not provide Py_GIL_DISABLED=1 for the Windows free-threaded build")
endif()
execute_process(
  COMMAND "${Python_EXECUTABLE}" -m nanobind --include_dir
  OUTPUT_STRIP_TRAILING_WHITESPACE
  OUTPUT_VARIABLE NANOBIND_INCLUDE_DIR
  COMMAND_ERROR_IS_FATAL ANY)

set(DUCC_ROOT "${CMAKE_CURRENT_SOURCE_DIR}/third_party/ducc")
set(SPHARM_CPP_DUCC_TAG "ducc0_0_42_0")
set(SPHARM_CPP_DUCC_REVISION "42ce8ab5c5d79e0e2fd9af4a09c2ea4693520e1b")
message(STATUS "Pinned DUCC revision: ${SPHARM_CPP_DUCC_TAG} (${SPHARM_CPP_DUCC_REVISION})")
set(DUCC_SHT_SOURCES
  "${DUCC_ROOT}/src/ducc0/fft/fft_inst1.cc"
  "${DUCC_ROOT}/src/ducc0/fft/fft_inst2.cc"
  "${DUCC_ROOT}/src/ducc0/infra/mav.cc"
  "${DUCC_ROOT}/src/ducc0/infra/string_utils.cc"
  "${DUCC_ROOT}/src/ducc0/infra/threading.cc"
  "${DUCC_ROOT}/src/ducc0/math/geom_utils.cc"
  "${DUCC_ROOT}/src/ducc0/math/gridding_kernel.cc"
  "${DUCC_ROOT}/src/ducc0/math/gl_integrator.cc"
  "${DUCC_ROOT}/src/ducc0/sht/sht.cc")
set(SPHARM_CPP_IMPLEMENTATION_SOURCES
  src/cpp/bindings.cpp
  src/cpp/transforms.cpp
  src/cpp/numerics.cpp)
if(SPHARM_CPP_ENABLE_COVERAGE)
  set_source_files_properties(
    src/cpp/module.cpp
    ${SPHARM_CPP_IMPLEMENTATION_SOURCES}
    PROPERTIES COMPILE_OPTIONS
      "-O0;-g;--coverage;-fprofile-arcs;-ftest-coverage;-fprofile-abs-path;-fprofile-update=atomic")
endif()

# These DUCC translation units contain dynamic initialization that can run
# before CPU dispatch (the thread pool and gridding-kernel database). Keep
# their initialization code at the x86-64 baseline in every namespace.
if(_profile STREQUAL "multiarch" AND CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
  set_source_files_properties(
    "${DUCC_ROOT}/src/ducc0/infra/threading.cc"
    "${DUCC_ROOT}/src/ducc0/math/gridding_kernel.cc"
    PROPERTIES COMPILE_OPTIONS "-march=x86-64")
endif()

function(add_spharm_backend target profile namespace_name)
  add_library(${target} OBJECT
    ${DUCC_SHT_SOURCES}
    ${SPHARM_CPP_IMPLEMENTATION_SOURCES})
  target_include_directories(${target} PRIVATE
    "${DUCC_ROOT}/src"
    "${DUCC_ROOT}/python"
    "${NANOBIND_INCLUDE_DIR}"
    ${Python_INCLUDE_DIRS}
    "${CMAKE_CURRENT_SOURCE_DIR}/src/cpp")
  target_compile_definitions(${target} PRIVATE
    DUCC0_USE_NANOBIND
    DUCC0_NAMESPACE=${namespace_name}
    SPHARM_CPP_NAMESPACE=${namespace_name}
    SPHARM_CPP_PROFILE_NAME=\"${profile}\")
  if(SPHARM_CPP_STABLE_ABI)
    # nanobind applies this to its linked library and the module target. The
    # profile object libraries also compile bindings and must use the same
    # limited Python API to keep the final abi3 module valid.
    target_compile_definitions(${target} PRIVATE Py_LIMITED_API=0x030C0000)
  endif()
  if(SPHARM_CPP_IS_FREE_THREADED)
    target_compile_definitions(${target} PRIVATE
      NB_FREE_THREADED
      ${Python_DEFINITIONS})
  endif()
  nanobind_compile_options(${target})
  if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options(${target} PRIVATE
      -O3 -fno-fast-math -fno-finite-math-only
      -Wall -Wextra -Wpedantic)
    if(profile STREQUAL "portable" AND CMAKE_SYSTEM_PROCESSOR MATCHES "^(x86_64|amd64|AMD64)$")
      target_compile_options(${target} PRIVATE -march=x86-64)
    elseif(profile STREQUAL "native")
      target_compile_options(${target} PRIVATE -march=native)
    elseif(profile STREQUAL "x86-64-v1")
      target_compile_options(${target} PRIVATE -march=x86-64)
    elseif(profile STREQUAL "x86-64-v3")
      target_compile_options(${target} PRIVATE -march=x86-64-v3)
    elseif(profile STREQUAL "x86-64-v4")
      target_compile_options(${target} PRIVATE -march=x86-64-v4)
    endif()
    if(SPHARM_CPP_ENABLE_SANITIZERS)
      target_compile_options(${target} PRIVATE
        -fsanitize=address,undefined -fno-omit-frame-pointer)
    endif()
  endif()
endfunction()

set(SPHARM_CPP_NANOBIND_OPTIONS NB_STATIC FREE_THREADED)
if(SPHARM_CPP_STABLE_ABI)
  list(APPEND SPHARM_CPP_NANOBIND_OPTIONS STABLE_ABI)
endif()
nanobind_add_module(_spharm_cpp ${SPHARM_CPP_NANOBIND_OPTIONS}
  src/cpp/module.cpp
  third_party/ducc/python/multiarch.cc)
target_include_directories(_spharm_cpp PRIVATE
  "${DUCC_ROOT}/src"
  "${DUCC_ROOT}/python"
  "${CMAKE_CURRENT_SOURCE_DIR}/src/cpp")
target_link_libraries(_spharm_cpp PRIVATE Threads::Threads)
target_compile_definitions(_spharm_cpp PRIVATE
  SPHARM_CPP_BUILD_PROFILE=\"${_profile}\"
  SPHARM_CPP_DUCC_TAG=\"${SPHARM_CPP_DUCC_TAG}\"
  SPHARM_CPP_DUCC_REVISION=\"${SPHARM_CPP_DUCC_REVISION}\")
if(_profile STREQUAL "native")
  target_compile_definitions(_spharm_cpp PRIVATE SPHARM_CPP_IS_NATIVE)
endif()
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
  target_compile_options(_spharm_cpp PRIVATE
    -O2 -fno-fast-math -Wall -Wextra -Wpedantic)
  if(CMAKE_SYSTEM_PROCESSOR MATCHES "^(x86_64|amd64|AMD64)$")
    target_compile_options(_spharm_cpp PRIVATE -march=x86-64)
  endif()
  if(SPHARM_CPP_ENABLE_SANITIZERS)
    target_compile_options(_spharm_cpp PRIVATE
      -fsanitize=address,undefined -fno-omit-frame-pointer)
    target_link_options(_spharm_cpp PRIVATE -fsanitize=address,undefined)
  endif()
endif()
if(SPHARM_CPP_ENABLE_COVERAGE)
  target_link_options(_spharm_cpp PRIVATE --coverage)
endif()

if(_profile STREQUAL "multiarch")
  set(CMAKE_INTERPROCEDURAL_OPTIMIZATION OFF)
  target_compile_definitions(_spharm_cpp PRIVATE DUCC0_MULTIARCH)
  add_spharm_backend(spharm_backend_v1 x86-64-v1 spharm_cpp_v1)
  add_spharm_backend(spharm_backend_v3 x86-64-v3 spharm_cpp_v3)
  add_spharm_backend(spharm_backend_v4 x86-64-v4 spharm_cpp_v4)
  target_compile_definitions(spharm_backend_v1 PRIVATE DUCC0_MULTIARCH)
  target_compile_definitions(spharm_backend_v3 PRIVATE DUCC0_MULTIARCH)
  target_compile_definitions(spharm_backend_v4 PRIVATE DUCC0_MULTIARCH)
  target_sources(_spharm_cpp PRIVATE
    $<TARGET_OBJECTS:spharm_backend_v1>
    $<TARGET_OBJECTS:spharm_backend_v3>
    $<TARGET_OBJECTS:spharm_backend_v4>)
else()
  add_spharm_backend(spharm_backend ${_profile} spharm_cpp_single)
  target_compile_definitions(spharm_backend PRIVATE SPHARM_CPP_SINGLE_BACKEND)
  target_sources(_spharm_cpp PRIVATE $<TARGET_OBJECTS:spharm_backend>)
endif()

install(TARGETS _spharm_cpp LIBRARY DESTINATION spharm_cpp RUNTIME DESTINATION spharm_cpp)
