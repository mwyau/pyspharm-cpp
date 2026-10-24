<!-- SPDX-License-Identifier: GPL-2.0-or-later -->

# Repository constraints

These apply to all implementation and review work in this repository.

## Numerical boundary

All analysis must respect the longitude sampling limit `lmax <= floor(nlon/2)`.
Regular-grid analysis must additionally use DUCC's native transform and its
supported `lmax <= nlat-2` limit. Gaussian analysis may use `lmax=nlat-1`
when the longitude limit permits it. Explicit analysis requests above either
limit must be rejected, never silently clamped. Never implement or emulate
SPHEREPACK `zbar` projection, maximum-degree correction, pseudoinverse,
least-squares fallback, or another transform to reproduce regular-grid
`lmax=nlat-1` analysis. This exclusion is permanent unless the owner
explicitly changes it. Supported regular-grid synthesis retains its separate
valid limits.

## DUCC and native integration

- Use only the unmodified DUCC source in the `third_party/ducc` Git submodule;
  preserve its commit pin and upstream provenance.
- Do not vendor/copy its source into another tracked tree, add a development
  vendoring script, retrieve it during wheel builds, or depend on its Python
  package at runtime.
- Keep our C++ driver, PySpharm bindings and compatibility code outside the
  submodule. Reuse DUCC's SHT/FFT, SIMD, threading, CPU detection and multiarch
  dispatch.
- Do not add an independent FFT, transform backend, or CPU-dispatch framework.
- Multiarch release defaults to Linux x86-64 v1/v3/v4 with baseline driver and
  isolated namespaces; use portable builds by default elsewhere. Native builds
  are explicitly host-specific. Any initialization that can run before CPU
  dispatch, including DUCC's thread-pool initialization, must use baseline
  instructions in every backend namespace.

## Documentation and verification

- Follow [write-like-a-scientist](https://github.com/mwyau/write-like-a-scientist/blob/main/skills/write-like-a-scientist/SKILL.md) when editing scientific text. Preserve equations and claims; distinguish tested behavior from expected behavior.
- Keep README and CONTRIBUTING focused on users and contributors. Put API explanations in `docs/` and implementation constraints here. Do not add work-log or agent-handoff files.
- Preserve the parity manifest: every public contract and approved exception must have a behavior test. Keep validation manual-only and GitHub Actions pinned to exact patch releases.

## Python build and release

- Use uv and a committed `uv.lock` for development environments, dependency
  groups, tests, editable installs, and ordinary project wheel/sdist builds.
- Keep scikit-build-core as the only PEP 517 backend, CMake for native builds,
  nanobind for bindings, and NumPy as the only required Python runtime
  dependency. Do not add Hatchling, setuptools, uv_build or the `ducc0` Python
  package without a demonstrated requirement.
- Keep ordinary `pip install` and isolated PEP 517 builds independent of uv.
- Include project and bundled-source license/attribution notices in both
  source distributions and wheels.
- Release preparation may create CI and artifacts, but do not publish to PyPI,
  create release tags, make the repository public, or force-push.
