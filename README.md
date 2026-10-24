<!-- SPDX-License-Identifier: GPL-2.0-or-later -->

# pyspharm-cpp

`pyspharm-cpp` implements PySpharm's NumPy interface for scalar and vector
spherical harmonic transforms (SHTs), regridding, Legendre functions, and
spectral interpolation. The extension uses C++17, nanobind, and the
[DUCC](https://gitlab.mpcdf.mpg.de/mtr/ducc) transform implementation.
The [compatibility comparison](docs/pyspharm-compatibility.md) documents
differences in precision, spectral truncation, and Legendre evaluation.

For DUCC's native C++ SHT API or its `ducc0` Python API, use
[DUCC](https://gitlab.mpcdf.mpg.de/mtr/ducc) directly. For atmospheric
filtering, regridding, and wind diagnostics with labeled Xarray fields, use
[spharmgrid](https://github.com/mwyau/spharmgrid). Its Xarray operations use
DUCC on the CPU. Optional JAX/S2FFT and PyTorch/torch-harmonics interfaces
can run on GPUs with an appropriate accelerator installation.

`pyspharm-cpp` is a separate project from PySpharm and DUCC. Its runtime
Python dependency is NumPy; SPHEREPACK and the `ducc0` Python package are
not required.

## Install

Install a wheel when one is available for your Python version and platform:

```sh
python -m pip install pyspharm-cpp
```

The only required Python runtime dependency is NumPy. Python 3.11 or newer and
a C++17 compiler are required when pip builds from source. See the
[installation and build notes](docs/quickstart.md) for build profiles and
ordinary PEP 517 installs.

The wheel configuration uses `cp311` for CPython 3.11 and
`cp312-abi3` for CPython 3.12 and newer. Free-threaded CPython 3.14 and
3.15 use separate `cp314t` and `cp315t` wheels. See
[the build and validation notes](CONTRIBUTING.md) for supported platforms and
test results.

## Quick example

```python
import numpy as np
import spharm_cpp as spharm

nlon, nlat = 32, 16
transform = spharm.Spharmt(nlon, nlat, gridtype="gaussian")
latitudes, _ = spharm.gaussian_lats_wts(nlat)
longitudes = np.arange(nlon) * (360.0 / nlon)

# A smooth, degree-one scalar field on a Gaussian latitude-longitude grid.
field = np.cos(np.deg2rad(latitudes))[:, None] * np.cos(np.deg2rad(longitudes))[None, :]
spectrum = transform.grdtospec(field, ntrunc=1)
reconstructed = transform.spectogrd(spectrum)
np.testing.assert_allclose(reconstructed, field, rtol=2e-12, atol=2e-12)
```

See the [quickstart](docs/quickstart.md), [API reference](docs/api.md), and
[numerical conventions](docs/conventions.md) for grid, coefficient, precision,
and transform details. The [PySpharm comparison](docs/pyspharm-compatibility.md)
records tested operations, observed differences, and comparison tolerances.

## Numerical behavior

The default map and coefficient types are `numpy.float64` and
`numpy.complex128`. Pass `dtype=numpy.float32` to `Spharmt` for
`float32`/`complex64` transforms. Grid arrays have shape `(nlat, nlon)`;
independent fields can use `(nlat, nlon, nfields)`.

Analysis respects both the latitude grid and longitude sampling limits:

- Regular Clenshaw–Curtis analysis supports `ntrunc <= nlat - 2`.
- Gaussian analysis supports `ntrunc <= nlat - 1`.
- Both grids require `ntrunc <= floor(nlon / 2)`.

The default analysis truncation is the minimum of these limits. An explicit
request above a limit raises `ValueError`; it is not silently clamped.
Synthesis accepts triangular spectra through degree `nlat - 1`, including
some modes that analysis cannot uniquely recover. `legfunc="computed"` is the
default; `legfunc="stored"` raises `ValueError`.

## Build profiles and threading

`SPHARM_CPP_BUILD_PROFILE` selects `portable`, `native`, or `multiarch`. The
default `auto` selects `multiarch` on Linux x86-64 and `portable` elsewhere.
Linux x86-64 multiarch builds include DUCC x86-64-v1, v3 and v4 variants and
dispatch using DUCC CPU detection. Set `DUCC0_MAX_PSABI_LEVEL` before import to
limit dispatch during fallback checks. `cpu_info()` matches DUCC's CPU and
dispatch metadata schema.

Transform calls use one DUCC worker by default. Set `DUCC0_NUM_THREADS` (or
`OMP_NUM_THREADS`) before importing to choose a worker count; `0` delegates to
DUCC's host default. Transform calls release the Python GIL while native work
runs. Inputs must not be changed until a call returns.

## Development

The repository uses `uv` and a committed `uv.lock` for development and CI.
Install the tools and editable package, then run the Python tests:

```sh
uv python install 3.11 3.12 3.13 3.14 3.15
uv sync --locked --python 3.14 --all-groups
uv run --locked python -m pytest -q
```

Free-threaded CPython 3.14 and 3.15 builds are also included in concurrency
and wheel validation. Their platform coverage is tracked separately from the
regular CPython matrix.

Build a portable wheel and source distribution with the locked build tools:

```sh
uv run --locked --no-sync python -m build --wheel --sdist --no-isolation \
  --config-setting cmake.define.SPHARM_CPP_BUILD_PROFILE=portable
```

The source distribution includes DUCC headers and compile sources used by the
extension, plus only the required binding and NUFFT headers, its multiarch
detector, and license. It excludes unrelated DUCC source files and components
such as Healpix and wgridder. An ordinary PEP 517 installation remains
independent of `uv`:

```sh
python -m pip install .
```

The `native` profile targets the build host and is not portable. `multiarch` is
currently supported only on Linux x86-64. Portable builds do not request a
host-specific instruction set. See [contributing](CONTRIBUTING.md) for
formatting, lint, typing, documentation, coverage and native-code checks.

## Source and license

DUCC is included as an unmodified Git submodule at release tag
`ducc0_0_42_0` (commit `42ce8ab5c5d79e0e2fd9af4a09c2ea4693520e1b`). See
[`THIRD_PARTY_NOTICES.md`](THIRD_PARTY_NOTICES.md) for included components and
license details. Redistributions must retain the applicable license notices.

This project is licensed under GPL-2.0-or-later. The public API and
conventions are based on PySpharm documentation; this package is not an
official PySpharm or DUCC project.
