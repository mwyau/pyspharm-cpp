<!-- SPDX-License-Identifier: GPL-2.0-or-later -->

# Quickstart

## Install a wheel

After a package release is available, install it with pip when a wheel is
available for your platform and Python version:

```console
python -m pip install pyspharm-cpp
```

The runtime dependency is NumPy. A source build requires Python 3.11 or newer,
a C++17 compiler, CMake and the build tools declared by the project.

## Analyze and synthesize a field

The example analyzes a degree-one spherical harmonic sampled on a Gaussian
grid, then reconstructs the field from its coefficients:

```python
import numpy as np
import spharm_cpp as spharm

nlon, nlat = 32, 16
transform = spharm.Spharmt(nlon, nlat, gridtype="gaussian")
latitudes, _ = spharm.gaussian_lats_wts(nlat)
longitudes = np.arange(nlon) * (360.0 / nlon)

field = np.cos(np.deg2rad(latitudes))[:, None] * np.cos(np.deg2rad(longitudes))[None, :]
spectrum = transform.grdtospec(field, ntrunc=1)
reconstructed = transform.spectogrd(spectrum)
np.testing.assert_allclose(reconstructed, field, rtol=2e-12, atol=2e-12)
```

`Spharmt` accepts either `gridtype="gaussian"` or `gridtype="regular"`.
Regular grids use Clenshaw–Curtis latitude rings. Gaussian latitude and
quadrature arrays are returned in north-to-south order by
`gaussian_lats_wts`.

## Vector operations

The vector methods return eastward (zonal) and northward (meridional)
components, in that order. This example differentiates the degree-one,
order-one scalar harmonic and analyzes the resulting wind components:

```python
order, degree = spharm.getspecindx(1)
potential = np.zeros(degree.size, dtype=np.complex128)
potential[(order == 1) & (degree == 1)] = 1.0
eastward, northward = transform.getgrad(potential)

vorticity, divergence = transform.getvrtdivspec(eastward, northward, ntrunc=1)
eastward_again, northward_again = transform.getuv(vorticity, divergence)
```

`getpsichi(ugrid, vgrid)` returns streamfunction and velocity potential,
in that order. Spherical derivatives use the radius `rsphere`; use consistent
units for the input field, wind components, and radius.

## Precision and batched fields

The default `dtype` is `numpy.float64`, paired with `numpy.complex128` spectral
coefficients. Select single precision with `dtype=numpy.float32`; real and
complex arrays must use the matching precision. Batched map inputs have shape
`(nlat, nlon, nfields)` and batched spectra have shape `(ncoeff, nfields)`.

## Install from a checkout

For a local editable installation, initialize the pinned DUCC submodule and
sync the locked environment:

```sh
git submodule update --init --recursive
uv sync --locked --python 3.14 --all-groups
uv run --locked python -m pytest -q
```

The default build profile is selected by platform. To build a portable wheel,
use:

```sh
uv run --locked --no-sync python -m build --wheel --no-isolation \
  --config-setting cmake.define.SPHARM_CPP_BUILD_PROFILE=portable
```

For a CPython 3.12+ standard build, the linked nanobind stable ABI can be
requested explicitly:

```sh
uv python install 3.12
uv sync --locked --python 3.12 --group build
uv run --locked --no-sync --python 3.12 python -m build \
  --wheel --no-isolation \
  --config-setting cmake.define.SPHARM_CPP_BUILD_PROFILE=portable \
  --config-setting cmake.define.SPHARM_CPP_STABLE_ABI=ON \
  --config-setting wheel.py-api=cp312
```

This produces a `cp312-abi3` wheel. Keep a version-specific `cp311` wheel for
Python 3.11. Free-threaded builds use the separate `cp314t` and `cp315t` ABI
tags; the linked `abi3` build mode does not apply to them.

See [numerical conventions](conventions.md) for the transform limits and
coefficient layout.
