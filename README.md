<!-- SPDX-License-Identifier: GPL-2.0-or-later -->

# pyspharm-cpp

[![CI](https://github.com/mwyau/pyspharm-cpp/actions/workflows/ci.yml/badge.svg)](https://github.com/mwyau/pyspharm-cpp/actions/workflows/ci.yml)
[![Docs](https://readthedocs.org/projects/pyspharm-cpp/badge/?version=latest)](https://pyspharm-cpp.readthedocs.io/en/latest/)
[![Codecov](https://codecov.io/github/mwyau/pyspharm-cpp/graph/badge.svg)](https://codecov.io/github/mwyau/pyspharm-cpp)

`pyspharm-cpp` implements the PySpharm Python API for scalar and vector spherical harmonic transforms (SHTs), spectral regridding, Legendre functions, and interpolation. It uses DUCC's C++ transforms with nanobind and requires only NumPy at runtime.

For DUCC's native C++ or Python transform API, use [DUCC](https://gitlab.mpcdf.mpg.de/mtr/ducc) directly. For higher-level atmospheric-science operations, including vorticity, divergence, streamfunction, velocity potential, filtering, and regridding, use [spharmgrid](https://github.com/mwyau/spharmgrid). Spharmgrid supports Xarray fields with a DUCC CPU backend and offers separate JAX/S2FFT and PyTorch/torch-harmonics interfaces for GPU calculations.

## Install

Install from PyPI when a distribution is available for your platform:

```sh
python -m pip install pyspharm-cpp
```

Python 3.11 or newer is required. Source builds need a C++17 compiler. The Python import is `spharm_cpp`, not `spharm`.

## Example

```python
import numpy as np
import spharm_cpp as spharm

nlon, nlat = 32, 16
sht = spharm.Spharmt(nlon, nlat, gridtype="gaussian")
lat, _ = spharm.gaussian_lats_wts(nlat)
lon = np.arange(nlon) * (360.0 / nlon)

field = np.cos(np.deg2rad(lat))[:, None] * np.cos(np.deg2rad(lon))[None, :]
coefficients = sht.grdtospec(field, ntrunc=1)
reconstructed = sht.spectogrd(coefficients)
np.testing.assert_allclose(reconstructed, field, atol=2e-12, rtol=2e-12)
```

## Differences from PySpharm

- **Backend:** the `spharm_cpp` module uses DUCC for transforms.
- **Defaults:** float64/complex128, with optional float32/complex64; computed Legendre functions only (`legfunc="stored"` is unsupported).
- **Truncation:** regular analysis stops at `min(nlat-2, nlon//2)`, Gaussian analysis at `min(nlat-1, nlon//2)`. Vector analysis also requires `ntrunc <= (nlon-1)//2`.
- **Highest regular degree:** synthesis evaluates odd-order harmonics at `l=nlat-1` that the tested SPHEREPACK synthesis omits. Regular-grid analysis cannot uniquely recover all coefficients at that degree.
- **Zero truncation:** coefficient arrays retain one dimension; `getpsichi(..., ntrunc=0)` returns zero maps instead of raising an exception.
- **Other changes:** `getgeodesicpts` calculates in float64 (with optional float32 output), `getspecindx` returns NumPy `intp`, and `regrid` supports precision conversion. Batched Legendre/interpolation and stricter input checks are additional features.

For even longitude counts, the imaginary scalar Nyquist coefficient is unobservable. See [numerical conventions](docs/conventions.md) and the [PySpharm comparison](docs/pyspharm-compatibility.md) for details, or the [quickstart](docs/quickstart.md) and [API reference](docs/api.md) for usage.

## License

GPL-2.0-or-later. DUCC is included as an unmodified submodule; see [third-party notices](THIRD_PARTY_NOTICES.md). DUCC is used under its respective license.
