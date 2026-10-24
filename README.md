<!-- SPDX-License-Identifier: GPL-2.0-or-later -->

# pyspharm-cpp

[![CI](https://github.com/mwyau/pyspharm-cpp/actions/workflows/ci.yml/badge.svg)](https://github.com/mwyau/pyspharm-cpp/actions/workflows/ci.yml)
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

## Numerical conventions

The default precision is float64 for grids and complex128 for coefficients; pass `dtype=np.float32` to `Spharmt` for float32/complex64. Arrays have shape `(nlat, nlon)` or `(nlat, nlon, nfields)` for multiple fields.

The maximum analysis truncation is `min(nlat-2, nlon//2)` on regular Clenshaw–Curtis grids and `min(nlat-1, nlon//2)` on Gaussian grids. Synthesis accepts degree `nlat-1` where the grid supports it, but not every coefficient is recoverable. Legendre functions are computed rather than stored.

See the [quickstart](docs/quickstart.md), [API reference](docs/api.md), [numerical conventions](docs/conventions.md), and [PySpharm comparison](docs/pyspharm-compatibility.md) for details.

## License

GPL-2.0-or-later. DUCC is included as an unmodified submodule; see [third-party notices](THIRD_PARTY_NOTICES.md). This project is independent of PySpharm and DUCC.
