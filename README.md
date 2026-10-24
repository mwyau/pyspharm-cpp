<!-- SPDX-License-Identifier: GPL-2.0-or-later -->

# pyspharm-cpp

[![CI](https://github.com/mwyau/pyspharm-cpp/actions/workflows/ci.yml/badge.svg)](https://github.com/mwyau/pyspharm-cpp/actions/workflows/ci.yml)
[![Docs](https://readthedocs.org/projects/pyspharm-cpp/badge/?version=latest)](https://pyspharm-cpp.readthedocs.io/en/latest/?badge=latest)
[![Codecov](https://codecov.io/github/mwyau/pyspharm-cpp/graph/badge.svg)](https://codecov.io/github/mwyau/pyspharm-cpp)
[![PyPI version](https://img.shields.io/pypi/v/pyspharm-cpp)](https://pypi.org/project/pyspharm-cpp/)
[![Conda version](https://img.shields.io/conda/vn/conda-forge/pyspharm-cpp)](https://anaconda.org/channels/conda-forge/packages/pyspharm-cpp/overview)
[![License](https://img.shields.io/github/license/mwyau/pyspharm-cpp)](https://github.com/mwyau/pyspharm-cpp/blob/main/LICENSE)

`pyspharm-cpp` implements the PySpharm Python API for scalar and vector spherical harmonic transforms (SHTs), spectral regridding, Legendre functions, and interpolation. It uses DUCC's C++ transforms with nanobind and requires only NumPy at runtime.

For DUCC's native C++ or Python transform API, use [DUCC](https://gitlab.mpcdf.mpg.de/mtr/ducc) directly. For higher-level atmospheric-science operations, including vorticity, divergence, streamfunction, velocity potential, filtering, and regridding, use [spharmgrid](https://github.com/mwyau/spharmgrid). Spharmgrid supports Xarray fields with a DUCC CPU backend and offers separate JAX/S2FFT and PyTorch/torch-harmonics interfaces for GPU calculations.

## Install

Install from PyPI when a distribution is available for your platform:

```sh
python -m pip install pyspharm-cpp
```

Python 3.11 or newer is required. Source builds need a C++17 compiler. Import the Python module as `spharm_cpp`.

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

## Performance

![Measured SHT latencies on a Gaussian 1280 × 640 grid](benchmarks/timings-gl-640.svg)

Float32 timings for a Gaussian 1280 × 640 grid on a four-CPU GitHub Actions
runner with AVX-512 (`x86-64-v4`) active. Two warmups and seven timed
repetitions were used.

| Operation        |   PySpharm | pyspharm-cpp |     DUCC | Speedup vs PySpharm |
| ---------------- | ---------: | -----------: | -------: | ------------------: |
| Scalar analysis  | 184.648 ms |     3.450 ms | 3.199 ms |               53.5× |
| Scalar synthesis |  47.849 ms |     3.065 ms | 2.964 ms |               15.6× |
| Wind analysis    | 612.969 ms |     7.590 ms | 6.405 ms |               80.8× |

See [all five benchmark cases and reproduction instructions](benchmarks/README.md).

## License

GPL-2.0-or-later. DUCC is included as an unmodified submodule; see [third-party notices](THIRD_PARTY_NOTICES.md). DUCC is used under its respective license.
