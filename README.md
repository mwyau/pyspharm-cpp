# pyspharm-cpp

[![CI](https://img.shields.io/github/actions/workflow/status/mwyau/pyspharm-cpp/ci.yml?branch=main&label=CI)](https://github.com/mwyau/pyspharm-cpp/actions/workflows/ci.yml)
[![Docs](https://readthedocs.org/projects/pyspharm-cpp/badge/?version=latest)](https://pyspharm-cpp.readthedocs.io/en/latest/?badge=latest)
[![Codecov](https://codecov.io/github/mwyau/pyspharm-cpp/graph/badge.svg)](https://codecov.io/github/mwyau/pyspharm-cpp)
[![PyPI version](https://img.shields.io/pypi/v/pyspharm-cpp)](https://pypi.org/project/pyspharm-cpp/)
[![Conda version](https://img.shields.io/conda/vn/conda-forge/pyspharm-cpp)](https://anaconda.org/channels/conda-forge/packages/pyspharm-cpp/overview)
[![License](https://img.shields.io/github/license/mwyau/pyspharm-cpp)](https://github.com/mwyau/pyspharm-cpp/blob/main/LICENSE)

`pyspharm-cpp` implements a PySpharm-compatible Python API for scalar and vector spherical harmonic transforms (SHTs), spectral regridding, Legendre functions, and interpolation. It uses ducc0's C++ transforms with nanobind and requires only NumPy at runtime.

For lower-level spherical harmonic transforms, consider using [ducc0](https://pypi.org/project/ducc0/). For atmospheric science operations such as relative vorticity, divergence, streamfunction, velocity potential, spectral filtering, and regridding, consider using [spharmgrid](https://github.com/mwyau/spharmgrid). Spharmgrid supports Xarray fields with a ducc0 CPU backend and separate JAX/S2FFT and PyTorch/torch-harmonics interfaces for GPU calculations.

## Installation

Install from PyPI:

```sh
python -m pip install pyspharm-cpp
```

Or install from conda-forge with Conda or Mamba:

```sh
conda install -c conda-forge pyspharm-cpp
mamba install -c conda-forge pyspharm-cpp
```

Python 3.11 or newer is required. Source builds need a C++17 compiler; Git checkouts also need the ducc0 submodule. Source distributions include the required ducc0 files. The Python module imports as `spharm_cpp`; `import spharm_cpp as spharm` is convenient when adapting PySpharm code.

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

- **Backend:** ducc0 replaces the SPHEREPACK transforms.
- **Precision:** float64/complex128 by default; float32/complex64 can be selected with `dtype=np.float32`.
- **Legendre functions:** `legfunc="computed"` is the default. Stored Legendre functions (`legfunc="stored"`) are not supported.
- **Batched interpolation:** `legendre` and `specintrp` also accept one-dimensional coordinate arrays, using overloads of the original scalar functions.
- **Truncation:** Clenshaw–Curtis analysis supports degrees through `min(nlat-2, nlon//2)`; Gaussian analysis through `min(nlat-1, nlon//2)`. Vector analysis excludes the even-longitude Nyquist order.
- **Other differences:** zero truncation, optional float64 geodesic coordinates, index dtype, and mixed-precision regridding are documented in the [PySpharm comparison](https://github.com/mwyau/pyspharm-cpp/blob/main/docs/pyspharm-compatibility.md).

See [numerical conventions](https://github.com/mwyau/pyspharm-cpp/blob/main/docs/conventions.md) for coefficient ordering and normalization, and the [quickstart](https://github.com/mwyau/pyspharm-cpp/blob/main/docs/quickstart.md) and [API reference](https://github.com/mwyau/pyspharm-cpp/blob/main/docs/api.md) for usage.

## Performance

[![Measured SHT latencies on a Gaussian 1280 × 640 grid](https://raw.githubusercontent.com/mwyau/pyspharm-cpp/main/benchmarks/timings-gl-640.svg)](https://github.com/mwyau/pyspharm-cpp/blob/main/benchmarks/timings-gl-640.svg)

Median timings in milliseconds on the Gaussian 1280 × 640 grid, which has the largest measured speedup among the five cases (69.42× for wind analysis). ducc0 and `pyspharm-cpp` use eight threads; `pyspharm` is serial. All timings use float32 and the same truncation.

| Operation        | pyspharm f32 | pyspharm-cpp f32 | ducc0 f32 | Speedup vs pyspharm |
| ---------------- | -----------: | ---------------: | --------: | ------------------: |
| Scalar analysis  |   139.024 ms |         2.103 ms |  2.411 ms |              66.11× |
| Scalar synthesis |    32.081 ms |         2.227 ms |  2.304 ms |              14.41× |
| Wind analysis    |   339.515 ms |         4.890 ms |  3.389 ms |              69.42× |
| Wind synthesis   |    67.204 ms |         5.432 ms |  4.558 ms |              12.37× |

See [benchmark results and reproduction instructions](https://github.com/mwyau/pyspharm-cpp/blob/main/benchmarks/README.md) for all five grid cases.

## License

GPL-2.0-or-later. ducc0 is included as an unmodified submodule; see [third-party notices](https://github.com/mwyau/pyspharm-cpp/blob/main/THIRD_PARTY_NOTICES.md). ducc0 is used under its respective license.
