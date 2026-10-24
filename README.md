# pyspharm-cpp

[![CI](https://img.shields.io/github/actions/workflow/status/mwyau/pyspharm-cpp/ci.yml?branch=main&label=CI)](https://github.com/mwyau/pyspharm-cpp/actions/workflows/ci.yml)
[![Docs](https://readthedocs.org/projects/pyspharm-cpp/badge/?version=latest)](https://pyspharm-cpp.readthedocs.io/en/latest/?badge=latest)
[![Codecov](https://codecov.io/github/mwyau/pyspharm-cpp/graph/badge.svg)](https://codecov.io/github/mwyau/pyspharm-cpp)
[![PyPI version](https://img.shields.io/pypi/v/pyspharm-cpp)](https://pypi.org/project/pyspharm-cpp/)
[![Conda version](https://img.shields.io/conda/vn/conda-forge/pyspharm-cpp)](https://anaconda.org/channels/conda-forge/packages/pyspharm-cpp/overview)
[![License](https://img.shields.io/github/license/mwyau/pyspharm-cpp)](https://github.com/mwyau/pyspharm-cpp/blob/main/LICENSE)

`pyspharm-cpp` implements the PySpharm Python API for scalar and vector spherical harmonic transforms (SHTs), spectral regridding, Legendre functions, and interpolation. It uses ducc0's C++ transforms with nanobind and requires only NumPy at runtime.

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

Python 3.11 or newer is required. Source builds need a C++17 compiler and the ducc0 submodule. The Python module imports as `spharm_cpp`; `import spharm_cpp as spharm` is convenient when adapting PySpharm code.

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
- **Truncation:** Clenshaw–Curtis analysis supports degrees through `min(nlat-2, nlon//2)`; Gaussian analysis through `min(nlat-1, nlon//2)`. Vector analysis excludes the even-longitude Nyquist order.
- **Other differences:** zero truncation, geodesic coordinate precision, index dtype, and mixed-precision regridding are documented in the [PySpharm comparison](docs/pyspharm-compatibility.md).

See [numerical conventions](docs/conventions.md) for coefficient ordering and normalization, and the [quickstart](docs/quickstart.md) and [API reference](docs/api.md) for usage.

## Performance

![Measured SHT latencies on a Gaussian 1280 × 640 grid](benchmarks/timings-gl-640.svg)

Median timings in milliseconds on a Gaussian 1280 × 640 grid, with four ducc0 threads and serial PySpharm. All timings use float32 and the same truncation.

| Operation        | PySpharm f32 | pyspharm-cpp f32 | ducc0 f32 | Speedup vs PySpharm |
| ---------------- | -----------: | ---------------: | --------: | ------------------: |
| Scalar analysis  |   142.171 ms |         3.361 ms |  2.959 ms |              42.30× |
| Scalar synthesis |    45.133 ms |         2.908 ms |  2.901 ms |              15.52× |
| Wind analysis    |   349.341 ms |         7.901 ms |  6.260 ms |              44.21× |
| Wind synthesis   |    74.312 ms |         7.571 ms |  6.445 ms |               9.81× |

See [benchmark results and reproduction instructions](benchmarks/README.md) for all five grid cases.

## License

GPL-2.0-or-later. ducc0 is included as an unmodified submodule; see [third-party notices](THIRD_PARTY_NOTICES.md). ducc0 is used under its respective license.
