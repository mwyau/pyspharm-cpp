<!-- SPDX-FileCopyrightText: 2026 Albert M. W. Yau -->
<!-- SPDX-License-Identifier: GPL-2.0-or-later -->

# Quickstart

## Installation

Install from PyPI:

```sh
python -m pip install pyspharm-cpp
```

Or install from conda-forge with Conda or Mamba:

```sh
conda install -c conda-forge pyspharm-cpp
# Or:
mamba install -c conda-forge pyspharm-cpp
```

The package imports as `spharm_cpp`; the examples use `import spharm_cpp as spharm`. Python 3.11 or newer and NumPy are required. For source builds, see [contributing](https://github.com/mwyau/pyspharm-cpp/blob/main/CONTRIBUTING.md).

## Scalar analysis and synthesis

```python
import numpy as np
import spharm_cpp as spharm

nlon, nlat = 32, 16
sht = spharm.Spharmt(nlon, nlat, gridtype="gaussian")
lat, _ = spharm.gaussian_lats_wts(nlat)
lon = np.arange(nlon) * 360.0 / nlon

field = np.cos(np.deg2rad(lat))[:, None] * np.cos(np.deg2rad(lon))[None, :]
spectrum = sht.grdtospec(field, ntrunc=1)
reconstructed = sht.spectogrd(spectrum)
np.testing.assert_allclose(reconstructed, field, rtol=2e-12, atol=2e-12)
```

Latitude rows run north to south; longitudes increase eastward without a duplicate endpoint. `gridtype="regular"` selects a Clenshaw–Curtis grid including both poles; `"gaussian"` selects Gaussian latitudes.

## Wind and derivative operations

`getgrad` returns eastward and northward gradient components. `getvrtdivspec` analyzes wind components into relative vorticity and divergence spectra; `getuv` synthesizes wind from those spectra.

```python
order, degree = spharm.getspecindx(1)
potential = np.zeros(degree.size, dtype=np.complex128)
potential[(order == 1) & (degree == 1)] = 1.0

eastward, northward = sht.getgrad(potential)
vorticity, divergence = sht.getvrtdivspec(eastward, northward, ntrunc=1)
u, v = sht.getuv(vorticity, divergence)
```

`getpsichi(u, v)` returns streamfunction and velocity-potential maps. Derivative and wind operations use the sphere radius `rsphere`, so the radius and input field units must be consistent.

## Precision and multiple fields

The default grid/coefficients are float64/complex128. Select float32/complex64 with `Spharmt(..., dtype=np.float32)`. Input precision must match the instance.

Multiple grid fields use shape `(nlat, nlon, nfields)`; spectral arrays use `(ncoeff, nfields)`. For example:

```python
sht32 = spharm.Spharmt(nlon, nlat, gridtype="gaussian", dtype=np.float32)
fields = np.stack((field, 2 * field), axis=-1).astype(np.float32)
spectra = sht32.grdtospec(fields, ntrunc=1)
assert spectra.shape == (3, 2)
reconstructed = sht32.spectogrd(spectra)
assert reconstructed.shape == fields.shape
```

## Regridding

`regrid` analyzes data on the source grid and synthesizes it on the target grid. The truncation must fit both grids.

```python
source = spharm.Spharmt(32, 16, gridtype="gaussian")
target = spharm.Spharmt(64, 32, gridtype="gaussian")
interpolated = spharm.regrid(source, target, field, ntrunc=1)
assert interpolated.shape == (32, 64)
```

See [numerical conventions](conventions.md) for truncation and normalization, [API reference](api.md) for all operations, and [PySpharm comparison](pyspharm-compatibility.md) for tested differences.
