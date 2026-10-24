<!-- SPDX-License-Identifier: GPL-2.0-or-later -->

# pyspharm-cpp

`pyspharm-cpp` computes scalar and vector spherical harmonic transforms on regular Clenshaw–Curtis and Gaussian latitude–longitude grids. It also implements gradients, vorticity, divergence, streamfunction, velocity potential, regridding, and spectral interpolation.

Its NumPy API follows PySpharm's coefficient ordering. DUCC computes the transforms; PySpharm/SPHEREPACK is used only for independent reference comparisons, not at runtime.

```{toctree}
---
maxdepth: 2
---
quickstart
api
conventions
pyspharm-compatibility
```

For development instructions, see [CONTRIBUTING.md](https://github.com/mwyau/pyspharm-cpp/blob/main/CONTRIBUTING.md). For the bundled DUCC license, see [third-party notices](https://github.com/mwyau/pyspharm-cpp/blob/main/THIRD_PARTY_NOTICES.md).
