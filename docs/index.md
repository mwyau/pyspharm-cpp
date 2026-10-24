<!-- SPDX-License-Identifier: GPL-2.0-or-later -->

# pyspharm-cpp

`pyspharm-cpp` computes scalar and vector spherical harmonic transforms,
gradients, vorticity, divergence, streamfunction, and velocity potential on
regular or Gaussian latitude–longitude grids. Its NumPy interface follows
PySpharm's function names and coefficient layout, with documented differences
in precision, truncation, and Legendre evaluation. DUCC performs the transforms;
the Python package does not require SPHEREPACK or the `ducc0` Python package.

The [quickstart](quickstart.md) shows scalar and vector operations.
The [API reference](api.md) lists functions and methods, the
[numerical conventions](conventions.md) define coordinates and spectral
limits, and the [PySpharm comparison](pyspharm-compatibility.md) records tests
against PySpharm/SPHEREPACK.

```{toctree}
---
maxdepth: 2
---
quickstart
api
conventions
pyspharm-compatibility
```

For development and contribution instructions, see
[CONTRIBUTING.md](https://github.com/mwyau/pyspharm-cpp/blob/main/CONTRIBUTING.md).
For licensing and bundled DUCC sources, see
[THIRD_PARTY_NOTICES.md](https://github.com/mwyau/pyspharm-cpp/blob/main/THIRD_PARTY_NOTICES.md).
