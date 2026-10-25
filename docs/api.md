# API reference

The API is generated from the extension's docstrings. Functions and methods
accept NumPy arrays and use PySpharm's triangular coefficient ordering.

## Spherical harmonic transforms

```{eval-rst}
.. autoclass:: spharm_cpp.Spharmt
   :members: grdtospec, spectogrd, getvrtdivspec, getuv, getgrad, getpsichi, specsmooth
```

## Grids and coefficients

```{eval-rst}
.. autofunction:: spharm_cpp.gaussian_lats_wts

.. autofunction:: spharm_cpp.getspecindx

.. autofunction:: spharm_cpp.getgeodesicpts
```

## Spectral interpolation and regridding

`legendre(lat, ntrunc)` returns a one-dimensional array of normalized Legendre
values. With a one-dimensional latitude array, `legendre(latitudes, ntrunc)`
returns shape `(npoints, ncoeff)`.

`specintrp(lon, dataspec, legfuncs)` returns a Python `float` for a scalar
longitude and one-dimensional Legendre array. With a one-dimensional longitude
array and Legendre values of shape `(npoints, ncoeff)`, it returns a real array
of shape `(npoints,)`. Batched output is float32 for complex64 coefficients
and float64 for complex128 coefficients; Legendre precision must match.

```{eval-rst}
.. autofunction:: spharm_cpp.legendre

.. autofunction:: spharm_cpp.specintrp

.. autofunction:: spharm_cpp.regrid
```

## CPU information

```{eval-rst}
.. autofunction:: spharm_cpp.cpu_info
```

Grid shapes, normalization, and truncation are defined in
[numerical conventions](conventions.md). The [PySpharm comparison](pyspharm-compatibility.md)
describes tested differences.
