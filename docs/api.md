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
