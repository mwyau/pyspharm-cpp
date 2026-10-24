<!-- SPDX-License-Identifier: GPL-2.0-or-later -->

# Numerical conventions

## Grids

Grid arrays use shape `(nlat, nlon)` for one field and
`(nlat, nlon, nfields)` for multiple fields. Latitude rows are stored from
north to south. Longitude samples are equally spaced from zero eastward,
without a duplicate endpoint:

```{math}
\lambda_j = 360^\circ \frac{j}{n_{\mathrm{lon}}}, \qquad j=0,\ldots,n_{\mathrm{lon}}-1.
```

Regular grids use Clenshaw–Curtis latitudes including the poles. Gaussian grids
use Gauss–Legendre latitudes and do not include the poles. The
`gaussian_lats_wts(nlat)` function returns latitudes in degrees and quadrature
weights for integration over sine of latitude; the weights sum to 2 up to
floating-point rounding.

## Spectral arrays

Coefficients are packed in triangular order. Zonal wavenumber `m` is the outer
index, from 0 through `ntrunc`; degree `l` is the inner index, from `m` through
`ntrunc`. The number of coefficients is

```{math}
N_{\mathrm{coeff}} = \frac{(n_{\mathrm{trunc}}+1)(n_{\mathrm{trunc}}+2)}{2}.
```

For a coefficient pair `(m, l)`, the zero-based packed index is

```{math}
i(m,l) = \frac{m(2n_{\mathrm{trunc}}-m+3)}{2} + (l-m).
```

`getspecindx(ntrunc)` returns the `m` and `l` arrays for this layout. For
example, the first entries are `(m,l)=(0,0), (0,1), (0,2), (1,1), (1,2)`.

The associated Legendre functions use latitude
:math:`\phi` and exclude the Condon–Shortley phase:

```{math}
\bar P_{lm}(\phi) =
\sqrt{\frac{2l+1}{2}\frac{(l-m)!}{(l+m)!}}\,P_{lm}(\sin\phi).
```

Here :math:`P_{lm}` is the associated Legendre function without the
Condon–Shortley factor. The normalization gives
:math:`\bar P_{00}=1/\sqrt{2}`. For positive `m`, real-field synthesis
includes both conjugate longitude orders and therefore uses twice the real
part of the complex contribution. The `m=0` coefficients of a real field are
real.

`legendre(lat, ntrunc)` and `legendre_many(latitudes, ntrunc)` return these
values in the same triangular order. Latitudes are in degrees.

## Analysis and synthesis limits

Scalar analysis is constrained by both latitude and longitude sampling:

| Grid                    | Latitude limit       | Longitude limit             |
| ----------------------- | -------------------- | --------------------------- |
| Regular Clenshaw–Curtis | `ntrunc <= nlat - 2` | `ntrunc <= floor(nlon / 2)` |
| Gaussian                | `ntrunc <= nlat - 1` | `ntrunc <= floor(nlon / 2)` |

The default analysis truncation is the minimum of these limits. Explicit
requests beyond a limit raise `ValueError`. The regular-grid `nlat - 1`
analysis case is not supported because the highest degree is not uniquely
recoverable under DUCC's native Clenshaw–Curtis analysis. Pole samples vanish
for every order `m > 0`, leaving only `nlat - 2` interior latitude samples.
At `ntrunc=nlat-1`, order `m=1` has `nlat-1` associated-Legendre
coefficients but only `nlat-2` informative latitude samples. With colatitude
:math:`\theta=\pi/2-\phi`, the sampled null field satisfies

```{math}
\sin[(n_{\mathrm{lat}}-1)\theta_j]
= \sin\theta_j\,U_{n_{\mathrm{lat}}-2}(\cos\theta_j) = 0,
\qquad \theta_j = \frac{j\pi}{n_{\mathrm{lat}}-1},
```

where `U` denotes the Chebyshev polynomial of the second kind.
Orders `m=0` and `m>=2` can have full-rank latitude matrices when the
longitude samples resolve them. At `ntrunc=nlat-2`, the latitude matrices
are full column rank across the supported analysis band. Full coefficient recovery also
requires resolved longitude modes; an even-longitude Nyquist order has an
unobservable imaginary coefficient component. The implementation rejects the
whole regular-grid `nlat-1` analysis request rather than applying an
order-specific correction.

`getvrtdivspec` and `getpsichi` use DUCC spin-1 analysis, whose longitude
limit is `ntrunc <= floor((nlon - 1) / 2)`. Their default truncation uses this
limit, so even-longitude vector analysis stops below Nyquist. An explicit
request that includes Nyquist raises `ValueError`. `getuv` and `getgrad` can
synthesize Nyquist modes when the latitude synthesis limit permits them.

Synthesis accepts triangular spectra through degree `nlat - 1`. This
does not guarantee unique recovery of every coefficient from sampled values. At even `nlon`, the scalar
longitude Nyquist order `m=nlon/2` has real samples proportional to
`(-1)**j`; its imaginary component is identically zero on the sampled grid.
Scalar analysis applies the half normalization required when positive and
negative Fourier orders coincide at Nyquist. Imaginary Nyquist coefficients
cannot be recovered. Vector operations have additional spin-one Nyquist
degeneracy; direct comparisons use orders below Nyquist unless a result
specifically studies that boundary.

`regrid` uses an analysis limit valid for the input grid and a synthesis limit
valid for the output grid. Its default is the minimum of those limits.

## Precision and physical units

`Spharmt` defaults to `numpy.float64` maps and `numpy.complex128` spectra.
`dtype=numpy.float32` selects `float32` maps and `complex64` spectra. Input
arrays must match the selected precision. `gaussian_lats_wts` returns `float64` coordinates and
quadrature weights; `getspecindx` returns NumPy `intp` index arrays. `getgeodesicpts(m)` defaults to `float64`
coordinates calculated in double precision. Passing `dtype=numpy.float32`
selects the legacy-precision float32 calculation and output used for direct
PySpharm comparisons.

`rsphere` is positive and finite, in the length unit used by the caller. Scalar
analysis and synthesis do not apply a radius factor. `getgrad` divides scalar
derivatives by `rsphere`; wind vorticity and divergence use the corresponding
spherical derivative factors. Use consistent units for radius, wind and
potential fields.

## Vector component order

The public methods use the following order:

| Method                        | Inputs or outputs                          |
| ----------------------------- | ------------------------------------------ |
| `getgrad(chispec)`            | eastward gradient, northward gradient      |
| `getvrtdivspec(ugrid, vgrid)` | vorticity spectrum, divergence spectrum    |
| `getuv(vrtspec, divspec)`     | eastward wind, northward wind              |
| `getpsichi(ugrid, vgrid)`     | streamfunction map, velocity-potential map |

`specintrp(lon, dataspec, legfuncs)` evaluates a real field at a longitude
using one row of precomputed Legendre values. `specintrp_many` evaluates a
vector of longitudes and matching Legendre rows. Longitudes are in degrees and
may be wrapped by multiples of 360 degrees.
