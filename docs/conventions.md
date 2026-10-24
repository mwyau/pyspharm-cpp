<!-- SPDX-License-Identifier: GPL-2.0-or-later -->

# Numerical conventions

## Grid coordinates

Grid arrays have shape `(nlat, nlon)`, or `(nlat, nlon, nfields)` for multiple fields. Latitude rows run from north to south. Longitude increases eastward from zero without a repeated endpoint:

```{math}
\lambda_j = \frac{360^\circ j}{n_{\mathrm{lon}}},
\qquad j=0,\ldots,n_{\mathrm{lon}}-1.
```

Regular Clenshaw–Curtis grids include the poles. Gaussian grids use Gauss–Legendre latitudes, which exclude the poles. `gaussian_lats_wts(nlat)` returns latitudes in degrees and Gaussian quadrature weights for integration with respect to the sine of latitude; the weights sum to 2 up to rounding.

## Coefficients

Spectral coefficients are packed in triangular order: longitudinal order `m` is the outer index and spherical harmonic degree `l` varies from `m` to `ntrunc`. The coefficient count and zero-based index are

```{math}
N_{\mathrm{coeff}} = \frac{(n_{\mathrm{trunc}}+1)(n_{\mathrm{trunc}}+2)}{2},
\qquad
i(m,l) = \frac{m(2n_{\mathrm{trunc}}-m+3)}{2}+(l-m).
```

`getspecindx(ntrunc)` returns the order and degree arrays, `(m,l)`, in this order.

The associated Legendre functions exclude the Condon–Shortley phase and use the normalization

```{math}
\bar P_{lm}(\phi) =
\sqrt{\frac{2l+1}{2}\frac{(l-m)!}{(l+m)!}}\,
P_{lm}(\sin\phi),
```

where latitude `\phi` is in radians in this formula and `P_{lm}` excludes the Condon–Shortley phase. In particular, `\bar P_{00}=1/\sqrt{2}`. The public `legendre` and `legendre_many` functions accept latitudes in **degrees**.

For real fields, positive orders contribute twice the real part of the complex harmonic. At `m=0` the spectral coefficient is real.

## Truncation and recoverability

| Grid | Maximum scalar analysis truncation |
| --- | --- |
| Regular Clenshaw–Curtis | `min(nlat-2, nlon//2)` |
| Gaussian | `min(nlat-1, nlon//2)` |

The default analysis uses the applicable maximum. Larger explicit requests raise `ValueError`. Vector analysis (`getvrtdivspec` and `getpsichi`) also requires `ntrunc <= (nlon-1)//2` because DUCC's spin-1 analysis excludes the even-longitude Nyquist order.

The highest regular-grid degree is not uniquely recoverable for every longitudinal order. At `ntrunc=nlat-1`, the `m=1` latitude matrix has `nlat-1` coefficients but only `nlat-2` informative interior latitude samples. With colatitude `\theta=\pi/2-\phi`, one sampled null field is

```{math}
\sin[(n_{\mathrm{lat}}-1)\theta_j] =
\sin\theta_j\,U_{n_{\mathrm{lat}}-2}(\cos\theta_j) = 0,
\quad \theta_j = \frac{j\pi}{n_{\mathrm{lat}}-1},
```

where `U` denotes a Chebyshev polynomial of the second kind. Restricting regular-grid analysis to `nlat-2` avoids this deficiency for the supported latitude band.

Synthesis accepts triangular spectra through degree `nlat-1`; acceptance does not imply that analysis can uniquely recover every coefficient.

For even `nlon`, the scalar Nyquist order `m=nlon/2` samples `(-1)^j` in longitude. Its imaginary component is unobservable. Scalar analysis applies a factor of one half to its positive-order coefficient to match PySpharm's normalization. Spin-1 analysis does not accept this order.

## Precision, units, and component order

`Spharmt` uses float64 maps and complex128 coefficients by default; `dtype=np.float32` selects float32/complex64. Inputs must match the selected precision. `getspecindx` returns NumPy `intp` indices. `getgeodesicpts` returns float64 coordinates by default, with a float32 option.

`rsphere` is a positive radius in the user's length units. Scalar transforms are independent of radius; gradients divide by radius, and wind-derived vorticity and divergence use spherical derivatives. Component order is eastward (`u`) then northward (`v`).

`getuv` returns eastward and northward wind; `getvrtdivspec` returns vorticity then divergence spectra; `getpsichi` returns streamfunction then velocity potential. Longitude arguments for `specintrp` and `specintrp_many` are in degrees.
