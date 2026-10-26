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

where latitude {math}`\phi` is in radians in this formula and {math}`P_{lm}` excludes the Condon–Shortley phase. In particular, {math}`\bar P_{00}=1/\sqrt{2}`. `legendre` accepts scalar or one-dimensional array latitudes in **degrees**.

For real fields, positive orders contribute twice the real part of the complex harmonic. At `m=0` the spectral coefficient is real.

## Truncation and recoverability

| Grid                    | Maximum scalar analysis truncation |
| ----------------------- | ---------------------------------- |
| Regular Clenshaw–Curtis | `min(nlat-2, nlon//2)`             |
| Gaussian                | `min(nlat-1, nlon//2)`             |

The default analysis uses the applicable maximum. Larger explicit requests raise `ValueError`. Vector analysis (`getvrtdivspec` and `getpsichi`) also requires `ntrunc <= (nlon-1)//2` because ducc0's spin-1 analysis excludes the even-longitude Nyquist order.

The highest regular-grid degree is not uniquely recoverable for every longitudinal order. At `ntrunc=nlat-1`, the `m=1` latitude matrix has `nlat-1` coefficients but only `nlat-2` informative interior latitude samples. With colatitude {math}`\theta=\pi/2-\phi`, one sampled null field is

```{math}
\sin[(n_{\mathrm{lat}}-1)\theta_j] =
\sin\theta_j\,U_{n_{\mathrm{lat}}-2}(\cos\theta_j) = 0,
\quad \theta_j = \frac{j\pi}{n_{\mathrm{lat}}-1},
```

where `U` denotes a Chebyshev polynomial of the second kind. Restricting regular-grid analysis to `nlat-2` avoids this deficiency for the supported latitude band.

Synthesis accepts triangular spectra through degree `nlat-1`; acceptance does not imply that analysis can uniquely recover every coefficient.

For even `nlon`, the scalar Nyquist order `m=nlon/2` samples `(-1)^j` in longitude. Its imaginary component is unobservable. Scalar analysis applies a factor of one half to its positive-order coefficient to match PySpharm's normalization. Spin-1 analysis uses longitudinal orders below Nyquist.

## Precision, units, and component order

`Spharmt` uses float64 maps and complex128 coefficients by default; `dtype=np.float32` selects float32/complex64. Inputs must match the selected precision. `getspecindx` returns NumPy `intp` indices. `getgeodesicpts` computes its geometry and returns coordinates in the selected precision, float32 by default (matching PySpharm), or float64 with `dtype=np.float64`.

`legendre` accepts one latitude or a one-dimensional latitude array and returns packed associated Legendre values. Latitude conversion and trigonometric setup use binary64 before converting once to the selected recurrence type. The recurrence itself uses binary32 for `dtype=np.float32` and binary64 for `dtype=np.float64`. For identical finite latitudes and truncation, each batched row is bitwise identical to the scalar result at the same dtype and build profile.

`specintrp` accepts scalar longitude and rank-one Legendre values for a Python `float` result, or one-dimensional longitude and rank-two Legendre arrays for batched interpolation. Coefficient accumulation uses binary32 for complex64 coefficients and binary64 for complex128 coefficients. Longitude reduction and phase recurrence use binary64, then phases are converted to the coefficient precision before accumulation. The phase recurrence is reset every 32 orders to limit accumulated phase error. The scalar result is returned as a Python `float` after native-precision accumulation. Batched output follows the coefficient precision regardless of longitude-array dtype. For identical finite longitudes, coefficients, and Legendre values, scalar and batched results are bitwise identical at the same dtype and build profile.

Float32 recurrence and summation have more rounding than float64 recurrence and summation. Independent tests compare high-degree Legendre values with Decimal references and check the degree norm; these checks cover their specified inputs and do not bound all latitudes or spectra. For calculations that require additional precision at high spectral truncation, use float64.

`rsphere` is a positive radius in the user's length units. Scalar transforms are independent of radius; gradients divide by radius, and wind-derived vorticity and divergence use spherical derivatives. Component order is eastward (`u`) then northward (`v`).

`getuv` returns eastward and northward wind; `getvrtdivspec` returns vorticity then divergence spectra; `getpsichi` returns streamfunction then velocity potential. Longitude arguments for `specintrp` are in degrees.

## Thread selection

Transforms and batched numerical kernels use DUCC's thread pool. Set
`DUCC0_NUM_THREADS` to select the thread count; when it is unset,
`OMP_NUM_THREADS` is used. If both are unset, operations use one thread. A
thread count of zero asks DUCC to use all threads available to its pool, and a
positive count is limited to the pool's capacity. The setting is read on the
first numerical call in each process, so set it before that call.
