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

`Spharmt` uses float64 maps and complex128 coefficients by default; `dtype=np.float32` selects float32/complex64. Inputs must match the selected precision. `getspecindx` returns NumPy `intp` indices. `getgeodesicpts` calculates coordinates in float64 and optionally returns float32 values.

For array latitudes, `legendre` evaluates the associated Legendre recurrence at the requested precision and returns shape `(npoints, ncoeff)`. `specintrp` accepts scalar longitude and rank-one Legendre values for a Python `float` result, or one-dimensional longitude and rank-two Legendre arrays for batched interpolation. Batched interpolation accepts float32 Legendre values with complex64 coefficients or float64 Legendre values with complex128 coefficients; its output follows that pair regardless of the longitude-array dtype.

At high spectral truncation, float32 recurrence and summation roundoff produce larger differences from float64-intermediate calculations. In tests, batched `legendre` differed from the scalar float32 result by at most 2.6e-2 absolute at `ntrunc=1024` near 89.999° latitude, where the coefficient magnitude was about 29.5. At `ntrunc=2048`, the degree norm differed by at most 3.5e-4 relative across eight tested latitudes. These results describe the tested inputs and do not bound other latitudes or spectra. For high-degree calculations sensitive to float32 roundoff, prefer float64.

`rsphere` is a positive radius in the user's length units. Scalar transforms are independent of radius; gradients divide by radius, and wind-derived vorticity and divergence use spherical derivatives. Component order is eastward (`u`) then northward (`v`).

`getuv` returns eastward and northward wind; `getvrtdivspec` returns vorticity then divergence spectra; `getpsichi` returns streamfunction then velocity potential. Longitude arguments for `specintrp` are in degrees.

## Thread selection

Transforms and batched numerical kernels use DUCC's thread pool. Set
`DUCC0_NUM_THREADS` to select the thread count; when it is unset,
`OMP_NUM_THREADS` is used. If both are unset, operations use one thread. A
thread count of zero asks DUCC to use all threads available to its pool, and a
positive count is limited to the pool's capacity. The setting is read on the
first numerical call in each process, so set it before that call.
