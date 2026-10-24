<!-- SPDX-License-Identifier: GPL-2.0-or-later -->

# PySpharm compatibility inventory

The comparison covers the public Python functions and methods of
`spharm_cpp` and PySpharm. Tests check array shapes, coefficient order,
signs, and numerical values on grids and truncations supported by both
implementations. Approved differences are listed separately. Agreement in these
tests is not a measurement of accuracy for every grid or spectrum.

## Comparison implementations

| Component         | Version or revision                                                                                                                   | Role                                                                                   |
| ----------------- | ------------------------------------------------------------------------------------------------------------------------------------- | -------------------------------------------------------------------------------------- |
| `jswhit/pyspharm` | commit [`a464ebfed334188fc0ef275fad4aafd2d9fa716e`](https://github.com/jswhit/pyspharm/tree/a464ebfed334188fc0ef275fad4aafd2d9fa716e) | Pinned Python source used to audit the original public signatures.                     |
| `pyspharm-syl`    | distribution `1.2.1`; installed `spharm.__version__ == 1.0.9`                                                                         | Binary reference module and compiled `_spherepack` backend used by differential tests. |
| NumPy             | `1.26.4`                                                                                                                              | Locked validation environment compatible with the reference wheel.                     |
| `pyspharm-cpp`    | installed package version; implementation commit is visible in the CI run                                                             | Implementation under test.                                                             |
| DUCC              | pinned submodule revision printed in the parity job log                                                                               | C++ transform implementation used by `spharm_cpp`.                                     |

The installed `pyspharm-syl==1.2.1` Python source was compared with the pinned
upstream module. Its material differences are:

- It removes upstream `Spharmt.__repr__`.
- It replaces several `numpy.expand_dims` calls with equivalent `numpy.reshape`
  calls for rank adaptation. The validation tests cover those rank-one and
  rank-two paths.
- No other public behavior changes were observed in the Python source. The
  wheel's compiled SPHEREPACK backend is separately tested and is not assumed
  to be byte-for-byte equivalent to source at the pinned commit.

The parity CI job imports `spharm` and `spharm_cpp` as separate modules and
prints both native-module paths, package and NumPy versions, Python version,
DUCC revision, and CPU profile in the job log. PySpharm/SPHEREPACK is not a
runtime dependency and is installed only in this dedicated test workflow.

## Reading compatibility results

Each tested comparison is classified as:

- **PASS** — observed behavior matches within the operation-specific tolerance.
- **ACCEPTED_DIFFERENCE** — a behavior difference has an explicit approved
  explanation and a passing behavior test.
- **FAIL** — a numerical mismatch, regression or unapproved difference.
- **BLOCKED** — the reference or required test did not run.

Missing manifest entries, failed comparisons, and blocked reference imports
fail the parity job. An accepted difference has its own passing behavior test;
it is not skipped. The job displays individual pytest results and a contract
status summary directly in its log.

## Legacy public API inventory

The original and C++ signatures are shown together. `None` in a C++ default
means the extension selects its documented native default. C++ arrays must
match the selected `Spharmt.dtype`; unsupported types are rejected instead of
being silently converted.

### Constructor and instance operations

| Contract                                   | Original PySpharm signature                                                                                                                                                             | Current C++ signature                                                                                                                                                | Outputs and compatibility                                                                                                                                                                      | Evidence / status                                                                                                                                                                                                                                            |
| ------------------------------------------ | --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | -------------------------------------------------------------------------------------------------------------------------------------------------------------------- | ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ |
| `Spharmt`                                  | `Spharmt(nlon, nlat, rsphere=6.3712e6, gridtype='regular', legfunc='stored')`                                                                                                           | `Spharmt(nlon, nlat, rsphere=6371200.0, gridtype='regular', legfunc='computed', *, dtype=None)`                                                                      | Same grid, radius and validation concepts. C++ defaults to float64 and only accepts `legfunc='computed'`. `dtype=np.float32` selects float32/complex64.                                        | `test_constructor_defaults_signature_attributes_and_storage_policy`; precision and stored-mode differences are explicitly approved.                                                                                                                          |
| `Spharmt` attributes                       | `nlon`, `nlat`, `rsphere`, `gridtype`, `legfunc`; core values cannot be rebound or deleted; arbitrary user attributes are allowed.                                                      | Same five read-only core values, plus read-only `dtype`; arbitrary user attributes are enabled by nanobind dynamic attributes.                                       | Compatible for the tested public state. Work arrays such as `wshaes` are backend workspace details and are not exposed by C++.                                                                 | `test_constructor_defaults_signature_attributes_and_storage_policy`.                                                                                                                                                                                         |
| `Spharmt.__repr__`                         | Pinned upstream: `Spharmt(nlon, nlat, rsphere, gridtype, legfunc)`. Removed by `pyspharm-syl==1.2.1`.                                                                                   | Rich diagnostic representation including grid and dtype.                                                                                                             | Intentional diagnostic representation difference; the legacy wheel's default object representation is also recorded.                                                                           | `test_representation_is_diagnostic_and_class_names_remain_distinct`; approved.                                                                                                                                                                               |
| `grdtospec(datagrid, ntrunc=None)`         | Rank-2 `(nlat,nlon)` float32 map gives rank-1 packed spectrum; rank-3 `(nlat,nlon,nfields)` gives rank-2 `(ncoeff,nfields)`. Default truncation is `nlat-1`.                            | Same ranks and coefficient order; output is complex64 or complex128 according to the instance. Default truncation is bounded by the grid's supported analysis limit. | Values agree on the shared band-limited domain. Regular-grid analysis intentionally stops at `min(nlat-2, floor(nlon/2))`; Gaussian analysis stops at `min(nlat-1, floor(nlon/2))`.            | `test_scalar_analysis_synthesis_parity_across_grid_parity`, `test_grid_analysis_handles_readonly_and_strided_maps_without_mutation`, `test_zero_truncation_rank_and_batch_axis_are_consistent`, maximum-degree and Nyquist investigations.                   |
| `spectogrd(dataspec)`                      | Rank-1 spectrum gives rank-2 map; rank-2 `(ncoeff,nfields)` gives rank-3 map. The map is float32. The packed coefficient count determines `ntrunc`; documented upper bound is `nlat-1`. | Same axis ordering and shapes; map dtype follows the instance. Synthesis supports through `nlat-1` where longitude sampling permits.                                 | Matches on the shared supported domain and against independent analytical basis samples. At the regular maximum degree the legacy backend has an additional omission described below.          | `test_scalar_analysis_synthesis_parity_across_grid_parity`, `test_float64_synthesis_matches_independent_decimal_basis_samples`, `test_spectral_synthesis_handles_readonly_strided_coefficients_without_mutation`, maximum-degree and Nyquist investigations. |
| `getvrtdivspec(ugrid, vgrid, ntrunc=None)` | Same-shaped rank-2 or rank-3 float32 eastward/northward winds; returns vorticity and divergence spectra with a matching field axis. Default `ntrunc=nlat-1`.                            | Same component order and instance precision. The default respects the grid's latitude analysis limit and DUCC's spin-1 limit `floor((nlon-1)/2)`.                    | Directly compared with SPHEREPACK for resolved modes below the vector Nyquist boundary. Degree zero is constrained to zero vector content.                                                     | `test_vector_operations_match_independent_spherepack`; regular analysis restriction is applied consistently.                                                                                                                                                 |
| `getuv(vrtspec, divspec)`                  | Same-shaped rank-1 or rank-2 vorticity/divergence spectra; returns eastward and northward maps, rank-2 or rank-3.                                                                       | Same component order and batch-axis placement, in instance precision.                                                                                                | Direct synthesis comparison verifies geographical signs and radius scaling. Even-longitude spin-one Nyquist modes have their own degeneracy and are outside the common comparison band.        | `test_vector_operations_match_independent_spherepack`.                                                                                                                                                                                                       |
| `getgrad(chispec)`                         | Rank-1 or rank-2 scalar-potential spectrum; returns eastward and northward gradient maps, rank-2 or rank-3.                                                                             | Same component order and shapes, in instance precision.                                                                                                              | Directly compared with the independent reference; derivative scaling uses `rsphere`.                                                                                                           | `test_vector_operations_match_independent_spherepack`; `test_float64_synthesis_matches_independent_decimal_basis_samples`.                                                                                                                                   |
| `getpsichi(ugrid, vgrid, ntrunc=None)`     | Same-shaped rank-2 or rank-3 winds; returns streamfunction and velocity-potential maps. Default `ntrunc=nlat-1`.                                                                        | Same map order and batch axis, bounded by the vector-analysis limit `floor((nlon-1)/2)` and the grid's latitude limit. At `ntrunc=0`, returns two zero maps.         | PySpharm-syl 1.2.1 raises `IndexError` at `ntrunc=0`; the C++ zero-map convention follows the constant-scalar differential convention. This maintainer-approved behavior difference is tested. | `test_vector_operations_match_independent_spherepack`; `test_degree_zero_velocity_potentials_follow_constant_scalar_convention`.                                                                                                                             |
| `specsmooth(datagrid, smooth)`             | Rank-2 or rank-3 float32 map and rank-1 smoothing factors of length `nlat`; returns a map of the same input rank. Legacy implementation analyzes at `nlat-1`.                           | Same rank and field-axis order; selected precision and supported analysis truncation apply.                                                                          | Compared with an explicit band-limited field and filter. The regular-grid default cannot request a degree that the C++ analysis intentionally rejects.                                         | `test_spectral_smoothing_matches_on_an_explicit_band_limited_field`.                                                                                                                                                                                         |

### Legacy module functions

| Contract            | Original PySpharm signature and behavior                                                                                                                                                                                                                            | Current C++ signature and behavior                                                                                                                                                                                                               | Compatibility, discrepancy and status                                                                                                                                                                                                    | Exact evidence                                                                                                                             |
| ------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ | ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------ |
| `regrid`            | `regrid(grdin, grdout, datagrid, ntrunc=None, smooth=None)`. Default truncation is `min(grdin.nlat-1, grdout.nlat-1)`; optional filter has length `grdout.nlat`.                                                                                                    | Same arguments. Uses source-grid analysis and target-grid synthesis bounds; optional filter; native input/output precision conversion is supported.                                                                                              | Direct parity uses a shared explicit truncation, both grid-type directions, two fields and an explicit filter. The higher legacy default is constrained by the approved analysis policy.                                                 | `test_regrid_matches_with_explicit_shared_band_and_optional_filter`.                                                                       |
| `gaussian_lats_wts` | `gaussian_lats_wts(nlat) -> (lats, weights)`, rank-1 float64 arrays.                                                                                                                                                                                                | Same signature, rank and float64 dtype.                                                                                                                                                                                                          | Compatible. Weights integrate over sine of latitude and sum to 2 within rounding.                                                                                                                                                        | `test_gaussian_latitudes_and_weights_match_spherepack`.                                                                                    |
| `getspecindx`       | `getspecindx(ntrunc) -> (m, l)`. Runtime dtype comes from `numpy.indices` (normally int64/`np.intp` on supported 64-bit NumPy), despite its historical Int32 docstring. Ordering is degree inner, zonal order outer. `numpy.squeeze` returns scalars at `ntrunc=0`. | Same signature; returns rank-1 `np.intp` arrays at every truncation, including zero.                                                                                                                                                             | Ordering and values match. Pointer-width integer dtype follows NumPy's native index type; rank-one zero-truncation output is an approved shape correction.                                                                               | `test_spectral_indices_match_numpy_native_reference_order`.                                                                                |
| `getgeodesicpts`    | `getgeodesicpts(m) -> (lat, lon)`, ordered rank-1 float32 arrays of `10*(m-1)**2+2` points.                                                                                                                                                                         | `getgeodesicpts(m, *, dtype=None)`, ordered rank-1 float64 arrays by default; `dtype=np.float32` requests the legacy-precision calculation and output.                                                                                           | Float32 outputs are compared in returned order, allowing only circular equivalence at the 0/360-degree seam. Float64 geometry is a newly approved precision improvement over the legacy float32 result.                                  | `test_geodesic_float32_order_and_double_precision_geometry`.                                                                               |
| `legendre`          | `legendre(lat, ntrunc) -> rank-1 float32 normalized associated Legendre values`.                                                                                                                                                                                    | `legendre(lat, ntrunc, *, dtype=None)`, float64 by default; `dtype=np.float32` gives float32.                                                                                                                                                    | Float32 matches the reference over the tested domain. Default/high-degree float64 is checked against an independent high-precision Decimal recurrence. Phase and triangular order match.                                                 | `test_legendre_float32_matches_reference_and_float64_matches_decimal`.                                                                     |
| `specintrp`         | `specintrp(lon, dataspec, legfuncs) -> Python float`; the spectrum and Legendre vector must imply the same triangular truncation. Longitude is in degrees.                                                                                                          | Same arguments and scalar return type; float32/complex64 or float64/complex128 inputs are supported when paired consistently. Longitude is reduced modulo 360 before radians conversion, preserving phase accuracy for very large finite values. | Direct reference comparison covers poles, equator, interior latitudes, negative longitude and longitude wrapping. Coefficient phase and order match. A separate Decimal argument-reduction regression covers longitudes through `1e308`. | `test_scalar_spectral_interpolation_matches_reference_at_lon_and_lat_points`, `test_scalar_interpolation_reduces_large_finite_longitudes`. |

The public contract manifest is `tests/parity/manifest.py`. A collection hook
fails if any listed contract or approved difference has no marked behavior
test.

## Additive APIs and precision property

These APIs do not replace or alter legacy numerical operations:

| API              | Signature and output                                                                                                                                                                                  | Independent validation                                                                                                                                                                                                             |
| ---------------- | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `legendre_many`  | `legendre_many(latitudes, ntrunc, *, dtype=None)`. Rank-1 latitude input produces `(nlatitudes,ncoeff)`; float64 by default and optional float32.                                                     | Each row is compared with repeated `pyspharm-syl` scalar `legendre` calls; a separate test compares it with repeated C++ scalar calls.                                                                                             |
| `specintrp_many` | `specintrp_many(longitudes, dataspec, legfuncs)`. Rank-1 float64 values in input longitude order; accepts one matching Legendre row per longitude.                                                    | Every element is compared with repeated reference `specintrp` calls.                                                                                                                                                               |
| `cpu_info`       | `cpu_info() -> dict[str, object]`. Returns DUCC's `architecture`, `features`, and `multiarch` entries; multiarch builds also return compiled/available profiles, configured limit and active profile. | The implementation follows the vendored upstream `ducc0.misc.cpu_info()` source. The test checks the exact key set and Python value types for the active profile. No DUCC version or source hash is included in this return value. |
| `Spharmt.dtype`  | Read-only string, `"float32"` or `"float64"`.                                                                                                                                                         | Test verifies analysis and synthesis dtypes for each selection.                                                                                                                                                                    |

## Observed behavior, mathematical limits and API choices

### Truncation and regular-grid maximum degree

The supported C++ analysis bounds are:

```text
regular:  lmax <= min(nlat-2, floor(nlon/2))
Gaussian: lmax <= min(nlat-1, floor(nlon/2))
```

The regular-grid constraint is based on observability, not a tolerance choice.
For `lmax=nlat-1`, every `m>0` vanishes at the two poles, leaving `nlat-2`
interior latitude samples. At `m=1`, there are `nlat-1` associated-Legendre
columns and only `nlat-2` rows, with one null direction. One representation of
that null field is
`sin((nlat-1)*theta) = sin(theta)*U_(nlat-2)(cos(theta))`; for `m>=2`, the
latitude columns are full rank when their longitude modes are resolved. At
`lmax=nlat-2`, the latitude matrices are full column rank across the analysis
band. Full coefficient recovery also requires resolved longitude modes; at an
even-longitude Nyquist order, the imaginary coefficient component is not
observable. The C++ analysis rejects the whole regular `nlat-1` request rather
than making order-specific promises.

Tests for even and odd `nlat` construct the `m=1` null vector,
compare SPHEREPACK and DUCC synthesis at the maximum regular-grid degree, and
check analysis at the lower supported truncation. Gaussian tests cover
`lmax=nlat-1` with resolved longitude modes, including individual
highest-degree harmonics and mixed spectra.

**Accepted difference:** for regular-grid synthesis at `l=nlat-1`, the
tested SPHEREPACK `SHSES`/`SHSEC` implementation omits odd-`m`
contributions, although its Python wrapper accepts the coefficients.
`spharm_cpp` synthesizes the corresponding analytical harmonics. The tests
compare both outputs with an independent Legendre recurrence and verify
SPHEREPACK's even-/odd-order behavior. The difference concerns synthesis at
the highest regular-grid degree; it does not alter the C++ analysis limit.

Synthesis acceptance does not imply an invertible analysis for every
coefficient vector. For regular grids, `getvrtdivspec` and `getpsichi` use the
`nlat-2` latitude limit and DUCC's spin-1 longitude limit. `specsmooth` and
`regrid` use the scalar analysis limits. Regridding and smoothing comparisons
therefore use an explicit shared band. The legacy regrid default `min(nlat-1)`
is not treated as a same-band parity comparison.

### Longitude Nyquist

With even `nlon`, scalar order `m=nlon/2` samples the real mode as
`2*Pbar*(-1)**j`. The imaginary component samples as zero because
`sin(pi*j)=0`; it cannot be recovered. Direct comparison found that DUCC's
analysis coefficients at the positive Nyquist order require a factor of one
half to match SPHEREPACK's positive/negative-frequency convention. The C++
scalar wrapper applies that narrow correction; synthesis remains unchanged.
The test covers regular and Gaussian grids. It compares C++ and SPHEREPACK
float32 synthesis, both implementations' analysis of both grids, and samples
from an independent Decimal-based associated-Legendre recurrence. It also
checks the C++ float64 synthesis and recovery against that independent field.
A pure imaginary scalar Nyquist coefficient produces an all-zero sampled field
and all-zero recovered spectra; its imaginary component is unobservable.
DUCC spin-1 analysis excludes the even-longitude Nyquist order, so
`getvrtdivspec` and `getpsichi` reject truncations that include it. Vector
synthesis at and below Nyquist is compared with analytic gradient and wind
fields, while the SPHEREPACK parity cases cover orders below Nyquist.

### Shapes, ownership and validation

PySpharm uses `numpy.squeeze` in several paths. At `ntrunc=0`, this can turn a
single coefficient/index into a scalar. C++ consistently preserves the
rank-one coefficient axis. Rank-two maps, rank-three map batches, rank-one and
rank-two spectra, singleton field batches, read-only arrays and strided input
maps are covered. C++ outputs own their storage and inputs are not modified.

The approved zero-truncation shape correction applies to coefficient and index
arrays. For `getpsichi(..., ntrunc=0)`, PySpharm-syl 1.2.1 raises `IndexError`,
while C++ returns zero maps following the constant-scalar differential
convention. The maintainer approved this difference; the behavior test verifies
both results and records it as `ACCEPTED_DIFFERENCE`.

The C++ wrapper rejects unsupported precision, malformed coefficient counts,
negative or out-of-grid truncations, nonfinite parameters, invalid radius and
invalid array shapes explicitly. Error wording is not treated as numerical
parity. Core transform settings remain read-only; arbitrary Python instance
attributes remain supported.

### Precision, tolerances and accuracy limits

`Spharmt` defaults to float64 maps and complex128 spectra. Set
`dtype=np.float32` for float32 maps and complex64 spectra. Legendre functions
also default to float64 and support explicit float32. Gaussian latitudes and
weights remain float64. `getspecindx` returns `np.intp`. Geodesic coordinates
now default to float64 and are calculated in double precision; the optional
`dtype=np.float32` path matches the legacy precision and point order.

Tolerances are operation-specific and tied to observed precision and scale.
The float32 transform and vector comparisons use few-parts-per-million
relative bounds with absolute floors appropriate to each deterministic test
field. High-degree float32 Legendre comparison allows `1e-5` absolute error;
the largest measured residual through degree 32 is about `8.1e-6`. Ordered
float32 geodesic longitude comparisons are circular at the 0/360 seam and use
an absolute bound of `7e-5` degrees, close to a few float32 ULPs at the largest
longitudes tested. Float64 synthesis basis fields are compared to Decimal
Legendre recurrence samples at `5e-14` relative and absolute tolerance; the
degree-192 Legendre test checks selected coefficients against an independent
90-digit recurrence. No bitwise equality or dynamic tolerance inflation is
used. Failed numerical assertions print maximum absolute and meaningful
relative errors for the affected operation and precision.

## Parity CI workflow

The `.github/workflows/validation.yml` workflow runs through manual
`workflow_dispatch`. It uses Linux x86-64, CPython 3.12,
`pyspharm-syl==1.2.1`, NumPy 1.26.4, and the locked validation group with a
binary reference wheel. It builds the C++ package from the checked-out source and compares it with
the installed reference. The
[October 24, 2026 run](https://github.com/mwyau/pyspharm-cpp/actions/runs/40060787028)
at commit `cdbee944` passed 37 tests: 20 public API contracts, 11 approved
differences, no failures, and no blocked tests. Results are available in the
job log. Machine-readable report artifacts are not yet generated.

Local equivalent after syncing the locked validation group:

```sh
uv sync --locked --python 3.12 --group build --group test --group validation \
  --no-build-package pyspharm-syl
uv run --locked --python 3.12 --no-sync python -m pytest -v -s \
  -o addopts=-ra tests/parity
```
