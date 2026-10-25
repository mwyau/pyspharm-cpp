# Comparison with PySpharm

`spharm_cpp` provides a PySpharm-compatible Python API using ducc0 transforms. The implementations agree in the tested cases, with the differences described below.

The PySpharm comparison uses `pyspharm-syl==1.2.1` (module version 1.0.9 and a compiled SPHEREPACK backend), NumPy 1.26.4, and CPython 3.12 on Linux x86-64. The original Python interface was compared with [`jswhit/pyspharm` at `a464ebf`](https://github.com/jswhit/pyspharm/tree/a464ebfed334188fc0ef275fad4aafd2d9fa716e).

## Tested operations

The tests cover the following interfaces, including batched inputs where applicable.

| Operation                                            | Comparison                                                                                   |
| ---------------------------------------------------- | -------------------------------------------------------------------------------------------- |
| `Spharmt`, attributes, `dtype`                       | Constructor, parameters, and output precision                                                |
| `grdtospec`, `spectogrd`                             | Scalar synthesis/analysis on regular and Gaussian grids; coefficient order; boundary degrees |
| `getvrtdivspec`, `getuv`, `getgrad`, `getpsichi`     | Eastward/northward signs, derivative scaling, vorticity, divergence, and potentials          |
| `specsmooth`, `regrid`                               | Smoothing and regridding at an explicit shared truncation                                    |
| `gaussian_lats_wts`, `getspecindx`, `getgeodesicpts` | Quadrature, index ordering, and geodesic point coordinates                                   |
| `legendre`, `specintrp`                              | Normalization, phase, polar behavior, and longitude                                          |
| `legendre_many`, `specintrp_many`                    | Batched results against scalar operations using precision-specific tolerances                |
| `cpu_info`                                           | CPU metadata keys and types                                                                  |

Direct PySpharm comparisons cover the shared API. The additional `legendre_many` and `specintrp_many` functions are compared with scalar calls; `cpu_info` is checked for metadata keys and types.

## Differences

| Behavior                          | `spharm_cpp`                                                                                 |
| --------------------------------- | -------------------------------------------------------------------------------------------- |
| Module and backend                | `spharm_cpp` uses ducc0 for spherical harmonic transforms.                                   |
| Legendre evaluation               | Stored Legendre functions (`legfunc="stored"`) are not supported.                            |
| Precision                         | float64/complex128 by default, float32/complex64 when requested.                             |
| Regular-grid analysis             | Maximum degree `nlat-2`, also bounded by longitude sampling; PySpharm accepts `nlat-1`.      |
| Scalar outputs at truncation zero | Coefficient and index arrays keep an explicit length-one axis instead of becoming scalars.   |
| Geodesic coordinates              | Calculated in float64; `dtype=np.float32` casts the output coordinates.                      |
| Regridding                        | Supports conversion between the two native precision modes.                                  |
| Input validation                  | Invalid shapes, dtypes, parameter values, and unsupported truncations raise explicit errors. |
| `getpsichi` at `ntrunc=0`         | Returns two zero maps; the reference raises `IndexError`.                                    |

See [numerical conventions](conventions.md) for the supported grid and vector-analysis limits.

## Regular-grid truncation

`spharm_cpp` limits Clenshaw–Curtis analysis to degree `nlat-2`, also subject to the longitude sampling limit. PySpharm accepts `nlat-1`. On the regular grid, not every coefficient at that highest degree can be recovered uniquely. Synthesis and analysis have separate degree limits.

At `ntrunc=nlat-1`, the tested SPHEREPACK regular-grid synthesis path omits highest-degree coefficients at odd longitudinal orders; `spharm_cpp` evaluates these modes. This difference is tested separately from analysis recoverability.

## Longitude Nyquist

For even `nlon`, the scalar order `m=nlon/2` samples {math}`e^{im\lambda_j}=(-1)^j`. Its imaginary coefficient component cannot be determined from the samples. The scalar analysis conversion applies a half-factor at this order, matching PySpharm's positive/negative-frequency convention in direct comparisons. Analytic basis tests separately check the sampled field.

Vector analysis excludes the even-longitude Nyquist order; `getvrtdivspec` and `getpsichi` use the upper longitude limit `(nlon-1)//2`. Tests of vector synthesis include analytic gradient and wind fields at the boundary, while direct SPHEREPACK vector comparisons use resolved orders below Nyquist.

## Precision and validation evidence

Float32 numerical comparisons use operation-specific absolute and relative tolerances, rather than bitwise identity. Selected float64 Legendre functions are checked against an independent high-precision recurrence; basis synthesis uses Decimal reference samples with a `5e-14` tolerance. Geodesic float32 longitude comparisons account for the 0°/360° seam and test the rounded float64 result.

Run the PySpharm comparison with the [`validation.yml` workflow](https://github.com/mwyau/pyspharm-cpp/actions/workflows/validation.yml). It installs the binary reference wheel only for testing and uploads JSON, Markdown, and JUnit results where available. The workflow runs on pushes to `main` and can also be started manually. Missing contracts, blocked reference imports, and unapproved differences fail validation.
