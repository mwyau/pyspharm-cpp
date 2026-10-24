<!-- SPDX-License-Identifier: GPL-2.0-or-later -->

# Comparison with PySpharm

`spharm_cpp` provides a PySpharm-style Python API using DUCC transforms. The implementations agree on the tested common input domain, subject to the differences below. Agreement with another implementation is evidence about compatibility, not proof of accuracy for every grid and degree.

The reference tests use `pyspharm-syl==1.2.1` (module version 1.0.9 and a compiled SPHEREPACK backend), NumPy 1.26.4, and CPython 3.12 on Linux x86-64. The original Python interface was compared with [`jswhit/pyspharm` at `a464ebf`](https://github.com/jswhit/pyspharm/tree/a464ebfed334188fc0ef275fad4aafd2d9fa716e). The `pyspharm-syl` Python source removes the upstream `Spharmt.__repr__` and changes several equivalent rank-adaptation operations.

## Tested operations

The independent tests exercise the following interfaces, including scalar and batched shapes where applicable.

| Operation                                            | Comparison                                                                                   |
| ---------------------------------------------------- | -------------------------------------------------------------------------------------------- |
| `Spharmt`, attributes, `dtype`, representation       | Constructor, defaults, configuration, and output precision                                   |
| `grdtospec`, `spectogrd`                             | Scalar synthesis/analysis on regular and Gaussian grids; coefficient order; boundary degrees |
| `getvrtdivspec`, `getuv`, `getgrad`, `getpsichi`     | Eastward/northward signs, derivative scaling, vorticity, divergence, and potentials          |
| `specsmooth`, `regrid`                               | Smoothing and regridding at an explicit shared truncation                                    |
| `gaussian_lats_wts`, `getspecindx`, `getgeodesicpts` | Quadrature, index ordering, and geodesic point coordinates                                   |
| `legendre`, `specintrp`                              | Normalization, phase, polar behavior, and longitude                                          |
| `legendre_many`, `specintrp_many`                    | Results against repeated scalar operations                                                   |
| `cpu_info`                                           | DUCC CPU metadata keys and types                                                             |

The full public contract inventory is enforced by `tests/parity/manifest.py`. It contains 20 public API items and 11 accepted differences. Each item requires a marked test; missing tests fail validation.

## Supported differences

| Behavior                          | `spharm_cpp`                                                                                 |
| --------------------------------- | -------------------------------------------------------------------------------------------- |
| Module and backend                | `spharm_cpp` provides a DUCC-based implementation.                                           |
| Legendre evaluation               | `legfunc="computed"` is the default and only supported mode; `"stored"` raises `ValueError`. |
| Precision                         | float64/complex128 by default, float32/complex64 when requested.                             |
| Regular-grid analysis             | Maximum degree `nlat-2`, also bounded by longitude sampling; PySpharm accepts `nlat-1`.      |
| Scalar outputs at truncation zero | Coefficient and index arrays keep an explicit length-one axis instead of becoming scalars.   |
| Geodesic coordinates              | Calculated in float64; `dtype=np.float32` casts the output coordinates.                      |
| Regridding                        | Supports conversion between the two native precision modes.                                  |
| Representation                    | `Spharmt.__repr__` includes grid and precision information.                                  |
| Input validation                  | Invalid shapes, dtypes, parameter values, and unsupported truncations raise explicit errors. |
| Regular maximum-degree synthesis  | Evaluates odd-order harmonics at degree `nlat-1` that the tested SPHEREPACK synthesis omits. |
| `getpsichi` at `ntrunc=0`         | Returns two zero maps; the reference raises `IndexError`.                                    |

These exceptions are tested rather than skipped. See [numerical conventions](conventions.md) for the supported grid and vector-analysis limits.

## Highest regular-grid degree

For regular grids with `nlat` latitude samples, spherical harmonics with `m>0` vanish at the poles. At degree `nlat-1`, the `m=1` expansion has more coefficients than independent interior latitude samples. Coefficient recovery is therefore not unique.

PySpharm/SPHEREPACK accepts the full coefficient array for synthesis, but the tested `SHSES`/`SHSEC` routines omit odd-`m` contributions at degree `nlat-1`. DUCC evaluates these harmonics. Independent Legendre tests confirm that the corresponding basis functions need not vanish; this is an observed implementation difference, separate from the `m=1` sampling nullspace.

`spharm_cpp` limits regular-grid analysis to `nlat-2`, while synthesis accepts degree `nlat-1`. Gaussian grids support analysis through `nlat-1` when adequately resolved in longitude.

## Longitude Nyquist

For even `nlon`, the scalar order `m=nlon/2` samples `e^{im\lambda_j}=(-1)^j`. Its imaginary coefficient component cannot be determined from the samples. The scalar analysis conversion applies a half-factor at this order, matching PySpharm's positive/negative-frequency convention in direct comparisons. Analytic basis tests separately check the sampled field.

Vector analysis excludes the even-longitude Nyquist order; `getvrtdivspec` and `getpsichi` use the upper longitude limit `(nlon-1)//2`. Tests of vector synthesis include analytic gradient and wind fields at the boundary, while direct SPHEREPACK vector comparisons use resolved orders below Nyquist.

## Precision and validation evidence

Float32 numerical comparisons use operation-specific absolute and relative tolerances, rather than bitwise identity. Selected float64 Legendre functions are checked against an independent high-precision recurrence; basis synthesis uses Decimal reference samples with a `5e-14` tolerance. Geodesic float32 longitude comparisons account for the 0°/360° seam and test the rounded float64 result.

Run the PySpharm comparison with the manual [`validation.yml` workflow](https://github.com/mwyau/pyspharm-cpp/actions/workflows/validation.yml). It installs the binary reference wheel only for testing and uploads JSON, Markdown, and JUnit results where available. Missing contracts, blocked reference imports, and unapproved differences fail validation.

An [earlier hosted comparison](https://github.com/mwyau/pyspharm-cpp/actions/runs/40060787028) passed 37 cases at commit `cdbee944` on October 24, 2026. That result does not establish parity for later commits until they have been tested separately.
