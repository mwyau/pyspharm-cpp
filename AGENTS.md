# Repository constraints

These rules apply to implementation, tests, and documentation.

## Numerical methods

- Scalar analysis must obey the longitude sampling limit `lmax <= nlon//2`.
  Clenshaw–Curtis analysis uses ducc0's native transform through
  `lmax <= nlat-2`; Gaussian analysis permits `lmax <= nlat-1`.
  Reject explicit truncations above these limits.
- Vector analysis also excludes the even-longitude Nyquist order:
  `lmax <= (nlon-1)//2`. Scalar and vector synthesis have their own
  truncation limits.
- Do not introduce SPHEREPACK `zbar` projection, a maximum-degree correction,
  a pseudoinverse, or another fallback to recover regular-grid
  `lmax=nlat-1` analysis. Change this boundary only with the owner's approval.
- Preserve coefficient normalization, packing, sign conventions, input
  precision, and the distinction between scalar and vector transforms.
  Test changed numerical behavior against analytic results and, where useful,
  the identified PySpharm implementation.

### Exact scalar/batch agreement for auxiliary numerical APIs

This is a required numerical contract for `legendre` and `specintrp`,
not an approximate-accuracy or performance target. It takes precedence over
the batch optimizations documented below.

- For the same finite coordinate values, coefficients, truncation, and chosen
  computation dtype, every row/element of a batched call must be **bitwise
  identical in that dtype** to the corresponding scalar call. This holds for
  one-element batches, full SIMD batches, partial SIMD tails, varied batch
  lengths and row positions, and supported input strides. Test under every
  applicable portable/multiarch dispatch profile and supported threading
  configuration. Different platforms or dispatch profiles need not produce
  the same bits as each other, but scalar and batch must agree within each.
- Test equality of the actual `float32`/`float64` output bit patterns,
  including signed zeros; a tolerance-based comparison (`allclose`) is
  insufficient. Scalar `specintrp` returns a Python `float`: convert it
  back to the input coefficients' real dtype before comparing its bits with
  the batched result, because Python represents that result as binary64.
  Keep separate, tolerance-based tests against independent scientific
  references and PySpharm; bitwise scalar/batch equality does not establish
  numerical accuracy.
- Preserve *native-precision recurrence and accumulation*: `float32` uses
  binary32 and `float64` uses binary64 for these operations in **both**
  overloads. Do not evaluate an entire `float32` result in binary64 merely
  to narrow it on return. Shared binary64 angle reduction, trigonometry,
  recurrence-constant preparation, or other explicitly justified setup is
  permitted only when both overloads use the same conversions and numerical
  sequence; it must not silently widen their recurrence/accumulation.
- `legendre` scalar and batch must agree on coordinate preprocessing,
  exact-pole handling, recurrence constants, multiply/add or FMA semantics,
  scaling/exponent updates, and final rounding. `specintrp` scalar and batch
  must agree on longitude reduction, phase evaluation and the 32-order
  recurrence-reset rule, coefficient grouping, summation order, multiply/add
  or FMA semantics, and final rounding. SIMD lanes, scalar tails, and thread
  scheduling must not change the answer for a given point. Prefer sharing
  per-point numerical logic over keeping independently evolving algorithms;
  never assume mathematically equivalent SIMD and scalar expressions are
  bitwise equivalent without testing.
- `getgeodesicpts` has a scalar size argument and returns arrays; it has no
  scalar-versus-batch overload. Its `dtype` still selects **native**
  `float32` or `float64` geometric computation and output, rather than
  always computing binary64 coordinates then casting. Match PySpharm's float32
  default; keep explicit `dtype=np.float64` as an extension. Validate geodesic
  geometry and ordering against independent references at appropriate
  dtype-specific tolerances; do not require its float32 values to equal
  float64 values narrowed to float32. Keep its public signature and point
  ordering unchanged.
- Do not extend bitwise scalar/batch requirements to unrelated SHT analysis
  or synthesis APIs without explicit approval. Preserve the existing public
  signatures, default dtypes, validation, return shapes, and coefficient
  conventions for these three auxiliary functions.

## ducc0 and native builds

- Use the pinned, unmodified Git submodule at `third_party/ducc`.
  Keep the C++ driver and Python bindings outside the submodule. Reuse
  ducc0's FFT, SHT, SIMD, threading, and CPU-dispatch implementations.
- Do not vendor a second copy of ducc0 or add a separate transform,
  FFT, or CPU-dispatch implementation.
- Linux x86-64 multiarch builds use v1/v3/v4 code paths with a baseline
  dispatcher. Other platforms default to portable builds. Code that executes
  before CPU dispatch, including thread-pool initialization, must use
  baseline instructions.
- Keep scikit-build-core, CMake, and nanobind as the build stack. NumPy is
  the only required Python runtime dependency. Isolated pip builds must not
  need uv, the ducc0 Python package, or network retrieval of the submodule.

## Performance decisions

These choices were measured on a Ryzen 9 5950X with the Linux x86-64-v3
profile (2 warmups, 7 timing repetitions). They are not general speedup
guarantees, particularly for ARM64, where batch scheduling is unmeasured.

- Prefer point-lane SIMD for array inputs to `specintrp`. Coefficient-lane
  SIMD was approximately 4.7× slower in float32 and 6.4× slower in float64.
  Point-lane SIMD improved selected float32 cases by about 1.6–1.7×;
  float64 results were mixed. Previously, native-width SIMD plus a scalar
  tail outperformed a four-lane tail on small `legendre` batches. These
  measurements do not override exact scalar/batch agreement: retain any
  SIMD or scalar-tail optimization only if it satisfies the numerical
  contract above.
- Parallelize array inputs to `legendre` and `specintrp` only when there are at
  least four SIMD batches. Small inputs did not benefit from thread
  scheduling; larger cases improved by up to 3.39× with four threads.
- For scalar analysis/synthesis and wind analysis, parallelize across
  fields when there are at least four fields and 32,768 pixels; otherwise
  use DUCC's internal threading. The four-thread batched speedup medians
  were 1.50×, 1.34×, and 1.45×, respectively, on medium/large Gaussian grids.
- Keep wind synthesis sequential across fields, using DUCC's internal
  threading. The outer field-workspace variant was slower: 16-field
  float64 synthesis increased from 5.86 to 9.30 ms (256 × 128 Gaussian)
  and from 38.5 to 53.3 ms (512 × 256 Gaussian) with one thread.
- Write vector synthesis and gradient results directly into owner-backed
  NumPy component arrays; do not restore intermediate-map copying or
  separate small-batch fallback paths.
- Preserve the interpolation phase-recurrence reset every 32 orders to
  control accumulated error; do not treat this as a purely performance
  parameter. Keep all scheduling and phase-reset thresholds as local
  `constexpr` values, not centralized or publicly configurable.
- Retain `native` as an optional developer build profile without native CI
  tests. Linux x86-64 distribution wheels, including free-threaded wheels,
  use `multiarch`; other platforms use `portable`.

## Documentation and checks

- Follow [write-like-a-scientist](https://github.com/mwyau/write-like-a-scientist/blob/main/skills/write-like-a-scientist/SKILL.md)
  for scientific prose and docstrings. Use established terminology, including
  *atmospheric science*, *relative vorticity*, and *spectral truncation*.
  State scientific limits, units, and evidence precisely.
- Describe current behavior in README and `docs/`; avoid dated CI results,
  development history, and work logs. README and CONTRIBUTING cover usage
  and contribution instructions; API details belong in `docs/`.
- Record concise implementation rationale and profiling-backed decisions in
  this file. Do not add standalone plans, progress logs, or profiling notes
  in Markdown; keep published measurements and reproduction instructions
  in `benchmarks/README.md`.
- Keep the parity manifest and its behavior tests synchronized with the
  public API and documented differences.
- Update `src/spharm_cpp/_spharm_cpp.pyi` when nanobind signatures or
  docstrings change. CI regenerates the stub, formats it, and compares it
  with the checked-in file and `py.typed`.
- For changes to `legendre`, `specintrp`, or `getgeodesicpts`,
  add regression tests for the exact precision/bit-pattern contract above.
  Include dtype and SIMD-boundary coverage, and update any tests, parity
  manifest entries, or documentation that currently assert float32 results
  are rounded copies of float64. Record performance regressions or wins
  using the existing auxiliary benchmark; never loosen exact agreement
  to recover a speedup.
- Run `uv run --locked prek run --all-files` and
  `uv run --locked python -m pytest -q` before submitting changes.

## CI and releases

- `.github/workflows/ci.yml` runs on pull requests, pushes to `main`, and
  manual dispatch. It tests supported Python versions, free-threaded and
  ABI3 wheels, source builds, code quality, documentation, and C++ analysis.
  Build and test Linux, macOS, and Windows distribution wheels only on
  pushes to `main` or manual dispatch, not on pull requests.
- `.github/workflows/coverage.yml` instruments a portable C++ build.
  `.github/workflows/validation.yml` checks PySpharm parity on pushes
  to `main` and manual dispatch.
- `.github/workflows/python-publish.yml` builds release wheels and the
  source distribution from a version tag after checking the package
  version. `vX.Y.Z.devN` tags publish to TestPyPI; `vX.Y.Z` tags publish
  to PyPI, followed by a GitHub Release containing only the source
  distribution (no wheel assets). Publishing uses trusted publishers
  and requires approval from the corresponding `testpypi` or `pypi`
  environment. Configure required reviewers for both environments and
  register each as a Trusted Publisher; do not store publishing credentials.
  Protect `v*` tags against modification or deletion and enable release
  immutability in repository settings.
- Use uv and the committed `uv.lock` for development dependencies,
  tests, editable installations, and local builds.
- Keep the source distribution self-contained, with ducc0 source and
  project and third-party license notices. Preserve wheel and source
  distribution attribution.
- Do not publish releases, create tags, change repository visibility, or
  force-push without the owner's approval.
