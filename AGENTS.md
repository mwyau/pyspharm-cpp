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

- Use point-lane SIMD for `specintrp_many`. Coefficient-lane SIMD was
  approximately 4.7× slower in float32 and 6.4× slower in float64.
  Point-lane SIMD improved selected float32 cases by about 1.6–1.7×;
  float64 results were mixed. Keep native-width SIMD plus a scalar tail
  for `legendre_many`: a four-lane tail slowed small batches.
- Parallelize `legendre_many` and `specintrp_many` only when there are at
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
  public API and documented differences. PySpharm validation runs manually.
- Update `src/spharm_cpp/_spharm_cpp.pyi` when nanobind signatures or
  docstrings change. CI regenerates the stub, formats it, and compares it
  with the checked-in file and `py.typed`.
- Run `uv run --locked prek run --all-files` and
  `uv run --locked python -m pytest -q` before submitting changes.

## Packaging and releases

- Use uv and the committed `uv.lock` for development dependencies,
  tests, editable installations, and local builds.
- Keep the source distribution self-contained, with ducc0 source and
  project and third-party license notices. Preserve wheel and source
  distribution attribution.
- Do not publish releases, create tags, change repository visibility, or
  force-push without the owner's approval.
