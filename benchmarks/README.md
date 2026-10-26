# Benchmarks

These charts compare scalar and wind spherical harmonic transform (SHT)
analysis and synthesis latencies for PySpharm `float32`, `ducc0` `float32`, and
`pyspharm-cpp` `float32` and `float64`. Each operation has its own zero-based linear
time scale; compare bar lengths within an operation and use the printed values
across operations. Smaller is better. Speedup is the PySpharm median divided
by the `pyspharm-cpp` `float32` median.

The measurements were collected on 2026-10-25 on an AMD Ryzen 9 5900XT
16-Core Processor. Linux reported 16 logical CPUs; process affinity exposed
CPUs 0–15. The benchmark found 16 effective CPUs. No cgroup CPU quota
interface was available in this WSL2 environment. The OS was Linux
6.6.87.2-microsoft-standard-WSL2 x86-64 with glibc 2.43. The processor exposes
AVX2 but not AVX-512. The `pyspharm-cpp` extension was built in Release mode
with GNU C++ 15.2.0 and CMake 4.4.4 using the `multiarch` profile. The active
SIMD profile was `x86-64-v3`; x86-64-v4 was unavailable and was not benchmarked.

The software versions were CPython 3.14.9, NumPy 2.5.4, PySpharm 1.0.9,
`ducc0` 0.42.0, and `pyspharm-cpp` 0.1.0.dev0. PySpharm was installed
from conda-forge, and `pyspharm-cpp` from the current checkout.

The 2.5° one-thread case uses one DUCC thread; the other four cases use eight
DUCC threads. PySpharm is serial. The main speedup comparison therefore
compares serial PySpharm with eight-thread `ducc0` and `pyspharm-cpp` for those
four cases. The `ducc0` spin-1 transforms form the wind baseline;
`pyspharm-cpp` wind operations also convert between wind fields and spectral
coefficients for relative vorticity and divergence. Each median uses seven
timed repetitions after two warmups without CPU pinning. The p90 is the
nearest-rank 90th percentile, which is the maximum of seven samples. Float64
runs use casts of the same `float32` input fields, and synthesis inputs are
computed before timing. Inputs use NumPy seed 20261024.

The same spectral truncation is used for all backends and operations:
min(nlat - 2, (nlon - 1) // 2) for Clenshaw–Curtis grids and
min(nlat - 1, (nlon - 1) // 2) for Gaussian grids.

## Clenshaw–Curtis · 2.5° · 144 × 73 (one thread)

![Clenshaw–Curtis · 2.5° · 144 × 73 (one thread) SHT benchmark timings](timings-cc-2p5-1t.svg)

| Operation | PySpharm f32 | `pyspharm-cpp` f32 | `ducc0` f32 | `pyspharm-cpp` f64 | Speedup vs PySpharm |
| --- | ---: | ---: | ---: | ---: | ---: |
| Scalar analysis | 0.090 ms | 0.254 ms | 0.251 ms | 0.250 ms | 0.36× |
| Scalar synthesis | 0.047 ms | 0.084 ms | 0.096 ms | 0.086 ms | 0.56× |
| Wind analysis | 0.202 ms | 0.536 ms | 0.512 ms | 0.509 ms | 0.38× |
| Wind synthesis | 0.112 ms | 0.190 ms | 0.192 ms | 0.185 ms | 0.59× |

## Clenshaw–Curtis · 2.5° · 144 × 73 (eight threads)

![Clenshaw–Curtis · 2.5° · 144 × 73 (eight threads) SHT benchmark timings](timings-cc-2p5.svg)

| Operation | PySpharm f32 | `pyspharm-cpp` f32 | `ducc0` f32 | `pyspharm-cpp` f64 | Speedup vs PySpharm |
| --- | ---: | ---: | ---: | ---: | ---: |
| Scalar analysis | 0.121 ms | 0.475 ms | 0.470 ms | 0.510 ms | 0.25× |
| Scalar synthesis | 0.044 ms | 0.295 ms | 0.328 ms | 0.345 ms | 0.15× |
| Wind analysis | 0.206 ms | 0.742 ms | 0.652 ms | 0.874 ms | 0.28× |
| Wind synthesis | 0.111 ms | 0.379 ms | 0.392 ms | 0.592 ms | 0.29× |

## Clenshaw–Curtis · 0.25° · 1440 × 721 (eight threads)

![Clenshaw–Curtis · 0.25° · 1440 × 721 (eight threads) SHT benchmark timings](timings-cc-0p25.svg)

| Operation | PySpharm f32 | `pyspharm-cpp` f32 | `ducc0` f32 | `pyspharm-cpp` f64 | Speedup vs PySpharm |
| --- | ---: | ---: | ---: | ---: | ---: |
| Scalar analysis | 111.258 ms | 7.496 ms | 6.817 ms | 7.288 ms | 14.84× |
| Scalar synthesis | 57.051 ms | 3.412 ms | 2.832 ms | 3.325 ms | 16.72× |
| Wind analysis | 258.342 ms | 15.657 ms | 13.686 ms | 38.720 ms | 16.50× |
| Wind synthesis | 111.987 ms | 7.746 ms | 6.456 ms | 24.865 ms | 14.46× |

## Gaussian · 128 latitude rings · 256 × 128 (eight threads)

![Gaussian · 128 latitude rings · 256 × 128 (eight threads) SHT benchmark timings](timings-gl-128.svg)

| Operation | PySpharm f32 | `pyspharm-cpp` f32 | `ducc0` f32 | `pyspharm-cpp` f64 | Speedup vs PySpharm |
| --- | ---: | ---: | ---: | ---: | ---: |
| Scalar analysis | 1.843 ms | 0.314 ms | 0.398 ms | 0.307 ms | 5.88× |
| Scalar synthesis | 0.157 ms | 0.336 ms | 0.334 ms | 0.355 ms | 0.47× |
| Wind analysis | 3.201 ms | 0.425 ms | 0.356 ms | 0.502 ms | 7.53× |
| Wind synthesis | 0.447 ms | 0.409 ms | 0.344 ms | 0.427 ms | 1.09× |

## Gaussian · 640 latitude rings · 1280 × 640 (eight threads)

![Gaussian · 640 latitude rings · 1280 × 640 (eight threads) SHT benchmark timings](timings-gl-640.svg)

| Operation | PySpharm f32 | `pyspharm-cpp` f32 | `ducc0` f32 | `pyspharm-cpp` f64 | Speedup vs PySpharm |
| --- | ---: | ---: | ---: | ---: | ---: |
| Scalar analysis | 139.024 ms | 2.103 ms | 2.411 ms | 2.206 ms | 66.11× |
| Scalar synthesis | 32.081 ms | 2.227 ms | 2.304 ms | 2.682 ms | 14.41× |
| Wind analysis | 339.515 ms | 4.890 ms | 3.389 ms | 28.454 ms | 69.42× |
| Wind synthesis | 67.204 ms | 5.432 ms | 4.558 ms | 7.392 ms | 12.37× |

| Grid case | Largest speedup | Operation |
| --- | ---: | --- |
| Clenshaw–Curtis 144 × 73, one thread | 0.59× | Wind synthesis |
| Clenshaw–Curtis 144 × 73, eight threads | 0.29× | Wind synthesis |
| Clenshaw–Curtis 1440 × 721, eight threads | 16.72× | Scalar synthesis |
| Gaussian 256 × 128, eight threads | 7.53× | Wind analysis |
| Gaussian 1280 × 640, eight threads | 69.42× | Wind analysis |

Across the five cases, `pyspharm-cpp` `float32` was slower than serial PySpharm
for all four operations on the 144 × 73 grid at both thread counts. On the
larger 1440 × 721 and 1280 × 640 grids, `pyspharm-cpp` was faster for all four
operations. The 256 × 128 Gaussian case was mixed: scalar and wind analysis
were faster, wind synthesis was close to parity, and scalar synthesis was
slower. The largest measured ratio was 69.42× for wind analysis on the
1280 × 640 Gaussian grid. These ratios include the difference between serial
PySpharm and eight DUCC threads for the multithreaded cases.

## Legendre and interpolation overloads

The auxiliary benchmark compares repeated scalar calls with one batched call
for four cases: 8 points at `lmax` 16, 64 points at `lmax` 64, 128 points at
`lmax` 256, and 16 points at `lmax` 1024. Measurements below used `float32`
and `float64`, two warmups, seven repetitions, and separate processes at one,
two, and four threads on an AMD Ryzen 9 5950X with the active `x86-64-v3`
profile. Ratios are scalar-call median divided by batched-call median,
summarized as ranges across the four cases. Ratios above one favor batching.

| Threads | Operation | Float32 ratio range | Float64 ratio range |
| ---: | --- | ---: | ---: |
| 1 | Legendre | 2.26–4.92× | 1.09–4.37× |
| 1 | Spectral interpolation | 1.88–2.80× | 1.68–2.29× |
| 2 | Legendre | 2.24–6.24× | 1.74–4.64× |
| 2 | Spectral interpolation | 1.91–3.63× | 1.56–3.28× |
| 4 | Legendre | 2.23–7.84× | 2.96–6.29× |
| 4 | Spectral interpolation | 1.90–3.90× | 2.19–5.66× |

The current scalar and batched implementations are bitwise identical at each
tested dtype and thread count. The ratios above therefore compare equivalent
outputs. Independent accuracy checks remain separate from these timing tests.

The previously recorded full-harness `small-high` measurements (16 points,
`lmax` 1024) were:

| dtype | Operation | Before | After |
| --- | --- | ---: | ---: |
| float32 | Legendre scalar | 50.471 | 82.739 |
| float32 | Legendre batched | 35.936 | 36.659 |
| float32 | Interpolation scalar | 7.392 | 5.668 |
| float32 | Interpolation batched | 3.343 | 2.995 |
| float64 | Legendre scalar | 67.118 | 67.379 |
| float64 | Legendre batched | 53.540 | 52.269 |
| float64 | Interpolation scalar | 7.742 | 5.908 |
| float64 | Interpolation batched | 3.520 | 3.517 |

The float32 scalar Legendre kernel now performs its recurrence in binary32.
Before the precision correction, `legendre_values<float>` called the templated
kernel with `std::sin(phi)` and `std::cos(phi)`, both `double`, so template
deduction selected a binary64 recurrence and storage narrowed its results. The
current code explicitly converts the shared binary64 angle setup to `T` and
calls `legendre_scalar_kernel<T>`. The historical before/after comparison
therefore measures different arithmetic and does not establish equivalent
float32 outputs. The independent accuracy tests validate the native-precision
result separately.

### Native float32 Legendre scaling cost

Profiling of one scalar `legendre(78°, 1024, dtype=float32)` call found
481,663 `scalbnf` calls in the unoptimized native-float implementation. Of
those, 286,367 had a scale exponent at or below -270 and returned signed zero.
The recurrence normalization keeps each restored value at or below 2^120, so
these results are at most 2^-150 in magnitude and round to signed zero under
the default round-to-nearest mode. `restore_legendre_scale` now returns that
signed zero directly. This changes neither the recurrence nor its precision;
values above the cutoff still use `scalbnf`.

The focused stage benchmark used the existing `small-high` workload and timed
only its scalar Legendre operation. Each process used two warmups and nine
samples; three independent processes per implementation ran pinned to CPUs
0–3 in alternating order. Values report the median of the three process
medians; brackets give their range:

| Implementation | float32 median | vs. legacy | float64 median |
| --- | ---: | ---: | ---: |
| Legacy scalar recurrence (`dc11c38`) | 65.445 ms [64.880, 68.591] | — | 73.001 ms [71.622, 79.868] |
| Native float32 before cutoff (`fd32e26`) | 103.455 ms [97.929, 104.062] | +58.1% | 77.134 ms [72.291, 79.381] |
| Native float32 with cutoff | 72.279 ms [72.153, 77.088] | +10.4% | 75.886 ms [73.135, 78.200] |

The cutoff reduced the native-float scalar time by 30.1% versus its parent.
It recovered most, but not all, of the difference from the former binary64
recurrence. The float64 differences are within the observed process spread;
this change does not alter the float64 kernel. The per-process measurements
and metadata are in
[`auxiliary-legendre-stage-comparison.csv`](auxiliary-legendre-stage-comparison.csv).

The high-degree result comes from the same native-float source except for the
zero shortcut. A separate Callgrind run, started after Python imports, counted
61.6 million instructions for the unoptimized native-float scalar call and
43.9 million with the shortcut; the old binary64 implementation counted
44.0 million. The unoptimized profile attributed 26.1 million instructions
(42.4%) to `scalbnf`; the optimized profile attributed 9.8 million (22.4%).
This supports repeated float32 scale restoration as the main avoidable cost.
The complete benchmark harness showed large shifts in operations unchanged by
this source edit, so its before/after timings are too noisy to support a
whole-harness speedup claim. The
[`auxiliary-precision-comparison.csv`](auxiliary-precision-comparison.csv)
retains those measurements for transparency: three processes at one, two, and
four threads across the four existing workloads, both dtypes, and all four
operations. It includes source revision, compiler, active dispatch, CPU
affinity, cgroup quota, warmups, repetitions, median, and p90. The focused
scalar stage benchmark is the performance evidence for this optimization.

The measurements used GCC 15.2.0, the `multiarch` build with active
`x86-64-v3`, an AMD Ryzen 9 5950X, and four pinned CPUs (0–3), with no cgroup
CPU quota. They used CPython 3.11.15, NumPy 1.26.4, and one thread for the
focused stage benchmark. Its source identifiers are `dc11c38` (the previous
binary64 scalar recurrence), `fd32e26` (native float32 before the cutoff), and
`fd32e26` plus `src/cpp/numerics.cpp` blob
`cbddac632cc12df7a273a38cb089a347b188546e` (native float32 with the cutoff).
The full-harness comparison CSV preserves its earlier candidate snapshot,
`ece8cd6505799a07f82ebc61be2b2c82a8cf521f`. x86-64-v4 was compiled but
unavailable on this host. Benchmark runs record source revision, compiler,
CPU, affinity, dispatch profile, thread count, warmups, and repetitions in
`BENCHMARK_METADATA`.

The focused exact-bit and accuracy tests passed for portable, x86-64, and
x86-64-v3 builds with one, two, and four DUCC threads. The full test suite
passed 161 tests, and the PySpharm parity suite passed 39 tests with 18 passes,
11 accepted differences, and no failures or blocked cases. Windows, macOS, and
x86-64-v4 execution were unavailable locally.

## Run benchmarks

The committed [measurements.csv](measurements.csv) contains 80 measurements:
five cases × four operations × four timing series (PySpharm `float32`, `ducc0`
`float32`, `pyspharm-cpp` `float32`, and `pyspharm-cpp` `float64`). It records the
thread count, truncation, median, and p90 for every measurement. The benchmark
rejects a requested thread count above the CPUs available from process
affinity and the cgroup CPU quota when that quota is exposed.

Create the conda benchmark environment and install the checkout with the
multiarch build profile:

    conda env create -f .github/conda/pyspharm-reference.yml
    conda run --no-capture-output -n pyspharm-reference python -m pip install \
      --no-build-isolation --no-deps -e . \
      --config-settings=cmake.define.SPHARM_CPP_BUILD_PROFILE=multiarch

Run the four eight-thread cases and the one-thread case in separate Python
processes because DUCC's thread configuration is initialized on its first
numerical call:

    conda run --no-capture-output -n pyspharm-reference python benchmarks/benchmark.py \
      --threads 8 --warmups 2 --repeats 7 --expected-profile x86-64-v3 \
      --output /tmp/sht-eight-thread.csv
    conda run --no-capture-output -n pyspharm-reference python benchmarks/benchmark.py \
      --cases cc-2p5-1t --threads 1 --warmups 2 --repeats 7 \
      --expected-profile x86-64-v3 --output /tmp/sht-one-thread.csv

Combine the measurements and regenerate the five figures:

    {
      head -n 1 /tmp/sht-eight-thread.csv
      tail -n +2 /tmp/sht-eight-thread.csv
      tail -n +2 /tmp/sht-one-thread.csv
    } > benchmarks/measurements.csv
    conda run --no-capture-output -n pyspharm-reference python benchmarks/plot.py \
      --csv benchmarks/measurements.csv --output-dir benchmarks

Add `--print-samples` to a benchmark command to print individual timing samples
when checking measurement variability. For CPU dispatch and affinity
information, the benchmark prints a `BENCHMARK_METADATA` record at startup.

Measure the auxiliary kernels in one, two, and four-thread processes:

    for threads in 1 2 4; do
      conda run --no-capture-output -n pyspharm-reference python \
        benchmarks/benchmark_auxiliary.py --threads "$threads" \
        --warmups 2 --repeats 7 --output "/tmp/auxiliary-$threads.csv"
    done

Select the existing high-truncation scalar Legendre workload to measure that
stage without timing the other operations:

    conda run --no-capture-output -n pyspharm-reference python \
      benchmarks/benchmark_auxiliary.py --threads 1 --warmups 2 --repeats 9 \
      --case small-high --operation legendre_scalar \
      --output /tmp/auxiliary-legendre-stage.csv

The optional `--case` and `--operation` filters select existing workloads;
without them the benchmark still runs every case and operation.

The figures can also be regenerated directly from the committed CSV:

    conda run --no-capture-output -n pyspharm-reference python benchmarks/plot.py \
      --csv benchmarks/measurements.csv --output-dir benchmarks

The plotter derives thread-count labels from the CSV and checks timing-series
completeness and consistent multithreaded thread counts. Independent numerical
validation is in [tests/parity/](../tests/parity).
