# Benchmarks

These charts compare scalar and wind spherical harmonic transform (SHT)
analysis and synthesis latencies for `pyspharm` float32, ducc0 float32, and
`pyspharm-cpp` float32 and float64.
Each operation has its own zero-based linear time scale; compare bar lengths
within an operation and use the printed values across operations. Smaller is
better. Speedup is `pyspharm` time divided by `pyspharm-cpp` float32 time.

The measurements ran on an AMD Ryzen 9 5950X 16-Core Processor with 32 logical
CPUs. Process affinity exposed CPUs 0–31, and the cgroup CPU quota was unlimited,
so 32 CPUs were available. The OS was Linux 7.0.14-20-pve x86-64 with glibc
2.43. The
`pyspharm-cpp` active SIMD profile was `x86-64-v3` (AVX2); the build includes
the x86-64, x86-64-v3, and x86-64-v4 profiles. The ducc0 wheel's active CPU
dispatch profile was not inspected.

The software versions were CPython 3.11.15, NumPy 1.26.4, `pyspharm` 1.0.9,
ducc0 0.42.0, `pyspharm-cpp` 0.1.0.dev0, uv 0.12.23, GNU Fortran 15.2.0,
setuptools 59.8.0, and wheel 0.48.0. The `pyspharm` source distribution uses
`numpy.distutils`; it was built on CPython 3.11 with
`FFLAGS=-fallow-argument-mismatch`. These versions record the environment that
produced the committed measurements.

The 2.5° one-thread case uses one DUCC thread; the other four cases use eight
DUCC threads. `pyspharm` is serial. The DUCC spin-1 transforms form the wind
baseline; `getvrtdivspec` and `getuv` also convert between wind fields and
spectral coefficients for relative vorticity and divergence. Each median uses seven timed
repetitions after two warmups without CPU pinning. The p90 is the nearest-rank
90th percentile, which is the maximum of seven samples. Float64 runs use casts
of the same float32 input fields, and synthesis inputs are computed before
timing. Inputs use NumPy seed 20261024.

The same spectral truncation is used for all backends and operations:
`min(nlat - 2, (nlon - 1) // 2)` for Clenshaw–Curtis grids and
`min(nlat - 1, (nlon - 1) // 2)` for Gaussian grids.

## Clenshaw–Curtis · 2.5° · 144 × 73 (one thread)

![Clenshaw–Curtis · 2.5° · 144 × 73 (one thread) SHT benchmark timings](timings-cc-2p5-1t.svg)

| Operation | pyspharm f32 | pyspharm-cpp f32 | ducc0 f32 | pyspharm-cpp f64 | Speedup vs pyspharm |
| --- | ---: | ---: | ---: | ---: | ---: |
| Scalar analysis | 0.071 ms | 0.252 ms | 0.252 ms | 0.239 ms | 0.28× |
| Scalar synthesis | 0.041 ms | 0.084 ms | 0.086 ms | 0.098 ms | 0.48× |
| Wind analysis | 0.186 ms | 0.519 ms | 0.520 ms | 0.501 ms | 0.36× |
| Wind synthesis | 0.107 ms | 0.190 ms | 0.183 ms | 0.213 ms | 0.57× |

## Clenshaw–Curtis · 2.5° · 144 × 73 (eight threads)

![Clenshaw–Curtis · 2.5° · 144 × 73 (eight threads) SHT benchmark timings](timings-cc-2p5.svg)

| Operation | pyspharm f32 | pyspharm-cpp f32 | ducc0 f32 | pyspharm-cpp f64 | Speedup vs pyspharm |
| --- | ---: | ---: | ---: | ---: | ---: |
| Scalar analysis | 0.075 ms | 0.217 ms | 0.246 ms | 0.213 ms | 0.35× |
| Scalar synthesis | 0.041 ms | 0.065 ms | 0.078 ms | 0.064 ms | 0.63× |
| Wind analysis | 0.190 ms | 0.411 ms | 0.433 ms | 0.392 ms | 0.46× |
| Wind synthesis | 0.107 ms | 0.094 ms | 0.094 ms | 0.109 ms | 1.13× |

## Clenshaw–Curtis · 0.25° · 1440 × 721 (eight threads)

![Clenshaw–Curtis · 0.25° · 1440 × 721 (eight threads) SHT benchmark timings](timings-cc-0p25.svg)

| Operation | pyspharm f32 | pyspharm-cpp f32 | ducc0 f32 | pyspharm-cpp f64 | Speedup vs pyspharm |
| --- | ---: | ---: | ---: | ---: | ---: |
| Scalar analysis | 62.899 ms | 5.996 ms | 6.677 ms | 6.336 ms | 10.49× |
| Scalar synthesis | 57.463 ms | 2.774 ms | 2.772 ms | 3.228 ms | 20.71× |
| Wind analysis | 168.856 ms | 13.711 ms | 13.455 ms | 25.338 ms | 12.32× |
| Wind synthesis | 93.492 ms | 6.842 ms | 6.578 ms | 18.792 ms | 13.67× |

## Gaussian · 128 latitude rings · 256 × 128 (eight threads)

![Gaussian · 128 latitude rings · 256 × 128 (eight threads) SHT benchmark timings](timings-gl-128.svg)

| Operation | pyspharm f32 | pyspharm-cpp f32 | ducc0 f32 | pyspharm-cpp f64 | Speedup vs pyspharm |
| --- | ---: | ---: | ---: | ---: | ---: |
| Scalar analysis | 1.746 ms | 0.110 ms | 0.102 ms | 0.111 ms | 15.82× |
| Scalar synthesis | 0.143 ms | 0.105 ms | 0.099 ms | 0.114 ms | 1.37× |
| Wind analysis | 1.937 ms | 0.220 ms | 0.183 ms | 0.223 ms | 8.79× |
| Wind synthesis | 0.410 ms | 0.211 ms | 0.175 ms | 0.224 ms | 1.94× |

## Gaussian · 640 latitude rings · 1280 × 640 (eight threads)

![Gaussian · 640 latitude rings · 1280 × 640 (eight threads) SHT benchmark timings](timings-gl-640.svg)

| Operation | pyspharm f32 | pyspharm-cpp f32 | ducc0 f32 | pyspharm-cpp f64 | Speedup vs pyspharm |
| --- | ---: | ---: | ---: | ---: | ---: |
| Scalar analysis | 299.388 ms | 1.690 ms | 1.793 ms | 2.313 ms | 177.15× |
| Scalar synthesis | 31.642 ms | 1.977 ms | 1.916 ms | 2.250 ms | 16.01× |
| Wind analysis | 289.422 ms | 4.408 ms | 3.328 ms | 6.115 ms | 65.66× |
| Wind synthesis | 64.565 ms | 4.923 ms | 4.281 ms | 9.045 ms | 13.12× |

## Reproduce

The committed [measurements.csv](measurements.csv) contains 80 measurements:
five cases × four operations × four backends. It records the thread count,
truncation, median, and p90 for every measurement. The benchmark rejects a
requested thread count above the CPUs available from process affinity and the
cgroup CPU quota.

The manually dispatched GitHub Actions workflow requires an AVX-512-capable
`ubuntu-24.04` runner, but uses the maximum CPU count available through affinity
and cgroup quota instead of requiring eight CPUs. It runs the five cases
separately under each `pyspharm-cpp` dispatch profile (`x86-64-v1`,
`x86-64-v3`, and `x86-64-v4`) and uploads per-profile CSVs and figures.
Those runner measurements are separate from the committed eight-thread,
AVX2 measurements below.

The `pyspharm==1.0.9` distribution is a source release. Its legacy
Fortran build uses GNU Fortran, NumPy 1.26.4, patched setuptools 84.1.0, wheel
0.48.0, and `FFLAGS=-fallow-argument-mismatch` with current GNU Fortran. The
setuptools update changes the build tool while retaining the reference source
and compiler configuration; it does not replace or regenerate the committed
measurements above. Set up the locked benchmark environment from the repository
root:

```sh
uv python install 3.11
uv sync --locked --python 3.11 --group benchmark \
  --no-install-project --no-install-package pyspharm
FFLAGS="-fallow-argument-mismatch" uv pip install \
  --python .venv/bin/python --no-build-isolation --no-deps \
  --no-binary pyspharm "pyspharm==1.0.9"
FFLAGS="-fallow-argument-mismatch" uv sync --locked --python 3.11 \
  --group benchmark --no-build-isolation
```

Run the four eight-thread cases and the one-thread case in separate Python
processes because DUCC's thread configuration is initialized once:

```sh
uv run --locked --no-sync --python 3.11 python benchmarks/benchmark.py \
  --threads 8 --warmups 2 --repeats 7 \
  --output /tmp/sht-eight-thread.csv
uv run --locked --no-sync --python 3.11 python benchmarks/benchmark.py \
  --cases cc-2p5-1t --threads 1 --warmups 2 --repeats 7 \
  --output /tmp/sht-one-thread.csv
```

Combine the files to replace the committed measurements, then regenerate all
five figures:

```sh
{
  head -n 1 /tmp/sht-eight-thread.csv
  tail -n +2 /tmp/sht-eight-thread.csv
  tail -n +2 /tmp/sht-one-thread.csv
} > benchmarks/measurements.csv
uv run --locked --no-sync --python 3.11 python benchmarks/plot.py \
  --csv benchmarks/measurements.csv --output-dir benchmarks
```

Add `--print-samples` to a benchmark command to print the individual timing
samples when checking measurement variability. For CPU dispatch and affinity
information, the benchmark prints a `BENCHMARK_METADATA` record at startup.

Implementation choices for SIMD and batched threading are recorded in
[AGENTS.md](../AGENTS.md). To measure the auxiliary kernels at 1, 2, and 4 threads,
run one process per thread count because the setting is read on the first
numerical call:

```sh
for threads in 1 2 4; do
  uv run --locked python benchmarks/benchmark_auxiliary.py \
    --threads "$threads" --warmups 2 --repeats 7 \
    --output "/tmp/auxiliary-${threads}.csv"
done
```

The figures can also be regenerated directly from the committed CSV:

```sh
uv run --locked --no-sync --python 3.11 python benchmarks/plot.py \
  --csv benchmarks/measurements.csv --output-dir benchmarks
```

The plotter derives thread-count labels from the CSV and fails if any
timing series is missing or the multithreaded cases disagree on thread count.
Independent numerical validation is
in [`tests/parity/`](../tests/parity).
