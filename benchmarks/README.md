<!-- SPDX-License-Identifier: GPL-2.0-or-later -->

# Benchmarks

These charts compare scalar and wind analysis and synthesis latencies for
PySpharm float32, DUCC float32, and `pyspharm-cpp` float32 and float64.
Each operation has its own zero-based linear time scale; compare bar lengths
within an operation, and use the printed values across operations.
Smaller is better. Speedup is PySpharm time divided by `pyspharm-cpp` float32 time.

All five cases were measured in one [GitHub Actions benchmark run](https://github.com/mwyau/pyspharm-cpp/actions/runs/40110235638)
on October 24, 2026, using an AMD EPYC 9V45 runner with four allocated CPUs,
AVX-512 (`x86-64-v4`) selected, CPython 3.12.15, NumPy 1.26.4,
PySpharm 1.0.9, and DUCC 0.42.0.
The 2.5° one-thread case uses one DUCC thread; the other four cases use
four DUCC threads. PySpharm remains serial. Each median uses
7 timed repetitions after 2 warmups without CPU pinning.
Float64 runs use casts of the same float32 input fields, and synthesis inputs
are computed before timing.

## Clenshaw–Curtis · 2.5° · 144 × 73 (one thread)

![Clenshaw–Curtis · 2.5° · 144 × 73 (one thread) SHT benchmark timings](timings-cc-2p5-1t.svg)

| Operation | PySpharm f32 | pyspharm-cpp f32 | DUCC f32 | pyspharm-cpp f64 | Speedup vs PySpharm |
| --- | ---: | ---: | ---: | ---: | ---: |
| Scalar analysis | 0.068 ms | 0.217 ms | 0.215 ms | 0.221 ms | 0.32× |
| Scalar synthesis | 0.029 ms | 0.061 ms | 0.054 ms | 0.061 ms | 0.48× |
| Wind analysis | 0.255 ms | 0.435 ms | 0.428 ms | 0.581 ms | 0.59× |
| Wind synthesis | 0.117 ms | 0.126 ms | 0.108 ms | 0.250 ms | 0.93× |

## Clenshaw–Curtis · 2.5° · 144 × 73 (four threads)

![Clenshaw–Curtis · 2.5° · 144 × 73 (four threads) SHT benchmark timings](timings-cc-2p5.svg)

| Operation | PySpharm f32 | pyspharm-cpp f32 | DUCC f32 | pyspharm-cpp f64 | Speedup vs PySpharm |
| --- | ---: | ---: | ---: | ---: | ---: |
| Scalar analysis | 0.069 ms | 0.221 ms | 0.243 ms | 0.224 ms | 0.31× |
| Scalar synthesis | 0.030 ms | 0.069 ms | 0.064 ms | 0.075 ms | 0.44× |
| Wind analysis | 0.178 ms | 0.399 ms | 0.386 ms | 0.415 ms | 0.45× |
| Wind synthesis | 0.083 ms | 0.114 ms | 0.101 ms | 0.114 ms | 0.72× |

## Clenshaw–Curtis · 0.25° · 1440 × 721 (four threads)

![Clenshaw–Curtis · 0.25° · 1440 × 721 (four threads) SHT benchmark timings](timings-cc-0p25.svg)

| Operation | PySpharm f32 | pyspharm-cpp f32 | DUCC f32 | pyspharm-cpp f64 | Speedup vs PySpharm |
| --- | ---: | ---: | ---: | ---: | ---: |
| Scalar analysis | 93.747 ms | 12.657 ms | 13.958 ms | 12.699 ms | 7.41× |
| Scalar synthesis | 78.324 ms | 3.978 ms | 3.749 ms | 4.015 ms | 19.69× |
| Wind analysis | 243.113 ms | 28.321 ms | 28.858 ms | 31.550 ms | 8.58× |
| Wind synthesis | 108.406 ms | 10.806 ms | 9.274 ms | 14.694 ms | 10.03× |

## Gaussian · 128 latitude rings · 256 × 128 (four threads)

![Gaussian · 128 latitude rings · 256 × 128 (four threads) SHT benchmark timings](timings-gl-128.svg)

| Operation | PySpharm f32 | pyspharm-cpp f32 | DUCC f32 | pyspharm-cpp f64 | Speedup vs PySpharm |
| --- | ---: | ---: | ---: | ---: | ---: |
| Scalar analysis | 1.575 ms | 0.131 ms | 0.133 ms | 0.139 ms | 12.06× |
| Scalar synthesis | 0.126 ms | 0.125 ms | 0.113 ms | 0.146 ms | 1.01× |
| Wind analysis | 4.581 ms | 0.273 ms | 0.241 ms | 0.313 ms | 16.77× |
| Wind synthesis | 0.300 ms | 0.249 ms | 0.222 ms | 0.285 ms | 1.21× |

## Gaussian · 640 latitude rings · 1280 × 640 (four threads)

![Gaussian · 640 latitude rings · 1280 × 640 (four threads) SHT benchmark timings](timings-gl-640.svg)

| Operation | PySpharm f32 | pyspharm-cpp f32 | DUCC f32 | pyspharm-cpp f64 | Speedup vs PySpharm |
| --- | ---: | ---: | ---: | ---: | ---: |
| Scalar analysis | 142.171 ms | 3.361 ms | 2.959 ms | 3.574 ms | 42.30× |
| Scalar synthesis | 45.133 ms | 2.908 ms | 2.901 ms | 3.004 ms | 15.52× |
| Wind analysis | 349.341 ms | 7.901 ms | 6.260 ms | 8.908 ms | 44.21× |
| Wind synthesis | 74.312 ms | 7.571 ms | 6.445 ms | 11.460 ms | 9.81× |

## Reproduce

The benchmarks compare the same truncation on each backend. The regular-grid
analysis limit is `nlat-2`, Gaussian is `nlat-1`, and even-longitude
vector Nyquist modes are excluded. The benchmark uses seeded float32 inputs.

Install the project and the optional `pyspharm-syl==1.2.1` and
`ducc0==0.42.0` references into a CPython 3.12 environment with NumPy 1.26.4.
From the repository root, run:

```sh
python benchmarks/benchmark.py --threads 4 --warmups 2 --repeats 7 \
  --output /path/to/measurements.csv
# To reproduce the five charts, also run:
python benchmarks/benchmark.py --cases cc-2p5-1t --threads 1 \
  --warmups 2 --repeats 7 --output /path/to/one-thread.csv
python benchmarks/plot.py \
  --csv /path/to/measurements.csv /path/to/one-thread.csv \
  --output-dir benchmarks
```

The reproduced four-case series explicitly uses four threads. `--threads 0`
selects up to eight available CPUs, respecting Linux affinity and cgroup
quota; use an explicit thread count for controlled runs.
`python benchmarks/plot.py` requires both CSV files from one runner and
fails if any timing series is missing. Raw CSVs are not tracked. Independent numerical validation is in [`tests/parity/`](../tests/parity).
