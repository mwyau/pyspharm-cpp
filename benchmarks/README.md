<!-- SPDX-License-Identifier: GPL-2.0-or-later -->

# Benchmarks

These charts compare median timings for PySpharm float32, DUCC float32, and
`pyspharm-cpp` float32 and float64. Each operation has its own zero-based
linear time scale; compare bar lengths within an operation and use the printed
values across operations. Smaller is better. Speedup is PySpharm time divided
by `pyspharm-cpp` float32 time.

All five cases were measured in one [GitHub Actions benchmark run](https://github.com/mwyau/pyspharm-cpp/actions/runs/40109072587)
on October 24, 2026, using an AMD EPYC 9V45 runner with four allocated CPUs,
AVX-512 (`x86-64-v4`) selected, CPython 3.12.15, NumPy 1.26.4,
PySpharm 1.0.9, and DUCC 0.42.0.
The 2.5° one-thread case uses one DUCC thread; the other four cases use
four DUCC threads. PySpharm is serial in both configurations.
Each median uses 7 timed repetitions after 2 warmups, without CPU pinning.
Float64 measurements use casts of the float32 input fields.

## Clenshaw–Curtis · 2.5° · 144 × 73 (one thread)

![Clenshaw–Curtis · 2.5° · 144 × 73 (one thread) SHT benchmark timings](timings-cc-2p5-1t.svg)

| Operation | PySpharm f32 | pyspharm-cpp f32 | DUCC f32 | pyspharm-cpp f64 | Speedup vs PySpharm |
| --- | ---: | ---: | ---: | ---: | ---: |
| Scalar analysis | 0.068 ms | 0.218 ms | 0.219 ms | 0.228 ms | 0.31× |
| Scalar synthesis | 0.031 ms | 0.080 ms | 0.053 ms | 0.061 ms | 0.38× |
| Wind analysis | 0.179 ms | 0.449 ms | 0.424 ms | 0.461 ms | 0.40× |

## Clenshaw–Curtis · 2.5° · 144 × 73 (four threads)

![Clenshaw–Curtis · 2.5° · 144 × 73 (four threads) SHT benchmark timings](timings-cc-2p5.svg)

| Operation | PySpharm f32 | pyspharm-cpp f32 | DUCC f32 | pyspharm-cpp f64 | Speedup vs PySpharm |
| --- | ---: | ---: | ---: | ---: | ---: |
| Scalar analysis | 0.068 ms | 0.211 ms | 0.237 ms | 0.220 ms | 0.32× |
| Scalar synthesis | 0.029 ms | 0.070 ms | 0.060 ms | 0.064 ms | 0.42× |
| Wind analysis | 0.182 ms | 0.427 ms | 0.387 ms | 0.442 ms | 0.42× |

## Clenshaw–Curtis · 0.25° · 1440 × 721 (four threads)

![Clenshaw–Curtis · 0.25° · 1440 × 721 (four threads) SHT benchmark timings](timings-cc-0p25.svg)

| Operation | PySpharm f32 | pyspharm-cpp f32 | DUCC f32 | pyspharm-cpp f64 | Speedup vs PySpharm |
| --- | ---: | ---: | ---: | ---: | ---: |
| Scalar analysis | 95.925 ms | 12.349 ms | 13.911 ms | 12.257 ms | 7.77× |
| Scalar synthesis | 77.371 ms | 3.910 ms | 3.866 ms | 3.896 ms | 19.79× |
| Wind analysis | 243.901 ms | 31.281 ms | 28.891 ms | 30.546 ms | 7.80× |

## Gaussian · 128 latitude rings · 256 × 128 (four threads)

![Gaussian · 128 latitude rings · 256 × 128 (four threads) SHT benchmark timings](timings-gl-128.svg)

| Operation | PySpharm f32 | pyspharm-cpp f32 | DUCC f32 | pyspharm-cpp f64 | Speedup vs PySpharm |
| --- | ---: | ---: | ---: | ---: | ---: |
| Scalar analysis | 0.624 ms | 0.136 ms | 0.126 ms | 0.136 ms | 4.58× |
| Scalar synthesis | 0.109 ms | 0.127 ms | 0.104 ms | 0.130 ms | 0.86× |
| Wind analysis | 3.242 ms | 0.262 ms | 0.213 ms | 0.261 ms | 12.38× |

## Gaussian · 640 latitude rings · 1280 × 640 (four threads)

![Gaussian · 640 latitude rings · 1280 × 640 (four threads) SHT benchmark timings](timings-gl-640.svg)

| Operation | PySpharm f32 | pyspharm-cpp f32 | DUCC f32 | pyspharm-cpp f64 | Speedup vs PySpharm |
| --- | ---: | ---: | ---: | ---: | ---: |
| Scalar analysis | 183.330 ms | 3.165 ms | 3.021 ms | 3.280 ms | 57.92× |
| Scalar synthesis | 45.608 ms | 3.020 ms | 2.734 ms | 3.047 ms | 15.10× |
| Wind analysis | 620.010 ms | 7.681 ms | 6.158 ms | 11.640 ms | 80.72× |

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
