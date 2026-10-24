<!-- SPDX-License-Identifier: GPL-2.0-or-later -->

# Benchmarks

These figures show measured **float32 median execution times** in milliseconds.
Each operation has its own zero-based linear time scale; compare bar lengths
within an operation, and use the printed values across operations. Smaller is
better. The speedup column is PySpharm time divided by pyspharm-cpp time.

Measured on a four-CPU AMD EPYC Linux runner (CPython 3.12, NumPy 1.26.4):
serial PySpharm 1.0.9 versus four-thread DUCC 0.42.0 and pyspharm-cpp
0.1.0.dev0, with AVX-512 (`x86-64-v4`) selected. Medians use seven timed
repetitions after two warmups, without CPU pinning.

The fifth case repeats the 2.5° regular grid with one DUCC thread, on the
same runner as the four-thread measurements. The figures report the two
thread configurations separately.

## Clenshaw–Curtis · 2.5° · 144 × 73 (four threads)

![Timings for Clenshaw–Curtis · 2.5° · 144 × 73](timings-cc-2p5.svg)

| Operation        | PySpharm | pyspharm-cpp |     DUCC | Speedup vs PySpharm |
| ---------------- | -------: | -----------: | -------: | ------------------: |
| Scalar analysis  | 0.069 ms |     0.209 ms | 0.222 ms |               0.33× |
| Scalar synthesis | 0.029 ms |     0.067 ms | 0.063 ms |               0.44× |
| Wind analysis    | 0.178 ms |     0.385 ms | 0.381 ms |               0.46× |

## Clenshaw–Curtis · 0.25° · 1440 × 721

![Timings for Clenshaw–Curtis · 0.25° · 1440 × 721](timings-cc-0p25.svg)

| Operation        |   PySpharm | pyspharm-cpp |      DUCC | Speedup vs PySpharm |
| ---------------- | ---------: | -----------: | --------: | ------------------: |
| Scalar analysis  |  95.027 ms |    12.638 ms | 14.260 ms |               7.52× |
| Scalar synthesis |  76.561 ms |     4.047 ms |  3.875 ms |              18.92× |
| Wind analysis    | 244.674 ms |    27.712 ms | 29.122 ms |               8.83× |

## Gaussian · 128 latitude rings · 256 × 128

![Timings for Gaussian · 128 latitude rings · 256 × 128](timings-gl-128.svg)

| Operation        | PySpharm | pyspharm-cpp |     DUCC | Speedup vs PySpharm |
| ---------------- | -------: | -----------: | -------: | ------------------: |
| Scalar analysis  | 0.634 ms |     0.131 ms | 0.128 ms |               4.82× |
| Scalar synthesis | 0.109 ms |     0.128 ms | 0.118 ms |               0.85× |
| Wind analysis    | 1.161 ms |     0.265 ms | 0.227 ms |               4.38× |

## Gaussian · 640 latitude rings · 1280 × 640

![Timings for Gaussian · 640 latitude rings · 1280 × 640](timings-gl-640.svg)

| Operation        |   PySpharm | pyspharm-cpp |     DUCC | Speedup vs PySpharm |
| ---------------- | ---------: | -----------: | -------: | ------------------: |
| Scalar analysis  | 184.648 ms |     3.450 ms | 3.199 ms |              53.53× |
| Scalar synthesis |  47.849 ms |     3.065 ms | 2.964 ms |              15.61× |
| Wind analysis    | 612.969 ms |     7.590 ms | 6.405 ms |              80.76× |

## Clenshaw–Curtis · 2.5° · 144 × 73 (one thread)

The additional one-thread figure is produced from the paired runner
measurements. It is not included in the historical four-thread snapshot.

## Reproduce

The benchmarks compare the same truncation on each backend. The regular-grid
analysis limit is `nlat-2`, Gaussian is `nlat-1`, and even-longitude
vector Nyquist modes are excluded. The benchmark uses seeded float32 inputs.

Install the project and the optional `pyspharm-syl==1.2.1` and
`ducc0==0.42.0` references into a CPython 3.12 environment with NumPy 1.26.4.
From the repository root, run:

```sh
python benchmarks/benchmark.py --threads 0 --warmups 2 --repeats 7 \
  --output /path/to/measurements.csv
# To reproduce the five charts, also run:
python benchmarks/benchmark.py --cases cc-2p5-1t --threads 1 \
  --warmups 2 --repeats 7 --output /path/to/one-thread.csv
python benchmarks/plot.py \
  --csv /path/to/measurements.csv /path/to/one-thread.csv \
  --output-dir benchmarks
```

`--threads 0` selects up to eight available CPUs, respecting Linux
affinity and cgroup quota. Use an explicit thread count for controlled runs.
`python benchmarks/plot.py` without a CSV regenerates the historical
four-thread figures from the embedded medians until all five measurements
from a matched runner are recorded. Raw benchmark CSVs are not
tracked. Independent numerical validation is in [`tests/parity/`](../tests/parity).
