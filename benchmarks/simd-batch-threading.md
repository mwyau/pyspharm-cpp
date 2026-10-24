# SIMD and batched transform measurements

The SIMD comparison measured `legendre_many` and `specintrp_many` before and
after the SIMD changes, comparing baseline commit `884aaae` with candidate
commit `2e40674`. It covered 44 combinations of operation, precision, point
count, and truncation. Each timing is the median of 7 repetitions after 2
warmups. Inputs and outputs were prepared outside the timed region, and the
working set was capped at 256 MiB.

The runs used an AMD Ryzen 9 5950X, Linux x86-64, CPython 3.15.1, NumPy 2.5.4,
the `x86-64-v3` build profile, and CPU affinity 0–31. The table gives selected
128-point timings at two truncations. Speedup is baseline time divided by
candidate time; values above 1 indicate a faster candidate.

| Operation | Precision | `lmax=64`, baseline → candidate | Speedup | `lmax=256`, baseline → candidate | Speedup |
| --- | --- | ---: | ---: | ---: | ---: |
| `legendre_many` | float32 | 0.662 → 0.471 ms | 1.41× | 10.198 → 11.690 ms | 0.87× |
| `legendre_many` | float64 | 0.636 → 0.591 ms | 1.08× | 23.655 → 23.387 ms | 1.01× |
| `specintrp_many` | float32 | 0.0863 → 0.0527 ms | 1.64× | 1.141 → 0.674 ms | 1.69× |
| `specintrp_many` | float64 | 0.0735 → 0.0784 ms | 0.94× | 1.031 → 1.093 ms | 0.94× |

The measured changes improved float32 spectral interpolation in these cases.
Legendre timing and float64 interpolation varied less consistently with
precision and truncation. Disassembly of the `x86-64-v3` object contains YMM
float32 and float64 arithmetic; this confirms vector instructions are present,
but does not measure their contribution to runtime.

The coefficient-lane interpolation variant was discarded after its medians
were about 4.7× slower for float32 and 6.4× slower for float64. Point-lane SIMD is
used for interpolation. A four-lane Legendre tail also slowed small batches;
the final kernel uses native-width SIMD batches and a scalar tail.

The kernel thread-scaling pass measured the current implementation at 1, 2,
and 4 threads with [`benchmark_auxiliary.py`](benchmark_auxiliary.py). It used
four bounded cases: 8 points at `lmax=16`, 64 at `lmax=64`, 128 at `lmax=256`,
and 16 at `lmax=1024`. Each timing is the median of 7 repetitions after 2
warmups. Inputs were preallocated, and the estimated working set was below
256 MiB. The runs used the same processor and affinity, with CPython 3.12.13,
NumPy 1.26.4, and the `x86-64-v3` profile.

The table reports the median speedup across the four cases per operation and
precision. Speedup is the one-thread time divided by the 2- or 4-thread time.

| Operation | Precision | 2 threads: median (range) | 4 threads: median (range) |
| --- | --- | ---: | ---: |
| `legendre_many` | float32 | 1.25× (1.04–2.04) | 1.24× (1.08–2.66) |
| `legendre_many` | float64 | 1.56× (1.02–1.76) | 2.58× (1.07–2.86) |
| `specintrp_many` | float32 | 1.05× (1.01–2.02) | 1.18× (1.00–2.73) |
| `specintrp_many` | float64 | 1.49× (1.03–1.97) | 1.93× (1.04–3.39) |

The small-low case remained near one-thread latency because it has fewer than
four SIMD batches and bypasses DUCC scheduling. Larger cases showed speedups of
up to 3.39×. Thread-count equivalence for these kernels is covered by
[`test_threading.py`](../tests/test_threading.py).

The batched transform comparison used baseline commit `2e40674`. It covered
regular 64 × 32, Gaussian 256 × 128, and Gaussian 512 × 256 grids; 1, 4, and
16 fields; float32 and float64; four scalar and wind operations; and 1 or 4
DUCC threads. Inputs and outputs were allocated before timing, with 2 warmups
and 7 measured repetitions. The same processor, software versions, build
profile, and affinity were used as for the SIMD comparison.

For the medium and large Gaussian grids with 4 or 16 fields, the median speedup
across 8 cases per operation was:

| Operation | 1 thread: median (range) | 4 threads: median (range) |
| --- | ---: | ---: |
| Scalar analysis | 1.00× (0.96–1.01) | 1.50× (1.32–1.91) |
| Scalar synthesis | 0.97× (0.94–1.03) | 1.34× (1.02–1.77) |
| Wind analysis | 0.99× (0.95–1.11) | 1.45× (1.17–2.26) |

The speedup is baseline time divided by candidate time. Single-thread medians
were close to the baseline; wind analysis ranged from 5% slower to 10% faster.
The candidate assigned fields to DUCC workers for batches of at least four
fields on grids of at least 32,768 pixels, with one DUCC thread per field
transform. Wind synthesis uses DUCC's transform-level threading and is not
assigned an outer field schedule. Its field-batched workspace variant was
reverted after 16-field float64 single-thread timings rose from 5.86 to 9.30 ms
on the medium grid and from 38.5 to 53.3 ms on the large grid.
