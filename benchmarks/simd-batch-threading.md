<!-- SPDX-License-Identifier: GPL-2.0-or-later -->

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

The auxiliary timing pass used the default single-thread setting. It measures
the SIMD changes, not scaling from the kernel-level parallel loops added later;
thread-count equivalence for those loops is covered by
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
The candidate assigned independent fields to DUCC workers while each field
transform used one inner thread. Wind synthesis uses DUCC's transform-level
threading and is not assigned an outer field schedule.
