<!-- SPDX-License-Identifier: GPL-2.0-or-later -->

# Benchmarks and reference comparisons

The benchmark writes one row per backend, precision, operation and resolution.
It reports warmed median and p90 latency, gridpoint throughput, result size,
Python-traced peak allocation, and the requested DUCC thread count. The traced
peak does not include all native C++ workspace allocations. Run one process for
each thread count; DUCC reads its thread environment at the first transform.
The candidate also reports `legendre_many` and `specintrp_many` over the
requested field count, alongside scalar C++ calls repeated from Python for the
same points. `legendre_many` applies the scaled degree/order recurrence to
independent latitude lanes. `specintrp_many` vectorizes phase arithmetic and
coefficient accumulation across points. It reduces each longitude modulo 360
degrees, advances sine and cosine phase state per lane, and re-anchors with
direct scalar `sin`/`cos` every 32 orders. The phase vectors feed the existing
SIMD coefficient arithmetic. Incomplete batches use a scalar tail. Compare
paired rows at multiple batch sizes: a batch smaller than the SIMD width, or a
short spectral truncation, may not amortize setup and can be slower than the
scalar path.

`getgeodesicpts(m)` returns one ordered point set of size `10*(m-1)^2+2`, so
there is no regular collection of same-shaped independent outputs to pass as a
SIMD batch. The construction keeps the reference's ordered float32 coordinate
arithmetic and scalar `atan2` conversions; replacing those with vector math can
change point values and ordering. The utility benchmark compares this single
native call with PySpharm's reference call. Keep the measured output order and
reference accuracy in view when evaluating any future vector path.

For transforms, the candidate reports batched scalar analysis/synthesis beside
Python loops that call the same scalar API once per field. Those loops use
strided field views and stack their outputs, while the batched API allocates a
single result array. Both paths call DUCC for every field today; these rows
measure the Python dispatch and repeated allocation/conversion overhead before
considering any change to that loop. The benchmark checks the two paths agree
before timing them.

With locked development tools installed:

```sh
uv sync --locked --python 3.14 --all-groups
uv run --locked python benchmarks/benchmark.py \
  --sizes 256x128 512x256 --fields 8 --threads 1 --output bench-threads-1.csv
uv run --locked python benchmarks/benchmark.py \
  --sizes 256x128 512x256 --fields 8 --threads 4 --output bench-threads-4.csv
```

For focused utility comparisons, use the same command with one resolution,
`--fields 8`, both precisions, at least three warmups, and enough repetitions
for stable medians. The `*_scalar_batch` and `*_many` rows return the same
shapes and isolate the end-to-end difference between repeated scalar bindings
and a single batch call. Re-run with batch sizes below and above the active
SIMD width to identify where batching helps; do not treat a small smoke run as
a performance claim.

For comparisons, install PySpharm 1.0.9 from the separately pinned
`jswhit/pyspharm` checkout and `ducc0==0.42.0` into the comparison environment,
then install the candidate wheel into that same environment. They use distinct
module names (`spharm` and `spharm_cpp`). Run:

```sh
uv run --python /path/to/comparison-env/bin/python \
  benchmarks/compare_pyspharm.py
uv run --python /path/to/comparison-env/bin/python \
  benchmarks/benchmark.py --threads 1 --output comparison.csv
```

PySpharm and the candidate are compared in float32 because the original public
package has no float64 mode. The benchmark also measures float64 for the
candidate and DUCC. The direct DUCC rows use its public `sht` API and are not a
comparison of PySpharm coefficient conventions.

For build-profile comparisons, build portable, native and multiarch wheels with
`uv run --locked --no-sync python -m build --wheel --no-isolation
--config-setting cmake.define.SPHARM_CPP_BUILD_PROFILE=PROFILE`, then install
each wheel into a separate environment and save a separate CSV. Label
native-host measurements with the host CPU; they are not portable results.
`pyspharm-rs` is optional and is not included unless its package is already
installed with a known API.

## Recorded results

`results/` contains warmed runs from Linux x86-64, CPython 3.15.0rc2, on a host
with AVX2 available. Runs used 3 warmups, 10 measured repetitions, 8 fields,
and 256x128 and 512x256 Gaussian grids. The reference runs compare PySpharm
1.0.9, direct DUCC 0.42.0, and this package. `pyspharm-rs` was not installed.
The profile runs compare separately installed portable, native, and multiarch
wheels. They are single-host measurements, not portable performance claims.

At 512x256 in float32 with one thread, median times in milliseconds were:

| Backend/profile | Scalar analysis | Wind analysis | 8-field scalar analysis |
| --- | ---: | ---: | ---: |
| Portable wheel | 1.731 | 5.532 | 14.220 |
| Native wheel | 1.116 | 2.858 | 9.450 |
| Multiarch wheel | 1.105 | 2.866 | 9.385 |
| Direct DUCC 0.42.0 | 1.105 | 2.690 | 9.130 |
| PySpharm 1.0.9 | 28.545 | 20.102 | 225.174 |

The multiarch wheel selected x86-64-v3 on this host. With four threads, the
same three operations measured 0.371, 0.927, and 3.022 ms. CSV files include
median and p90 latency for every benchmarked operation, both precisions, and
both resolutions. Python-traced allocation excludes native C++ workspaces.

The recorded profile wheels used GCC 15.2.0 and locked CMake 4.4.4. DUCC
translation units use `-O3 -fno-fast-math -fno-finite-math-only`; portable and
multiarch v1 use `-march=x86-64`, native uses `-march=native`, and multiarch v3
and v4 use their corresponding `-march=x86-64-v3` and `-march=x86-64-v4`
flags. The small baseline driver uses `-O2 -march=x86-64`; interprocedural
optimization is disabled. Python-traced peak allocation is the only recorded
memory metric. Host load and system-wide memory pressure were not sampled
during the initial timing runs, so those results remain preliminary.

## Batch and profile validation run

`validation-2026-10-23-{portable,native,multiarch}-256x128.csv` contains a
fresh 256x128 Gaussian-grid run from each separately installed CPython 3.15
wheel. `batch-sweep-2026-10-23-multiarch-256x128.csv` contains utility and
multifield rows for batches of 1, 2, 4, 8 and 16. Each process used one DUCC
thread, 5 warmups and 31 measured repetitions. The run sampled host load at
4.41/5.53/4.85 on 32 logical CPUs and 68 GiB available memory; it was not CPU
pinned. The selected multiarch profile was x86-64-v3 (AVX2). GCC 15.2.0 and
CMake 4.4.4 built the wheels using the profile flags above. These are
single-host measurements, not acceptance thresholds or portable claims.

With 8 points and truncation 127, median milliseconds were:

| Profile | Scalar `legendre` calls | `legendre_many` | Scalar `specintrp` calls | `specintrp_many` | `getgeodesicpts(32)` |
| --- | ---: | ---: | ---: | ---: | ---: |
| Portable | 0.432 | 0.195 | 0.064 | 0.041 | 0.634 |
| Native | 0.423 | 0.100 | 0.066 | 0.029 | 0.603 |
| Multiarch v3 | 0.426 | 0.098 | 0.065 | 0.029 | 0.600 |

The multiarch utility batch sweep measured the scalar-call/batch speed ratio:

| Batch size | `legendre_many` | `specintrp_many` |
| ---: | ---: | ---: |
| 1 | 0.99x | 1.00x |
| 2 | 1.01x | 1.03x |
| 4 | 4.24x | 2.30x |
| 8 | 4.34x | 2.23x |
| 16 | 4.38x | 2.22x |

Batches below the 4-lane AVX2 width do not show a useful SIMD gain. `legendre`
uses a degree/order recurrence for one latitude, so its dependency chain is
serial; independent latitudes are the SIMD dimension. The table above records
the pre-recurrence interpolation baseline, which evaluated sine and cosine at
every order. The phase recurrence results below supersede those interpolation
timings. The ordered geodesic generator remains scalar to preserve float32
point values and ordering; one output set does not provide a same-shaped batch
dimension.

On the same multiarch wheel, 8-field transform time stayed near the scalar
loop time: float64 analysis was 2.243 ms for scalar calls and 2.221 ms for the
batch API; synthesis was 2.580 and 2.585 ms. Both call DUCC once per field, so
the measurement does not justify changing that loop. Python-traced peak
allocation fell from 1,064,552 to 880 bytes for analysis and from 2,105,336 to
920 bytes for synthesis; this tracer excludes native workspace memory and is
not a process-wide memory measurement.

The same multiarch run compared float32 transforms with installed DUCC 0.42.0
on identical 256x128 grids and one thread. Candidate/DUCC median milliseconds
were 0.243/0.234 for scalar analysis, 0.283/0.271 for scalar synthesis,
0.537/0.495 for wind analysis, and 0.691/0.643 for wind synthesis. PySpharm
was not installed in these fresh wheel environments; the earlier
`reference-threads-1.csv` retains the PySpharm 1.0.9 comparison. No current
PySpharm speed claim is based on this run.

## Scaled Legendre regression and batch measurements

The `legendre-scale-2026-10-23-*.csv` files record the first benchmark run
after adding exponent tracking to the degree/order recurrence and a no-op
scaling fast path to the SIMD batch kernel. The profile files use 256x128
Gaussian grids (`lmax=127`), eight latitudes, both precisions, one thread, five
warmups, and 31 measured repetitions. The multiarch batch files sweep batch
sizes 1, 2, 4, 8 and 16 over 32x16 (`lmax=15`) and 256x128 (`lmax=127`) grids
in float64. Each row compares stacked scalar calls with one `*_many` call and
checks their outputs before timing. The same files include scalar Legendre,
scalar interpolation, geodesic point generation, and multi-field transform
measurements.

The run used CPython 3.15.0rc2, GCC 15.2.0, CMake 4.4.4, and the host's AVX2
profile. The host was not pinned; immediately before the run it had 32 logical
CPUs, load averages 8.22/6.21/5.31, about 57 GiB available memory, and no
memory pressure. CPU pressure was low, though the system had about 7.9 GiB of
8 GiB swap in use. Treat these single-host medians as exploratory rather than
definitive.

At `lmax=127` and eight latitudes, median milliseconds and the scalar/batch
ratio were:

| Profile | Precision | Scalar Legendre calls | `legendre_many` | Ratio | Scalar interpolation calls | `specintrp_many` | Ratio |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: |
| Portable | float64 | 0.4467 | 0.2415 | 1.85x | 0.0664 | 0.0426 | 1.56x |
| Native | float64 | 0.4560 | 0.1537 | 2.97x | 0.0683 | 0.0343 | 1.99x |
| Multiarch v3 | float64 | 0.4339 | 0.1523 | 2.85x | 0.0672 | 0.0300 | 2.24x |

On multiarch v3, the batch sweep at `lmax=127` measured `legendre_many`
speedups of 1.11x, 1.01x, 2.81x, 2.85x and 2.89x for batch sizes 1, 2, 4, 8
and 16. `specintrp_many` measured 1.01x, 1.02x, 2.19x, 2.24x and 2.23x. At
`lmax=15`, batches of four or more helped both utilities, while one-point
interpolation was slightly slower through the batch API. Small batches near
the vector width therefore do not guarantee a speedup.

The 256x128 multiarch transform timings remained near the repeated scalar-loop
times: eight-field analysis was 2.245 ms batched versus 2.265 ms by scalar
calls, and synthesis was 2.604 versus 2.583 ms. The field loop still calls
DUCC once per field, so no transform batching change is justified by this
measurement. The single-call geodesic utility measured 0.615 ms on multiarch,
0.635 ms native, and 0.650 ms portable; it remains scalar to preserve ordering
and reference rounding.

## SIMD interpolation phase experiment

`results/specintrp-kernel-sweep-2026-10-23.csv` compares a native scalar loop
with six isolated C++ SIMD kernels using the same coefficients and real
`legendre_many` outputs: current coefficient-wise phase arithmetic (A),
per-order real/imaginary accumulation (B), phase recurrence with 32- and
64-order re-anchoring (C), and each combination (D). The sweep covers both
precisions, batches 1–64 at `lmax=31` and 127, and practical batch sizes at
`lmax=511`, 1023 and 2048. The separate confirmation file uses 31 repetitions
for selected representative cases. Both record median/p90 latency, memory
footprints, and error against the direct scalar implementation and a long
double phase reference.

Runs used one process per case pinned to CPU 8 on an AMD Ryzen 9 5950X, with
four double-precision SIMD lanes, GCC 15.2.0, and `-O3 -fno-fast-math
-march=x86-64-v3`. The broad sweep used five warmups and 31 repetitions through
`lmax=127`, then 11, 7 and 5 repetitions at `lmax=511`, 1023 and 2048. The
confirmation run used five warmups and 31 repetitions. CPU pressure was low at
the start and memory pressure was zero; the host's swap was full. Hardware
performance counters were unavailable because the host sets
`perf_event_paranoid=4`, so this evidence comes from paired wall-clock kernel
measurements rather than sampled instruction or cache counters.

The combined variant D with 32-order re-anchoring was retained. At `lmax=127`
and eight points, its native kernel was about 2.5x faster than A; at `lmax=511`
and 1023 it was about 1.2x and 1.1x faster. For `lmax=2048`, results ranged
from no useful gain to about 1.2x, consistent with the larger Legendre gather
traffic dominating phase work. Transposing or copying the full Legendre array
was not adopted: counters were unavailable, and the measurements do not
justify adding a large temporary array.

The installed-wheel comparison in
`results/specintrp-api-comparison-2026-10-23.csv` ran before and after the
change on identical inputs, with both wheels installed in separate CPython
3.15.0rc2 environments on the same pinned AVX2 host. Median `specintrp_many`
times in float64 were:

| Truncation | Batch | Previous wheel | Phase recurrence | Improvement |
| ---: | ---: | ---: | ---: | ---: |
| 31 | 8 | 5.009 µs | 2.574 µs | 1.95x |
| 127 | 8 | 31.418 µs | 20.138 µs | 1.56x |
| 127 | 16 | 66.653 µs | 36.618 µs | 1.82x |
| 511 | 8 | 324.790 µs | 245.270 µs | 1.32x |
| 1023 | 8 | 1294.945 µs | 1101.037 µs | 1.18x |
| 2048 | 8 | 5808.712 µs | 5832.839 µs | 1.00x |

At `lmax=127`, batch 8, the new batch call is 3.7x faster than repeated
Python scalar calls end to end. At batch sizes below the four-lane v3 width,
the kernel uses its scalar tail and showed little or inconsistent speedup.
The 80-bit long-double reference reported a maximum global relative error of
1.3e-13 across the measured combined-kernel cases, including longitudes
shifted by large multiples of 360 degrees. The unchanged scalar `specintrp`
API remains the compatibility baseline for ordinary longitudes; reduction in
`specintrp_many` keeps periodic behavior stable at high order.

### Hosted high-order phase regression

Actions runs `40006774214` and `40007236641` exposed float64 errors in the
8-point, degree-2048 periodic-phase case. Run `40012569947` later failed the
degree-8 interpolation test with v4 active, but its `longitude_storage[::2]`
input contained four points. The regression using eight distinct lane inputs
was added after the fix and passed on v4 in runs `40019503881` and
`40020759865`. The degree-8 failure still demonstrates that the defect was not
limited to long phase-recurrence drift. The test reference uses 80-digit
`Decimal` arithmetic to multiply the exact binary64 longitude by order before
reducing modulo 360. It checks both precisions, mixed/real-only/imaginary-only
coefficients, and periodic equivalents. Float64 absolute tolerance remains
`1e-20`.

The degree-2048 v4 output `-1.5641369526414072e-10` matched its real-only
contribution, while the full complex reference was
`-1.0863811359287616e-10`. The compiler mechanism was subsequently located in
GCC 13.3 disassembly of the full pre-fix `spectral_interpolate_many` binding
under `-O3 -march=x86-64-v4`: the scalar `m > 0` guard became a `k1=0x01`
mask on the imaginary `vmulpd`, updating only lane 0. The unconditional
zero-factor update in `4313d91` removes that mask and preserves all eight v4
lanes. The compact four-variant harness did not preserve the binding's code
generation context and did not reproduce the issue; it has been removed from
CI and the repository. Existing installed-package regressions passed on v4 in
runs `40019503881` and `40020759865`. The exact assembly excerpt and
configuration are recorded in
`results/avx512-interpolation-failure-evidence-2026-10-23.txt`. This is
evidence for the tested GCC 13.3 binding, not a claim of a general GCC,
libstdc++ or DUCC defect. No upstream issue was filed and DUCC remains
unmodified.

The fix in `4313d91` uses DUCC `native_simd<double>` at full width in every
profile, including eight lanes for v4. Instead of conditionally skipping the
imaginary vector update for `m == 0`, it sets a scalar imaginary factor to zero
for that order and always performs the vector operation. The output vector is
stored with DUCC `storeu`. The scalar `specintrp` path, phase recurrence,
coefficient order and v1/v3 SIMD widths are unchanged. A regression now uses a
single nonzero complex coefficient with distinct data in all eight lanes at
degree 8, plus mixed, imaginary-only and real-only coefficients at degree
2048.

Actions run `40019503881` passed all nine jobs. Its CPython 3.11 test runner
reported `active_profile: x86-64-v4` and passed all 61 tests, including the two
interpolation regressions. Its multiarch wheel passed all 61 tests and the
focused regressions with caps 1 and 3; its runner lacked v4 and explicitly
reported that cap 4 was untested. The focused profile sweep in run
`40020759865` passed: one v4-capable host ran all three profile caps, confirmed
each with `cpu_info()`, and passed four parametrized regression cases per cap.
The separate diagnostic workflow was removed; primary CI retains numerical
package tests and CPU-profile reporting.

The post-correction local API matrix passed 240 sparse complex-reference
checks under each executable profile (v1 and v3), spanning five truncations,
eight batch sizes, both precisions and real/imaginary/mixed coefficients. It
also checked read-only arrays, strided longitudes and scalar-vs-batch results.
A fresh GCC 15.2 ASan/UBSan multiarch build passed all 61 tests under v1 and
v3. The local host does not execute v4; the hosted profile runs above provide
the v4 runtime evidence. The exact matrix and sanitizer conditions are in
`results/specintrp-avx512-final-validation-2026-10-23.txt`.

`results/specintrp-ci-fix-v3-2026-10-23.csv` remains an exploratory pre-fix
comparison of the phase recurrence and a scalar per-lane phase-state candidate
at `lmax=127`, 511 and 1023. Each row used five warmups and 41 calls in
separate CPython 3.15.0rc2 environments forced to v3. The machine was not CPU
pinned and load was not sampled, so small differences are not meaningful. The
v4 correction does not change the phase algorithm or v1/v3 widths. The paired
post-correction measurements are recorded below.

### v3 performance after the AVX-512 correction

`results/specintrp-avx512-fix-v3-2026-10-23.csv` compares the exact `ba9a133`
safe-checkpoint source with the corrected `4313d91` source. Both multiarch
wheels were built with GCC 15.2 `-O3 -fno-fast-math -fno-finite-math-only`
using the same pinned DUCC revision, Python 3.15.0rc2, NumPy 2.5.4 and
`x86-64-v3` dispatch. Runs used a Ryzen 9 5950X host, two DUCC threads, CPU 18
affinity, seven warmups and 41 measured calls. Each version used the same
coefficient and Legendre arrays for `lmax=127, 511, 1023` and batches 1, 4, 8
and 16. The CSV records median/p90 latency, scalar-loop comparisons, numerical
differences from repeated scalar calls and input-memory sizes.

For batch size 8, the safe checkpoint versus corrected median (p90) times
were 17.924 (18.084) versus 17.533 (17.663) µs at `lmax=127`, 241.404
(247.445) versus 242.716 (250.471) µs at `lmax=511`, and 1135.945 (1171.532)
versus 1151.244 (1182.843) µs at `lmax=1023`. The batch-8 median changes
ranged from 2.2% faster to 1.35% slower; across all 12 truncation/batch pairs,
the largest median increase was 1.58% (`lmax=1023`, batch 4). The
repeated-scalar speedups were 3.83x/3.93x, 3.94x/3.91x and 3.41x/3.35x
respectively. At all 12 pairs, the v3 Legendre arrays and interpolation outputs
from the two builds were bitwise identical. Per-pair results and SHA-256
identities of both wheel artifacts are in
`results/specintrp-avx512-fix-v3-equality-2026-10-23.txt`. No v3 SIMD
regression is apparent, and no speedup from the correctness fix is claimed.

An unrelated CPU test was active elsewhere on the host during these runs. The
benchmark was pinned to CPU 18, which was mostly idle, but the host was not
otherwise isolated. Treat the small timing deltas as noise and use these data
as a regression check rather than definitive performance claims. The local CPU
does not support v4, so this table contains no AVX-512 timing. The hosted
multiarch wheel has passed degree-8 and degree-2048 complex regressions under
all three profile caps on one v4-capable runner in Actions run `40020759865`.

### Linked stable ABI probe

`results/stable-abi-2026-10-23-cp312-linux-x86_64.csv` compares the same
portable C++ implementation built as a normal `cp312-cp312` wheel and as a
linked `cp312-abi3` wheel. Both were built with GCC 15.2, CPython 3.12.13 and
the pinned DUCC revision, then installed into separate clean environments with
NumPy 2.5.4. The benchmark used a Gaussian `256 x 128` grid, eight fields,
float64, one DUCC thread, five warmups and 25 measured calls per operation.
Three runs per wheel were alternated on an AMD Ryzen 9 5950X; the host was not
CPU-pinned or otherwise isolated.

Across the three run medians, abi3 latency ranged from 8.7% lower to 3.6%
higher depending on operation. The transform operations showed no consistent
shift, and the small deltas are within the limitations of this unpinned host.
This is evidence against a large local regression, not a general performance
claim. The source ABI matrix and installed-wheel tests are checked separately
from these timings.
