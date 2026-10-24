<!-- SPDX-License-Identifier: GPL-2.0-or-later -->

# Benchmarks

The [benchmark report](https://github.com/mwyau/pyspharm-cpp/blob/main/benchmarks/README.md) compares scalar and wind transforms from PySpharm, `spharm_cpp`, and ducc0 on five Clenshaw–Curtis and Gaussian grids. It includes numerical tables, figures, the measurement environment, and reproduction commands. The [measurement CSV](https://github.com/mwyau/pyspharm-cpp/blob/main/benchmarks/measurements.csv) records the results used for the figures.

The benchmark uses matching spectral truncation for each grid. Synthesis coefficients are prepared before timing, and timing values are medians from repeated measurements after warmup.

The direct ducc0 baseline measures spin-0 and spin-1 transforms. `spharm_cpp` uses ducc0 transforms with additional coefficient conversion for PySpharm-compatible outputs, especially vorticity and divergence. Therefore, the direct spin-1 DUCC result is a transform baseline, not an end-to-end timing of the same wind API.

Performance depends on grid resolution, thread count, processor, and precision. The figures show individual measurements rather than performance guarantees.
