# Benchmarks

The [benchmark report](https://github.com/mwyau/pyspharm-cpp/blob/main/benchmarks/README.md) compares scalar and wind transforms from `pyspharm==1.0.9`, `spharm_cpp`, and ducc0 on five Clenshaw–Curtis and Gaussian grids. It includes numerical tables, figures, the measurement environment, and reproduction commands. The [measurement CSV](https://github.com/mwyau/pyspharm-cpp/blob/main/benchmarks/measurements.csv) records the results used for the figures.

The benchmark uses matching spectral truncation for each grid. Synthesis coefficients are prepared before timing, and timing values are medians from repeated measurements after warmup.

The ducc0 baseline measures spin-0 and spin-1 transforms. The `spharm_cpp` wind operations also convert between wind fields and spectral coefficients for relative vorticity and divergence.

The four multithreaded cases use eight threads. The 2.5° grid has both a one-thread case and an eight-thread case. The benchmark script checks process affinity and cgroup quota before measuring.

The manual GitHub Actions workflow currently targets a four-CPU `ubuntu-24.04` runner, so its preflight stops before measurement. Configure an eight-CPU runner before using the workflow to reproduce the suite.

Performance depends on grid resolution, thread count, processor, and precision. The figures show individual measurements rather than performance guarantees.
