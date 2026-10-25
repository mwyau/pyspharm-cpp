# Benchmarks

The [benchmark report](https://github.com/mwyau/pyspharm-cpp/blob/main/benchmarks/README.md) compares scalar and wind transforms from `pyspharm==1.0.9`, `spharm_cpp`, and ducc0 on five Clenshaw–Curtis and Gaussian grids. It includes numerical tables, figures, the measurement environment, and reproduction commands. The [measurement CSV](https://github.com/mwyau/pyspharm-cpp/blob/main/benchmarks/measurements.csv) records the results used for the figures.

The benchmark uses matching spectral truncation for each grid. Synthesis coefficients are prepared before timing, and timing values are medians from repeated measurements after warmup.

The manual workflow builds `pyspharm==1.0.9` from its source distribution on Python 3.11 with NumPy 1.26.4, GNU Fortran, and patched `setuptools==84.1.0`. The setuptools update changes the build tool only; it keeps the reference implementation and Fortran build configuration used for the benchmark. The committed measurements remain the separately recorded results in [`benchmarks/README.md`](https://github.com/mwyau/pyspharm-cpp/blob/main/benchmarks/README.md) and [`measurements.csv`](https://github.com/mwyau/pyspharm-cpp/blob/main/benchmarks/measurements.csv).

The ducc0 baseline measures spin-0 and spin-1 transforms. The `spharm_cpp` wind operations also convert between wind fields and spectral coefficients for relative vorticity and divergence.

The four multithreaded cases use eight threads. The 2.5° grid has both a one-thread case and an eight-thread case. The benchmark script checks process affinity and cgroup quota before measuring.

The manual GitHub Actions workflow requires an AVX-512-capable `ubuntu-24.04` runner and uses the maximum CPUs available to its process, typically four. It benchmarks `x86-64-v1`, `x86-64-v3`, and `x86-64-v4` separately, uploading five figures and CSVs for each profile. These workflow artifacts are distinct from the committed eight-thread results.

Performance depends on grid resolution, thread count, processor, and precision. The figures show individual measurements rather than performance guarantees.
