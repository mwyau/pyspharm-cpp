# Benchmarks

The [benchmark report](https://github.com/mwyau/pyspharm-cpp/blob/main/benchmarks/README.md) compares scalar and wind transforms from `pyspharm==1.0.9`, `spharm_cpp`, and ducc0 on five Clenshaw–Curtis and Gaussian grids. It includes numerical tables, figures, the measurement environment, and reproduction commands. The [measurement CSV](https://github.com/mwyau/pyspharm-cpp/blob/main/benchmarks/measurements.csv) records the results used for the figures.

The benchmark uses matching spectral truncation for each grid. Synthesis coefficients are prepared before timing, and timing values are medians from repeated measurements after warmup.

The committed measurements were regenerated with the Python 3.14 conda environment in `.github/conda/pyspharm-reference.yml` and conda-forge PySpharm 1.0.9. They were collected on an AMD Ryzen 9 5900XT under Linux 6.6.87.2 WSL2 with the x86-64-v3 SIMD profile active. The benchmark report records the full environment and reproduction commands.

The ducc0 baseline measures spin-0 and spin-1 transforms. The `spharm_cpp` wind operations also convert between wind fields and spectral coefficients for relative vorticity and divergence.

The four multithreaded cases use eight threads. The 2.5° grid has both a one-thread case and an eight-thread case. The benchmark script checks process affinity and cgroup quota before measuring.

The manual GitHub Actions workflow requires an AVX-512-capable `ubuntu-24.04` runner and uses the maximum CPUs available to its process, typically four. It benchmarks `x86-64-v1`, `x86-64-v3`, and `x86-64-v4` separately, uploading five figures and CSVs for each profile. These workflow artifacts are distinct from the committed eight-thread results.

Performance depends on grid resolution, thread count, processor, and precision. The figures show individual measurements rather than performance guarantees.
