<!-- SPDX-License-Identifier: GPL-2.0-or-later -->

# References

`pyspharm-cpp` uses DUCC for spherical harmonic transforms and follows PySpharm's Python API, which was built around NCAR SPHEREPACK. The following publications describe these numerical methods and their development.

- Reinecke, M., Belkner, S., and Carron, J. (2023). “Improved cosmic microwave background (de-)lensing using general spherical harmonic transforms.” *Astronomy & Astrophysics*, **678**, A165. [doi:10.1051/0004-6361/202346717](https://doi.org/10.1051/0004-6361/202346717). Describes general spherical harmonic transforms and their implementation in DUCC, with an emphasis on nonuniform sampling.
- Reinecke, M., and Seljebotn, D. S. (2013). “Libsharp – spherical harmonic transforms revisited.” *Astronomy & Astrophysics*, **554**, A112. [doi:10.1051/0004-6361/201321494](https://doi.org/10.1051/0004-6361/201321494). Describes libsharp, from which DUCC's spherical harmonic transform implementation derives.
- Adams, J. C., and Swarztrauber, P. N. (1999). “SPHEREPACK 3.0: A model development facility.” *Monthly Weather Review*, **127**, 1872–1878. [Journal article](https://journals.ametsoc.org/view/journals/mwre/127/8/1520-0493_1999_127_1872_samdf_2.0.co_2.xml). Describes the spherical harmonic transforms and geophysical differential operators implemented in NCAR SPHEREPACK.

For the Python interface used in compatibility comparisons, see [PySpharm by Jeffrey Whitaker](https://github.com/jswhit/pyspharm). DUCC source and license information are available in its [upstream repository](https://gitlab.mpcdf.mpg.de/mtr/ducc) and this project's [third-party notices](https://github.com/mwyau/pyspharm-cpp/blob/main/THIRD_PARTY_NOTICES.md).
