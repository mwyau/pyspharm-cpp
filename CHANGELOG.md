<!-- SPDX-License-Identifier: GPL-2.0-or-later -->

# Changelog

## Unreleased

- Corrected mathematical equations in the documentation and recorded the
  tested API differences, numerical limits, and workflow configuration.
- Reorganized native code into focused C++ binding, transform and numerical
  modules while retaining the pinned DUCC implementation.
- Added scalar and batched numerical utilities, free-threaded CPython build
  support, typed signatures, and expanded numerical and concurrency tests.
- Added runner architecture checks, reproducible wheel metadata validation,
  development checks, documentation and native coverage configuration.
- Added standard-ABI wheel builds: CPython 3.11 keeps a `cp311` wheel, CPython
  3.12 and newer use a shared `cp312-abi3` wheel, and free-threaded versions
  remain separately tagged.
- No package has been published and no release tag has been created.
