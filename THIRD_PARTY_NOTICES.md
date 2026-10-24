<!-- SPDX-License-Identifier: GPL-2.0-or-later -->

# Third-party source and license notices

## DUCC

This project builds DUCC C++ sources from the unmodified Git submodule at
`third_party/ducc`, upstream repository
<https://gitlab.mpcdf.mpg.de/mtr/ducc.git>. The submodule is pinned to release
tag `ducc0_0_42_0`, commit
`42ce8ab5c5d79e0e2fd9af4a09c2ea4693520e1b`.

The build compiles DUCC's SHT/FFT and required support translation units from
`third_party/ducc/src/ducc0/` and reuses
`third_party/ducc/python/multiarch.cc` and `multiarch.h` for CPU detection and
dispatch. The source distribution includes the complete pinned
`third_party/ducc/src/**` tree so transitive headers retain their original
paths, the two multiarch files, and `third_party/ducc/LICENSE`. It excludes
DUCC Git metadata. A wheel built from the sdist therefore needs ordinary Python
build dependencies but does not retrieve DUCC or require Git metadata.

DUCC source code is licensed under GNU GPL version 2 or later. Files in its
FFT component and internal dependencies may additionally be used under the
3-clause BSD license when their headers provide both grants. Preserve the
upstream license file and per-file notices when redistributing source or
binaries. The top-level `LICENSE` supplies the GPL license text used by this
project's package metadata.

## PySpharm behavioral reference

PySpharm is used only as an API and numerical behavior reference. It is not
imported or required at runtime. No PySpharm source is incorporated in this
package.
