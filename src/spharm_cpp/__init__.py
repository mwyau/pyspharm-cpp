# SPDX-FileCopyrightText: 2026 Albert M. W. Yau
# SPDX-License-Identifier: GPL-2.0-or-later
"""C++17 spherical harmonic transforms with a PySpharm-compatible API."""

from importlib.metadata import version

from ._spharm_cpp import (
    Spharmt,
    cpu_info,
    gaussian_lats_wts,
    getgeodesicpts,
    getspecindx,
    legendre,
    regrid,
    specintrp,
)

__version__ = version("pyspharm-cpp")

__all__ = [
    "Spharmt",
    "cpu_info",
    "gaussian_lats_wts",
    "getgeodesicpts",
    "getspecindx",
    "legendre",
    "regrid",
    "specintrp",
]
