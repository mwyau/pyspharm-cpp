# SPDX-FileCopyrightText: 2026 Albert M. W. Yau
# SPDX-License-Identifier: GPL-2.0-or-later
"""Public PySpharm API inventory enforced by parity tests."""

LEGACY_CONTRACTS = {
    "Spharmt.__init__": "constructor signature, defaults, validation, and dtype",
    "Spharmt.attributes": "public read-only configuration and user attributes",
    "Spharmt.__repr__": "diagnostic representation",
    "Spharmt.grdtospec": "scalar analysis and batch behavior",
    "Spharmt.spectogrd": "scalar synthesis and batch behavior",
    "Spharmt.getvrtdivspec": "wind analysis",
    "Spharmt.getuv": "wind synthesis",
    "Spharmt.getgrad": "scalar gradient",
    "Spharmt.getpsichi": "streamfunction and velocity potential",
    "Spharmt.specsmooth": "spectral smoothing",
    "regrid": "grid conversion with explicit truncation and smoothing",
    "gaussian_lats_wts": "Gaussian latitudes and quadrature weights",
    "getspecindx": "triangular coefficient indices and dtype",
    "getgeodesicpts": "ordered geodesic coordinates and precision",
    "legendre": "normalized associated Legendre values for scalar and batched latitudes",
    "specintrp": "spectral interpolation at a scalar point or across batched points",
}

ADDITIVE_CONTRACTS = {
    "cpu_info": "DUCC CPU metadata schema and value types",
    "Spharmt.dtype": "selected map and coefficient precision",
}

APPROVED_DIFFERENCES = {
    "backend": ("The public spharm_cpp interface uses the compiled DUCC backend."),
    "legendre_storage": (
        "The default is legfunc='computed'; legfunc='stored' raises ValueError."
    ),
    "default_precision": (
        "The default is float64/complex128; dtype=numpy.float32 selects "
        "float32/complex64."
    ),
    "regular_analysis_limit": (
        "Regular-grid analysis stops at nlat-2, preserving DUCC's native "
        "Clenshaw-Curtis limit."
    ),
    "shape_at_zero_truncation": (
        "Rank-one coefficient outputs retain a consistent shape at ntrunc=0."
    ),
    "geodesic_float64_default": (
        "getgeodesicpts defaults to double-precision geometry and float64 "
        "coordinates; dtype=numpy.float32 casts the computed coordinates to float32."
    ),
    "regrid_native_precision_conversion": (
        "Regridding may analyze and synthesize at different supported native "
        "precisions, converting coefficients in the extension."
    ),
    "instance_representation": (
        "The C++ class has a richer diagnostic repr than original PySpharm."
    ),
    "safer_validation": (
        "Shape, dtype, count, truncation, finite-value, and grid configuration "
        "constraints receive explicit validation."
    ),
    "regular_high_degree_spherepack_omission": (
        "At regular-grid lmax=nlat-1, SPHEREPACK synthesis omits odd-m "
        "coefficients; DUCC evaluates the analytical basis."
    ),
    "getpsichi_zero_truncation": (
        "At ntrunc=0, the reference may raise IndexError while spharm_cpp "
        "returns zero streamfunction and velocity-potential maps, following the "
        "constant-scalar differential convention."
    ),
}

UNAPPROVED_DIFFERENCES = {}

PUBLIC_CONTRACTS = {**LEGACY_CONTRACTS, **ADDITIVE_CONTRACTS}
