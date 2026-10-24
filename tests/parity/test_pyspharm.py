# SPDX-License-Identifier: GPL-2.0-or-later
"""Differential PySpharm/SPHEREPACK compatibility tests."""

from __future__ import annotations

import inspect
import math
from decimal import Decimal, localcontext

import numpy as np
import pytest

import spharm_cpp as cpp

from . import conftest as parity
from .manifest import APPROVED_DIFFERENCES, PUBLIC_CONTRACTS


def _pair(gridtype: str, nlat: int, nlon: int, dtype=np.float32):
    reference = parity.REFERENCE.Spharmt(
        nlon, nlat, gridtype=gridtype, legfunc="stored"
    )
    current = cpp.Spharmt(nlon, nlat, gridtype=gridtype, dtype=dtype)
    return reference, current


def _indices(ntrunc: int):
    orders, degrees = cpp.getspecindx(ntrunc)
    return np.asarray(orders), np.asarray(degrees)


def _spectra(ntrunc: int, nfields: int, *, scale: float, seed: int):
    rng = np.random.default_rng(seed)
    nalm = (ntrunc + 1) * (ntrunc + 2) // 2
    values = (
        rng.standard_normal((nalm, nfields)) + 1j * rng.standard_normal((nalm, nfields))
    ).astype(np.complex64)
    m, _ = _indices(ntrunc)
    values[m == 0] = values[m == 0].real
    values *= np.float32(scale)
    return values


def _as_grid_shape(values):
    return np.asarray(values).reshape(np.asarray(values).shape[:2] + (-1,))


@pytest.mark.contract("Spharmt.__init__")
@pytest.mark.contract("Spharmt.attributes")
@pytest.mark.approved_difference("default_precision")
@pytest.mark.approved_difference("legendre_storage")
@pytest.mark.approved_difference("backend")
def test_constructor_defaults_signature_attributes_and_storage_policy():
    reference = parity.REFERENCE.Spharmt(18, 9, gridtype="gaussian")
    current = cpp.Spharmt(18, 9, gridtype="gaussian")

    signature = inspect.signature(parity.REFERENCE.Spharmt)
    assert tuple(signature.parameters) == (
        "nlon",
        "nlat",
        "rsphere",
        "gridtype",
        "legfunc",
    )
    assert signature.parameters["rsphere"].default == pytest.approx(6.3712e6)
    assert signature.parameters["gridtype"].default == "regular"
    assert signature.parameters["legfunc"].default == "stored"
    assert (
        "__init__(self, nlon: int, nlat: int, rsphere: float = 6371200.0, "
        "gridtype: str = 'regular', legfunc: str = 'computed', *, dtype: object | None = None)"
    ) in cpp.Spharmt.__init__.__doc__

    assert (current.nlon, current.nlat, current.rsphere, current.gridtype) == (
        reference.nlon,
        reference.nlat,
        reference.rsphere,
        reference.gridtype,
    )
    assert current.legfunc == "computed"
    assert current.dtype == "float64"
    current.user_note = "mutable Python instance state"
    assert current.user_note == "mutable Python instance state"
    for name, value in (
        ("nlon", 20),
        ("nlat", 10),
        ("rsphere", 1.0),
        ("gridtype", "regular"),
        ("legfunc", "stored"),
        ("dtype", "float32"),
    ):
        with pytest.raises((AttributeError, TypeError)):
            setattr(current, name, value)
        with pytest.raises((AttributeError, TypeError)):
            delattr(current, name)
    with pytest.raises(ValueError, match="only legfunc='computed'"):
        cpp.Spharmt(18, 9, legfunc="stored")
    del current.user_note
    assert not hasattr(current, "user_note")
    assert reference.legfunc == "stored"


@pytest.mark.contract("Spharmt.__repr__")
@pytest.mark.approved_difference("instance_representation")
def test_representation_is_diagnostic_and_class_names_remain_distinct():
    reference = parity.REFERENCE.Spharmt(18, 9, legfunc="stored")
    current = cpp.Spharmt(18, 9)
    assert type(current).__name__ == type(reference).__name__ == "Spharmt"
    assert "Spharmt(18, 9" in repr(current)
    assert "dtype=float64" in repr(current)
    # pyspharm-syl 1.2.1 has no Python __repr__; upstream commit a464ebf does.
    assert "Spharmt" in repr(reference)
    assert cpp is not parity.REFERENCE
    assert cpp.__name__ == "spharm_cpp"


@pytest.mark.contract("Spharmt.grdtospec")
@pytest.mark.contract("Spharmt.spectogrd")
@pytest.mark.approved_difference("regular_analysis_limit")
@pytest.mark.approved_difference("shape_at_zero_truncation")
@pytest.mark.parametrize(
    ("gridtype", "nlat", "nlon", "ntrunc"),
    [
        ("regular", 9, 17, 7),
        ("regular", 10, 18, 8),
        ("gaussian", 9, 16, 8),
        ("gaussian", 10, 17, 8),
    ],
)
def test_scalar_analysis_synthesis_parity_across_grid_parity(
    gridtype, nlat, nlon, ntrunc
):
    reference, current = _pair(gridtype, nlat, nlon)
    nalm = (ntrunc + 1) * (ntrunc + 2) // 2
    spec = _spectra(ntrunc, 2, scale=1.0, seed=nlat * 31 + nlon)

    expected_grid = reference.spectogrd(spec)
    actual_grid = current.spectogrd(spec)
    assert actual_grid.shape == expected_grid.shape == (nlat, nlon, 2)
    assert actual_grid.dtype == expected_grid.dtype == np.float32
    parity.assert_close(
        "spectogrd_float32", "float32", actual_grid, expected_grid, rtol=5e-6, atol=2e-6
    )

    expected_spec = reference.grdtospec(expected_grid, ntrunc)
    actual_spec = current.grdtospec(expected_grid, ntrunc)
    assert actual_spec.shape == expected_spec.shape == (nalm, 2)
    assert actual_spec.dtype == expected_spec.dtype == np.complex64
    parity.assert_close(
        "grdtospec_float32", "float32", actual_spec, expected_spec, rtol=3e-6, atol=3e-6
    )

    one_grid = expected_grid[:, :, 0]
    assert current.spectogrd(spec[:, 0]).shape == reference.spectogrd(spec[:, 0]).shape
    assert current.grdtospec(one_grid, ntrunc).ndim == 1
    assert reference.grdtospec(one_grid, ntrunc).ndim == 1

    # Rank-three singleton batches retain their field axis in both APIs.
    singleton = expected_grid[:, :, :1]
    assert current.grdtospec(singleton, ntrunc).shape == (nalm, 1)
    assert reference.grdtospec(singleton, ntrunc).shape == (nalm, 1)

    zero_reference = reference.grdtospec(one_grid, 0)
    zero_current = current.grdtospec(one_grid, 0)
    assert np.asarray(zero_reference).size == 1
    assert zero_current.shape == (1,)


@pytest.mark.contract("Spharmt.spectogrd")
def test_float64_synthesis_matches_independent_decimal_basis_samples():
    modes = ((0, 0, 0.75 + 0j), (1, 3, 0.5 + 0.25j), (2, 3, -0.4 + 0.3j))
    for gridtype, nlat, nlon in (
        ("regular", 8, 18),
        ("regular", 9, 17),
        ("gaussian", 8, 18),
        ("gaussian", 9, 17),
    ):
        transform = cpp.Spharmt(nlon, nlat, gridtype=gridtype)
        assert transform.dtype == "float64"
        latitudes = (
            np.linspace(90.0, -90.0, nlat)
            if gridtype == "regular"
            else cpp.gaussian_lats_wts(nlat)[0]
        )
        longitude = np.deg2rad(np.arange(nlon, dtype=np.float64) * 360.0 / nlon)
        ntrunc = 3
        index_m, index_l = _indices(ntrunc)
        for order, degree, coefficient in modes:
            index = int(np.flatnonzero((index_m == order) & (index_l == degree))[0])
            spectrum = np.zeros(index_m.size, dtype=np.complex128)
            spectrum[index] = coefficient
            expected = np.empty((nlat, nlon), dtype=np.float64)
            for row, latitude in enumerate(latitudes):
                pbar = _decimal_pbar(float(latitude), order, degree)
                if order == 0:
                    expected[row] = pbar * coefficient.real
                else:
                    expected[row] = (
                        2
                        * pbar
                        * (
                            coefficient.real * np.cos(order * longitude)
                            - coefficient.imag * np.sin(order * longitude)
                        )
                    )
            actual = transform.spectogrd(spectrum)
            assert actual.dtype == np.float64
            parity.assert_close(
                f"spectogrd_float64_analytic_{gridtype}_n{nlat}_m{order}_l{degree}",
                "float64",
                actual,
                expected,
                rtol=5e-14,
                atol=5e-14,
            )


@pytest.mark.contract("Spharmt.grdtospec")
def test_grid_analysis_handles_readonly_and_strided_maps_without_mutation():
    nlat, nlon, ntrunc = 9, 18, 6
    reference, current = _pair("gaussian", nlat, nlon)
    rng = np.random.default_rng(604)
    base = rng.standard_normal((nlat, nlon * 2)).astype(np.float32)
    strided = base[:, ::2]
    strided.flags.writeable = False
    before = strided.copy()
    expected = reference.grdtospec(strided, ntrunc)
    actual = current.grdtospec(strided, ntrunc)
    parity.assert_close(
        "grdtospec_readonly_strided_float32",
        "float32",
        actual,
        expected,
        rtol=5e-6,
        atol=2e-6,
    )
    np.testing.assert_array_equal(strided, before)


@pytest.mark.contract("Spharmt.spectogrd")
def test_spectral_synthesis_handles_readonly_strided_coefficients_without_mutation():
    nlat, nlon, ntrunc = 10, 18, 6
    reference, current = _pair("gaussian", nlat, nlon)
    rng = np.random.default_rng(338)
    storage = (
        rng.standard_normal(2 * ((ntrunc + 1) * (ntrunc + 2) // 2))
        + 1j * rng.standard_normal(2 * ((ntrunc + 1) * (ntrunc + 2) // 2))
    ).astype(np.complex64)
    strided = storage[::2]
    m, _ = _indices(ntrunc)
    strided[m == 0] = strided[m == 0].real
    strided.flags.writeable = False
    before = strided.copy()
    expected = reference.spectogrd(strided)
    actual = current.spectogrd(strided)
    parity.assert_close(
        "spectogrd_readonly_strided_float32",
        "float32",
        actual,
        expected,
        rtol=5e-6,
        atol=4e-6,
    )
    np.testing.assert_array_equal(strided, before)


@pytest.mark.contract("Spharmt.spectogrd")
@pytest.mark.contract("Spharmt.grdtospec")
def test_transform_outputs_own_storage_and_do_not_alias_inputs_or_fields():
    nlat, nlon, ntrunc, nfields = 9, 18, 5, 2
    transform = cpp.Spharmt(nlon, nlat, gridtype="gaussian", dtype=np.float64)
    nalm = (ntrunc + 1) * (ntrunc + 2) // 2
    spectra = np.zeros((nalm, nfields), dtype=np.complex128)
    spectra[0] = (1.25 + 0j, -0.75 + 0j)
    spectrum_before = spectra.copy()
    maps = transform.spectogrd(spectra)
    assert maps.shape == (nlat, nlon, nfields)
    assert not np.shares_memory(maps, spectra)
    np.testing.assert_array_equal(spectra, spectrum_before)
    second_field_before = maps[:, :, 1].copy()
    maps[:, :, 0] = 12.0
    np.testing.assert_array_equal(maps[:, :, 1], second_field_before)

    map_input = np.stack((np.ones((nlat, nlon)), np.full((nlat, nlon), 2.0)), axis=2)
    map_before = map_input.copy()
    coefficients = transform.grdtospec(map_input, ntrunc)
    assert coefficients.shape == (nalm, nfields)
    assert not np.shares_memory(coefficients, map_input)
    np.testing.assert_array_equal(map_input, map_before)


@pytest.mark.contract("Spharmt.getvrtdivspec")
@pytest.mark.contract("Spharmt.getuv")
@pytest.mark.contract("Spharmt.getgrad")
@pytest.mark.contract("Spharmt.getpsichi")
def test_vector_operations_match_independent_spherepack():
    nlat, nlon, ntrunc = 12, 24, 6
    reference, current = _pair("gaussian", nlat, nlon)
    vrt = _spectra(ntrunc, 2, scale=1e-6, seed=781)
    div = _spectra(ntrunc, 2, scale=1e-6, seed=782)
    vrt[0] = 0
    div[0] = 0

    expected_winds = reference.getuv(vrt, div)
    actual_winds = current.getuv(vrt, div)
    assert actual_winds[0].shape == expected_winds[0].shape == (nlat, nlon, 2)
    for component, (actual, expected) in enumerate(zip(actual_winds, expected_winds)):
        parity.assert_close(
            f"getuv_component_{component}_float32",
            "float32",
            actual,
            expected,
            rtol=3e-6,
            atol=1e-5,
        )

    expected_components = reference.getvrtdivspec(*expected_winds, ntrunc)
    actual_components = current.getvrtdivspec(*expected_winds, ntrunc)
    for component, (actual, expected) in enumerate(
        zip(actual_components, expected_components)
    ):
        parity.assert_close(
            f"getvrtdivspec_component_{component}_float32",
            "float32",
            actual,
            expected,
            rtol=4e-6,
            atol=2e-12,
        )

    potential = _spectra(ntrunc, 2, scale=1e6, seed=783)
    expected_gradient = reference.getgrad(potential)
    actual_gradient = current.getgrad(potential)
    for component, (actual, expected) in enumerate(
        zip(actual_gradient, expected_gradient)
    ):
        parity.assert_close(
            f"getgrad_component_{component}_float32",
            "float32",
            actual,
            expected,
            rtol=3e-6,
            atol=4e-6,
        )

    potential_scale = np.float32(1e-9)
    small_winds = reference.getuv(vrt * potential_scale, div * potential_scale)
    expected_potentials = reference.getpsichi(*small_winds, ntrunc)
    actual_potentials = current.getpsichi(*small_winds, ntrunc)
    for component, (actual, expected) in enumerate(
        zip(actual_potentials, expected_potentials)
    ):
        parity.assert_close(
            f"getpsichi_component_{component}_float32",
            "float32",
            actual,
            expected,
            rtol=3e-5,
            atol=2e-2,
        )


@pytest.mark.contract("Spharmt.getvrtdivspec")
@pytest.mark.contract("Spharmt.getuv")
@pytest.mark.contract("Spharmt.getgrad")
@pytest.mark.contract("Spharmt.getpsichi")
@pytest.mark.parametrize("gridtype", ["regular", "gaussian"])
def test_spin_one_operations_at_and_below_even_longitude_nyquist(gridtype):
    nlat, nlon = 10, 16
    reference, current = _pair(gridtype, nlat, nlon)
    nyquist = nlon // 2
    reference_boundary = {}
    current_boundary = {}

    for order in (nyquist - 1, nyquist):
        nalm = (order + 1) * (order + 2) // 2
        orders, degrees = _indices(order)
        index = int(np.flatnonzero((orders == order) & (degrees == order))[0])
        for label, coefficient in (("real", 1.0), ("imaginary", 1.0j)):
            potential = np.zeros(nalm, dtype=np.complex64)
            potential[index] = coefficient
            reference_gradient = reference.getgrad(potential)
            current_gradient = current.getgrad(potential)
            if order < nyquist:
                for component, (actual, expected) in enumerate(
                    zip(current_gradient, reference_gradient)
                ):
                    parity.assert_close(
                        f"getgrad_near_nyquist_{gridtype}_m{order}_{label}_component_{component}",
                        "float32",
                        actual,
                        expected,
                        rtol=4e-6,
                        atol=2e-12,
                    )
            else:
                assert max(np.max(np.abs(value)) for value in reference_gradient) == 0
                assert max(np.max(np.abs(value)) for value in current_gradient) > 0
                reference_boundary[f"getgrad_{label}_max_abs"] = 0.0
                current_boundary[f"getgrad_{label}_max_abs"] = float(
                    max(np.max(np.abs(value)) for value in current_gradient)
                )

            zero = np.zeros_like(potential)
            for component_name, vrt, div in (
                ("divergence", zero, potential),
                ("vorticity", potential, zero),
            ):
                reference_winds = reference.getuv(vrt, div)
                current_winds = current.getuv(vrt, div)
                if order < nyquist:
                    for component, (actual, expected) in enumerate(
                        zip(current_winds, reference_winds)
                    ):
                        parity.assert_close(
                            f"getuv_near_nyquist_{gridtype}_m{order}_{label}_{component_name}_component_{component}",
                            "float32",
                            actual,
                            expected,
                            # SPHEREPACK's regular-grid single-precision wind
                            # synthesis accumulates about 1e-4 relative error
                            # for this highest tested order.
                            rtol=1e-4,
                            atol=0.1,
                        )
                    reference_analysis = reference.getvrtdivspec(
                        *reference_winds, order
                    )
                    current_analysis = current.getvrtdivspec(*reference_winds, order)
                    for component, (actual, expected) in enumerate(
                        zip(current_analysis, reference_analysis)
                    ):
                        parity.assert_close(
                            f"getvrtdivspec_near_nyquist_{gridtype}_m{order}_{label}_{component_name}_component_{component}",
                            "float32",
                            actual,
                            expected,
                            rtol=5e-6,
                            atol=3e-6,
                        )
                else:
                    assert max(np.max(np.abs(value)) for value in reference_winds) == 0
                    assert max(np.max(np.abs(value)) for value in current_winds) > 0
                    ref_vrt, ref_div = reference.getvrtdivspec(
                        *reference_winds, nyquist
                    )
                    assert np.max(np.abs(ref_vrt)) == 0
                    assert np.max(np.abs(ref_div)) == 0
                    current_boundary[f"getuv_{label}_{component_name}_max_abs"] = float(
                        max(np.max(np.abs(value)) for value in current_winds)
                    )
                    reference_boundary[f"getuv_{label}_{component_name}_max_abs"] = 0.0
                    with pytest.raises(ValueError, match=r"floor\(\(nlon-1\)/2\)"):
                        current.getvrtdivspec(*current_winds, nyquist)
                    with pytest.raises(ValueError, match=r"floor\(\(nlon-1\)/2\)"):
                        current.getpsichi(*current_winds, nyquist)
                    default_analysis = current.getvrtdivspec(*current_winds)
                    assert (
                        default_analysis[0].shape
                        == default_analysis[1].shape
                        == ((nyquist * (nyquist + 1) // 2),)
                    )

    parity.record_investigation(
        f"spin_one_even_longitude_nyquist_{gridtype}",
        {
            "status": "PASS",
            "gridtype": gridtype,
            "nlat": nlat,
            "nlon": nlon,
            "resolved_order": nyquist - 1,
            "nyquist_order": nyquist,
            "ducc_spin1_analysis_limit": (nlon - 1) // 2,
            "reference_exact_nyquist_getgrad_and_getuv_max_abs": reference_boundary,
            "cpp_exact_nyquist_getgrad_and_getuv_max_abs": current_boundary,
            "reference_exact_nyquist_getvrtdivspec_max_abs": 0.0,
            "cpp_exact_nyquist_analysis": "explicit requests raise ValueError",
            "interpretation": (
                "The sampled SPHEREPACK spin-1 path omits the exact longitude "
                "Nyquist order. DUCC synthesis produces observable wind samples "
                "for the tested real and imaginary coefficient components, while "
                "DUCC analysis documents a maximum order of floor((nlon-1)/2). "
                "Orders immediately below Nyquist match the independent reference."
            ),
        },
    )


@pytest.mark.contract("Spharmt.specsmooth")
def test_spectral_smoothing_matches_on_an_explicit_band_limited_field():
    nlat, nlon, ntrunc = 12, 24, 6
    reference, current = _pair("gaussian", nlat, nlon)
    spec = _spectra(ntrunc, 1, scale=1e5, seed=109)
    field = reference.spectogrd(spec)
    smooth = np.linspace(0.35, 1.0, nlat, dtype=np.float32)
    smooth[ntrunc + 1 :] = 0
    expected = reference.specsmooth(field, smooth)
    actual = current.specsmooth(field, smooth)
    parity.assert_close(
        "specsmooth_float32", "float32", actual, expected, rtol=4e-6, atol=0.3
    )


@pytest.mark.contract("regrid")
@pytest.mark.parametrize(
    ("source_type", "target_type"),
    [("regular", "gaussian"), ("gaussian", "regular")],
)
def test_regrid_matches_with_explicit_shared_band_and_optional_filter(
    source_type, target_type
):
    source_ref, source_cpp = _pair(source_type, 12, 24)
    target_ref, target_cpp = _pair(target_type, 10, 20)
    ntrunc = 6
    spectrum = _spectra(ntrunc, 2, scale=1e5, seed=501)
    field = source_ref.spectogrd(spectrum)
    factors = np.linspace(0.2, 1.0, target_ref.nlat, dtype=np.float32)
    expected = parity.REFERENCE.regrid(
        source_ref, target_ref, field, ntrunc=ntrunc, smooth=factors
    )
    actual = cpp.regrid(source_cpp, target_cpp, field, ntrunc=ntrunc, smooth=factors)
    assert actual.shape == expected.shape == (target_ref.nlat, target_ref.nlon, 2)
    parity.assert_close(
        f"regrid_{source_type}_to_{target_type}_float32",
        "float32",
        actual,
        expected,
        rtol=4e-6,
        atol=0.4,
    )


@pytest.mark.contract("regrid")
@pytest.mark.approved_difference("regrid_native_precision_conversion")
def test_regrid_converts_precision_between_grid_instances():
    source_reference = parity.REFERENCE.Spharmt(
        24, 12, gridtype="regular", legfunc="stored"
    )
    target_reference = parity.REFERENCE.Spharmt(
        20, 10, gridtype="gaussian", legfunc="stored"
    )
    source_cpp = cpp.Spharmt(24, 12, gridtype="regular", dtype=np.float64)
    target_cpp = cpp.Spharmt(20, 10, gridtype="gaussian", dtype=np.float32)
    ntrunc = 6
    spectrum = _spectra(ntrunc, 2, scale=0.75, seed=211).astype(np.complex128)
    field64 = source_cpp.spectogrd(spectrum)
    factors = np.linspace(0.4, 1.0, target_reference.nlat, dtype=np.float32)
    for smooth in (None, factors):
        expected = parity.REFERENCE.regrid(
            source_reference,
            target_reference,
            field64.astype(np.float32),
            ntrunc=ntrunc,
            smooth=smooth,
        )
        actual = cpp.regrid(
            source_cpp,
            target_cpp,
            field64,
            ntrunc=ntrunc,
            smooth=smooth,
        )
        assert actual.shape == expected.shape == (10, 20, 2)
        assert actual.dtype == np.float32
        parity.assert_close(
            f"regrid_float64_to_float32_smooth_{smooth is not None}",
            "float32 output",
            actual,
            expected,
            rtol=6e-6,
            atol=5e-6,
        )


@pytest.mark.contract("gaussian_lats_wts")
def test_gaussian_latitudes_and_weights_match_spherepack():
    for nlat in (3, 4, 5, 8, 9, 16, 31):
        expected_lat, expected_weight = parity.REFERENCE.gaussian_lats_wts(nlat)
        actual_lat, actual_weight = cpp.gaussian_lats_wts(nlat)
        assert actual_lat.dtype == actual_weight.dtype == np.float64
        assert actual_lat.shape == actual_weight.shape == (nlat,)
        parity.assert_close(
            "gaussian_lats_wts_latitude_float64",
            "float64",
            actual_lat,
            expected_lat,
            rtol=0,
            atol=2e-12,
        )
        parity.assert_close(
            "gaussian_lats_wts_weight_float64",
            "float64",
            actual_weight,
            expected_weight,
            rtol=0,
            atol=2e-14,
        )


@pytest.mark.contract("getspecindx")
@pytest.mark.approved_difference("shape_at_zero_truncation")
def test_spectral_indices_match_numpy_native_reference_order():
    for ntrunc in (0, 1, 3, 16, 64):
        expected_m, expected_l = parity.REFERENCE.getspecindx(ntrunc)
        actual_m, actual_l = cpp.getspecindx(ntrunc)
        np.testing.assert_array_equal(actual_m, np.asarray(expected_m))
        np.testing.assert_array_equal(actual_l, np.asarray(expected_l))
        assert actual_m.dtype == actual_l.dtype == np.dtype(np.intp)
        assert actual_m.ndim == actual_l.ndim == 1
    expected_m0, expected_l0 = parity.REFERENCE.getspecindx(0)
    actual_m0, actual_l0 = cpp.getspecindx(0)
    assert np.asarray(expected_m0).item() == actual_m0[0] == 0
    assert np.asarray(expected_l0).item() == actual_l0[0] == 0


@pytest.mark.contract("getgeodesicpts")
@pytest.mark.approved_difference("geodesic_float64_default")
def test_geodesic_float32_order_and_double_precision_geometry():
    for m in (1, 2, 3, 4, 6, 10, 20):
        expected_lat, expected_lon = parity.REFERENCE.getgeodesicpts(m)
        actual_lat, actual_lon = cpp.getgeodesicpts(m, dtype=np.float32)
        assert actual_lat.shape == actual_lon.shape == (10 * (m - 1) ** 2 + 2,)
        assert actual_lat.dtype == actual_lon.dtype == np.float32
        assert expected_lat.shape == expected_lon.shape == actual_lat.shape
        longitude_delta = (
            actual_lon.astype(np.float64) - expected_lon + 180
        ) % 360 - 180
        parity.assert_close(
            "getgeodesicpts_latitude_float32",
            "float32",
            actual_lat,
            expected_lat,
            rtol=0,
            atol=7e-5,
        )
        parity.assert_close(
            "getgeodesicpts_longitude_float32",
            "float32",
            longitude_delta,
            np.zeros_like(longitude_delta),
            rtol=0,
            atol=7e-5,
        )

    for m in (2, 3, 10, 20):
        lat64, lon64 = cpp.getgeodesicpts(m)
        assert lat64.shape == lon64.shape == (10 * (m - 1) ** 2 + 2,)
        assert lat64.dtype == lon64.dtype == np.float64
        assert (lat64[0], lat64[1], lon64[0], lon64[1]) == (90.0, -90.0, 0.0, 0.0)
        phi64 = np.deg2rad(lat64)
        lam64 = np.deg2rad(lon64)
        norms = np.sqrt(
            (np.cos(phi64) * np.cos(lam64)) ** 2
            + (np.cos(phi64) * np.sin(lam64)) ** 2
            + np.sin(phi64) ** 2
        )
        np.testing.assert_allclose(norms, 1.0, rtol=0, atol=4e-16)

    lat64, lon64 = cpp.getgeodesicpts(2)
    expected_vertex_latitude = math.degrees(math.asin(1 / math.sqrt(5)))
    np.testing.assert_allclose(
        np.abs(lat64[2:]), expected_vertex_latitude, rtol=0, atol=5e-13
    )
    phi = np.deg2rad(lat64)
    lam = np.deg2rad(lon64)
    xyz_norm = np.sqrt(
        (np.cos(phi) * np.cos(lam)) ** 2
        + (np.cos(phi) * np.sin(lam)) ** 2
        + np.sin(phi) ** 2
    )
    np.testing.assert_allclose(xyz_norm, 1.0, rtol=0, atol=3e-16)
    for m in (1, 2, 3, 4, 6, 10, 20):
        lat64, lon64 = cpp.getgeodesicpts(m)
        lat32, lon32 = cpp.getgeodesicpts(m, dtype=np.float32)
        np.testing.assert_array_equal(lat32, lat64.astype(np.float32))
        np.testing.assert_array_equal(lon32, lon64.astype(np.float32))
    with pytest.raises(TypeError, match="float32 or numpy.float64"):
        cpp.getgeodesicpts(2, dtype=np.int32)
    with pytest.raises(TypeError):
        cpp.getgeodesicpts(2, np.float32)


def _decimal_pbar(latitude: float, m: int, ell: int) -> float:
    """Independent high-precision normalized associated-Legendre recurrence."""
    phi = math.radians(latitude)
    x = Decimal.from_float(math.sin(phi))
    cos_phi = Decimal.from_float(math.cos(phi))
    with localcontext() as context:
        context.prec = 90
        diagonal = Decimal(1) / Decimal(2).sqrt()
        for order in range(1, m + 1):
            diagonal *= cos_phi * (Decimal(2 * order + 1) / Decimal(2 * order)).sqrt()
        if ell == m:
            return float(diagonal)
        previous = diagonal
        current = x * Decimal(2 * m + 3).sqrt() * diagonal
        if ell == m + 1:
            return float(current)
        for degree in range(m + 2, ell + 1):
            ld, md = Decimal(degree), Decimal(m)
            lp, lm = ld + md, ld - md
            a = (((2 * ld - 1) * (2 * ld + 1)) / (lp * lm)).sqrt()
            b = (((2 * ld + 1) * (lp - 1) * (lm - 1)) / ((2 * ld - 3) * lp * lm)).sqrt()
            previous, current = current, a * x * current - b * previous
        return float(current)


@pytest.mark.contract("legendre")
def test_legendre_float32_matches_reference_and_float64_matches_decimal():
    for ntrunc in (0, 1, 8, 32):
        for latitude in (-90.0, -33.25, 0.0, 41.5, 90.0):
            expected = parity.REFERENCE.legendre(latitude, ntrunc)
            actual = cpp.legendre(latitude, ntrunc, dtype=np.float32)
            parity.assert_close(
                "legendre_float32",
                "float32",
                actual,
                expected,
                rtol=5e-6,
                atol=1e-5,
            )
    ntrunc = 192
    m_index, l_index = _indices(ntrunc)
    latitude = 37.25
    assert cpp.legendre(latitude, 4).dtype == np.float64
    actual64 = cpp.legendre(latitude, ntrunc, dtype=np.float64)
    references = np.asarray(
        [
            _decimal_pbar(latitude, int(order), int(degree))
            for order, degree in zip(m_index, l_index)
            if (int(order), int(degree))
            in {(0, 0), (0, 128), (1, 191), (32, 160), (96, 192)}
        ],
        dtype=np.float64,
    )
    select = np.isin(m_index * 1000 + l_index, [0, 128, 1191, 32160, 96192])
    # Exact integer pairs avoid any ambiguity in the triangular coefficient index.
    expected64 = np.asarray(
        [
            _decimal_pbar(latitude, int(order), int(degree))
            for order, degree in zip(m_index[select], l_index[select])
        ]
    )
    parity.assert_close(
        "legendre_float64_decimal",
        "float64",
        actual64[select],
        expected64,
        rtol=3e-12,
        atol=3e-14,
    )
    assert references.shape == expected64.shape


@pytest.mark.contract("specintrp")
def test_scalar_spectral_interpolation_matches_reference_at_lon_and_lat_points():
    ntrunc = 12
    nalm = (ntrunc + 1) * (ntrunc + 2) // 2
    rng = np.random.default_rng(815)
    spec = (rng.standard_normal(nalm) + 1j * rng.standard_normal(nalm)).astype(
        np.complex64
    )
    m, _ = _indices(ntrunc)
    spec[m == 0] = spec[m == 0].real
    for latitude in (-90.0, -11.25, 0.0, 32.75, 90.0):
        pbar = parity.REFERENCE.legendre(latitude, ntrunc)
        for longitude in (-720.0, -17.5, 0.0, 23.75, 179.0, 360.0):
            expected = parity.REFERENCE.specintrp(longitude, spec, pbar)
            actual = cpp.specintrp(longitude, spec, pbar)
            assert isinstance(expected, float)
            assert isinstance(actual, float)
            parity.assert_close(
                "specintrp_float32",
                "float32",
                np.asarray(actual),
                np.asarray(expected),
                rtol=5e-6,
                atol=5e-6,
            )


@pytest.mark.contract("legendre_many")
def test_legendre_many_matches_repeated_reference_scalar_calls():
    latitudes = np.asarray([-90.0, -33.25, 0.0, 41.5, 90.0], dtype=np.float32)
    ntrunc = 32
    expected = np.stack(
        [parity.REFERENCE.legendre(float(latitude), ntrunc) for latitude in latitudes]
    )
    actual = cpp.legendre_many(latitudes, ntrunc, dtype=np.float32)
    parity.assert_close(
        "legendre_many_float32",
        "float32",
        actual,
        expected,
        rtol=5e-6,
        atol=1e-5,
    )


@pytest.mark.contract("specintrp_many")
def test_specintrp_many_matches_repeated_reference_scalar_calls():
    ntrunc = 8
    nalm = (ntrunc + 1) * (ntrunc + 2) // 2
    rng = np.random.default_rng(912)
    spec = (rng.standard_normal(nalm) + 1j * rng.standard_normal(nalm)).astype(
        np.complex64
    )
    m, _ = _indices(ntrunc)
    spec[m == 0] = spec[m == 0].real
    latitudes = [-78.0, -10.0, 0.0, 24.5, 81.0]
    longitudes = np.asarray([-361.0, -30.0, 0.0, 75.25, 413.0])
    legendre_rows = np.stack(
        [parity.REFERENCE.legendre(latitude, ntrunc) for latitude in latitudes]
    )
    expected = np.asarray(
        [
            parity.REFERENCE.specintrp(lon, spec, pbar)
            for lon, pbar in zip(longitudes, legendre_rows)
        ]
    )
    actual = cpp.specintrp_many(longitudes, spec, legendre_rows)
    assert actual.dtype == np.float64
    parity.assert_close(
        "specintrp_many_float32",
        "float32",
        actual,
        expected,
        rtol=6e-6,
        atol=6e-6,
    )


@pytest.mark.contract("cpu_info")
def test_cpu_info_matches_upstream_ducc_schema_and_types():
    info = cpp.cpu_info()
    assert info["architecture"] in {"x86-64", "aarch64"}
    assert isinstance(info["features"], list)
    assert all(isinstance(feature, str) for feature in info["features"])
    assert isinstance(info["multiarch"], bool)
    base_keys = {"architecture", "features", "multiarch"}
    if info["multiarch"]:
        assert set(info) == base_keys | {
            "compiled_profiles",
            "available_profiles",
            "configured_limit",
            "active_profile",
        }
        assert info["active_profile"] in info["compiled_profiles"]
        assert info["active_profile"] in info["available_profiles"]
    else:
        assert set(info) == base_keys


@pytest.mark.contract("Spharmt.dtype")
def test_dtype_property_selects_consistent_array_precisions():
    for dtype, real_type, complex_type in (
        (np.float32, np.float32, np.complex64),
        (np.float64, np.float64, np.complex128),
    ):
        transform = cpp.Spharmt(18, 9, gridtype="gaussian", dtype=dtype)
        assert transform.dtype == np.dtype(dtype).name
        grid = np.zeros((9, 18), dtype=real_type)
        spectrum = transform.grdtospec(grid, 4)
        assert spectrum.dtype == np.dtype(complex_type)
        assert transform.spectogrd(spectrum).dtype == np.dtype(real_type)


@pytest.mark.contract("Spharmt.spectogrd")
@pytest.mark.contract("Spharmt.grdtospec")
@pytest.mark.approved_difference("shape_at_zero_truncation")
def test_zero_truncation_rank_and_batch_axis_are_consistent():
    reference, current = _pair("gaussian", 9, 18)
    spectrum = np.asarray([1.0 + 0j], dtype=np.complex64)
    assert current.spectogrd(spectrum).shape == reference.spectogrd(spectrum).shape
    coefficients = current.grdtospec(np.ones((9, 18), dtype=np.float32), 0)
    assert coefficients.shape == (1,)
    assert (
        np.asarray(reference.grdtospec(np.ones((9, 18), dtype=np.float32), 0)).shape
        == ()
    )


@pytest.mark.contract("Spharmt.getpsichi")
@pytest.mark.approved_difference("getpsichi_zero_truncation")
def test_degree_zero_velocity_potentials_follow_constant_scalar_convention():
    nlat, nlon = 9, 18
    current = cpp.Spharmt(nlon, nlat, gridtype="gaussian", dtype=np.float64)
    ugrid = np.zeros((nlat, nlon), dtype=np.float64)
    vgrid = np.zeros_like(ugrid)
    psi, chi = current.getpsichi(ugrid, vgrid, ntrunc=0)
    assert psi.shape == chi.shape == (nlat, nlon)
    assert psi.dtype == chi.dtype == np.float64
    np.testing.assert_array_equal(psi, np.zeros_like(psi))
    np.testing.assert_array_equal(chi, np.zeros_like(chi))

    reference = parity.REFERENCE.Spharmt(
        nlon, nlat, gridtype="gaussian", legfunc="stored"
    )
    with pytest.raises(IndexError):
        reference.getpsichi(
            np.zeros((nlat, nlon), dtype=np.float32),
            np.zeros((nlat, nlon), dtype=np.float32),
            ntrunc=0,
        )


@pytest.mark.contract("Spharmt.__init__")
@pytest.mark.contract("Spharmt.grdtospec")
@pytest.mark.approved_difference("safer_validation")
def test_explicit_validation_errors_for_invalid_configuration_and_requests():
    with pytest.raises(ValueError, match="nlon must be at least 4"):
        cpp.Spharmt(3, 9)
    with pytest.raises(ValueError, match="nlat must be at least 3"):
        cpp.Spharmt(18, 2)
    with pytest.raises(ValueError, match="finite and positive"):
        cpp.Spharmt(18, 9, rsphere=float("inf"))
    with pytest.raises(ValueError, match="gridtype"):
        cpp.Spharmt(18, 9, gridtype="unknown")
    with pytest.raises(TypeError, match="dtype must be numpy.float32 or numpy.float64"):
        cpp.Spharmt(18, 9, dtype=np.int32)
    transform = cpp.Spharmt(18, 9)
    with pytest.raises(ValueError, match="triangular"):
        transform.spectogrd(np.zeros(5, np.complex128))
    with pytest.raises(ValueError, match="ntrunc must be between"):
        transform.grdtospec(np.zeros((9, 18)), -1)


def test_zz_every_public_item_has_a_marked_validation():
    request = pytest.Config  # Keep pytest from interpreting this as a fixture.
    del request
    # Collection completeness is enforced in pytest_collection_finish below.
    assert PUBLIC_CONTRACTS


def _associated_legendre_matrix(latitude: np.ndarray, order: int, lmax: int):
    if order > lmax:
        return np.empty((latitude.size, 0), dtype=np.float64)
    values = []
    for lat in latitude:
        phi = math.radians(float(lat))
        x, cosine = math.sin(phi), math.cos(phi)
        pmm = 1.0 / math.sqrt(2.0)
        for m in range(1, order + 1):
            pmm *= cosine * math.sqrt((2 * m + 1) / (2 * m))
        row = [pmm]
        if order < lmax:
            previous = pmm
            current = x * math.sqrt(2 * order + 3) * pmm
            row.append(current)
            for degree in range(order + 2, lmax + 1):
                ld, md = float(degree), float(order)
                lp, lm = ld + md, ld - md
                a = math.sqrt(((2 * ld - 1) * (2 * ld + 1)) / (lp * lm))
                b = math.sqrt(
                    ((2 * ld + 1) * (lp - 1) * (lm - 1)) / ((2 * ld - 3) * lp * lm)
                )
                previous, current = current, a * x * current - b * previous
                row.append(current)
        values.append(row)
    return np.asarray(values, dtype=np.float64)


@pytest.mark.contract("Spharmt.spectogrd")
@pytest.mark.approved_difference("regular_high_degree_spherepack_omission")
@pytest.mark.parametrize("nlat", [8, 9])
def test_regular_maximum_degree_rank_and_spherepack_nullspace(nlat):
    nlon, lmax = 2 * nlat + 2, nlat - 1
    latitudes = np.linspace(90.0, -90.0, nlat)
    interior = latitudes[1:-1]
    for order in range(lmax + 1):
        full = _associated_legendre_matrix(latitudes, order, lmax)
        interior_matrix = _associated_legendre_matrix(interior, order, lmax)
        rank_high = int(
            np.linalg.matrix_rank(full if order == 0 else interior_matrix, tol=1e-12)
        )
        lower_matrix = _associated_legendre_matrix(interior, order, lmax - 1)
        lower_full = _associated_legendre_matrix(latitudes, order, lmax - 1)
        rank_lower = int(
            np.linalg.matrix_rank(lower_full if order == 0 else lower_matrix, tol=1e-12)
        )
        if order == 0:
            assert rank_high == nlat
            assert rank_lower == nlat - 1
        elif order == 1:
            assert interior_matrix.shape == (nlat - 2, nlat - 1)
            assert rank_high == nlat - 2
            assert rank_lower == nlat - 2
            _, _, right = np.linalg.svd(interior_matrix, full_matrices=True)
            null_vector = right[-1]
            assert np.linalg.norm(interior_matrix @ null_vector) < 2e-12
            assert np.linalg.norm(full @ null_vector) < 2e-12
        else:
            assert rank_high == nlat - order
            assert rank_lower == nlat - order - 1 if order < lmax else rank_lower == 0

    reference, current = _pair("regular", nlat, nlon)
    current64 = cpp.Spharmt(nlon, nlat, gridtype="regular")

    # At L=nlat-1 the m=1 latitude matrix has one null direction. Construct it
    # explicitly and show that the analytical DUCC basis combination vanishes
    # at every sampled regular-grid latitude (the poles add exact zeros).
    _, _, right = np.linalg.svd(
        _associated_legendre_matrix(interior, 1, lmax), full_matrices=True
    )
    null_vector = right[-1]
    nalm = (lmax + 1) * (lmax + 2) // 2
    mstart = lmax + 1
    null_spectrum = np.zeros(nalm, dtype=np.complex128)
    null_spectrum[mstart : mstart + null_vector.size] = null_vector
    null_grid = current64.spectogrd(null_spectrum)
    assert np.max(np.abs(null_grid)) < 2e-12

    # SPHEREPACK accepts this truncation, but its regular SHSES path drops the
    # highest-degree coefficients for odd zonal wavenumbers. Record the direct
    # differential results instead of hiding them as round-trip parity.
    coefficient_m, coefficient_l = _indices(lmax)
    unsupported = []
    longitude = np.arange(nlon, dtype=np.float64) * (360.0 / nlon)
    for order in range(lmax + 1):
        index = int(
            np.flatnonzero((coefficient_m == order) & (coefficient_l == lmax))[0]
        )
        spectrum32 = np.zeros(nalm, dtype=np.complex64)
        spectrum32[index] = 1
        spectrum64 = spectrum32.astype(np.complex128)
        reference_grid = reference.spectogrd(spectrum32)
        current_grid = current64.spectogrd(spectrum64)
        pbar = np.asarray(
            [_decimal_pbar(float(latitude), order, lmax) for latitude in latitudes]
        )
        if order == 0:
            analytical = np.broadcast_to(pbar[:, None], (nlat, nlon))
        else:
            analytical = (
                2.0 * pbar[:, None] * np.cos(np.deg2rad(order * longitude))[None, :]
            )
        parity.assert_close(
            f"regular_lmax_analytic_mode_m{order}_float64",
            "float64",
            current_grid,
            analytical,
            rtol=3e-13,
            atol=3e-13,
        )
        expected_legacy = analytical.astype(np.float32)
        if order % 2 and order > 0:
            expected_legacy = np.zeros_like(expected_legacy)
            unsupported.append(order)
        parity.assert_close(
            f"spherepack_regular_lmax_mode_m{order}_float32",
            "float32",
            reference_grid,
            expected_legacy,
            rtol=4e-6,
            atol=3e-6,
        )
        agrees = parity.record_comparison(
            f"cpp_vs_spherepack_regular_lmax_mode_m{order}",
            "float64_vs_float32",
            current_grid,
            reference_grid,
            rtol=4e-6,
            atol=3e-6,
            approved_difference=(
                "regular_high_degree_spherepack_omission"
                if order > 0 and order % 2 == 1
                else None
            ),
        )
        if agrees != (order == 0 or order % 2 == 0):
            raise AssertionError(
                f"unexpected PySpharm/DUCC regular maximum-degree behavior at m={order}"
            )
        recovered = reference.grdtospec(reference_grid, lmax)[index]
        expected_recovered = 0.0 if order > 0 and order % 2 else 1.0
        np.testing.assert_allclose(recovered, expected_recovered, rtol=0, atol=4e-6)

    rejected = np.zeros((nlat, nlon), dtype=np.float32)
    with pytest.raises(ValueError, match="nlat-2"):
        current.grdtospec(rejected, lmax)
    lower_spec = reference.grdtospec(rejected, lmax - 1)
    assert current.grdtospec(rejected, lmax - 1).shape == lower_spec.shape
    assert unsupported == list(range(1, lmax + 1, 2))


@pytest.mark.contract("Spharmt.grdtospec")
@pytest.mark.approved_difference("regular_analysis_limit")
@pytest.mark.parametrize("nlat", [8, 9])
def test_regular_analysis_limit_is_explicit_and_lower_band_remains_supported(nlat):
    nlon = 2 * nlat + 2
    maximum_supported = nlat - 2
    reference, current = _pair("regular", nlat, nlon)
    grid = np.zeros((nlat, nlon), dtype=np.float32)

    # SPHEREPACK accepts nlat-1; the public C++ API intentionally rejects it
    # because its supported transform backend requires nlat-2 for analysis.
    reference_result = reference.grdtospec(grid, nlat - 1)
    assert reference_result.shape == ((nlat) * (nlat + 1) // 2,)
    with pytest.raises(ValueError, match="nlat-2"):
        current.grdtospec(grid, nlat - 1)

    expected = reference.grdtospec(grid, maximum_supported)
    actual = current.grdtospec(grid, maximum_supported)
    assert actual.shape == expected.shape
    np.testing.assert_allclose(actual, expected, rtol=0, atol=0)


@pytest.mark.contract("Spharmt.getvrtdivspec")
@pytest.mark.contract("Spharmt.getpsichi")
@pytest.mark.contract("Spharmt.specsmooth")
@pytest.mark.contract("Spharmt.getuv")
@pytest.mark.contract("Spharmt.getgrad")
@pytest.mark.contract("regrid")
@pytest.mark.approved_difference("regular_analysis_limit")
def test_regular_analysis_limit_propagates_to_derived_operations():
    nlat, nlon, ntrunc = 8, 18, 6
    reference, current = _pair("regular", nlat, nlon)
    target_reference, target_current = _pair("gaussian", nlat + 1, nlon)
    vorticity = _spectra(ntrunc, 1, scale=1e-6, seed=603)[:, 0]
    divergence = _spectra(ntrunc, 1, scale=1e-6, seed=604)[:, 0]
    vorticity[0] = 0
    divergence[0] = 0
    expected_winds = reference.getuv(vorticity, divergence)
    actual_winds = current.getuv(vorticity, divergence)
    for component, (actual, expected) in enumerate(zip(actual_winds, expected_winds)):
        parity.assert_close(
            f"regular_getuv_component_{component}_float32",
            "float32",
            actual,
            expected,
            rtol=4e-6,
            atol=1e-5,
        )
    ugrid, vgrid = expected_winds

    actual_vrt, actual_div = current.getvrtdivspec(ugrid, vgrid)
    expected_vrt, expected_div = reference.getvrtdivspec(ugrid, vgrid, ntrunc=ntrunc)
    parity.assert_close(
        "regular_getvrtdivspec_default_float32",
        "float32",
        actual_vrt,
        expected_vrt,
        rtol=4e-6,
        atol=2e-12,
    )
    parity.assert_close(
        "regular_getvrtdivspec_divergence_default_float32",
        "float32",
        actual_div,
        expected_div,
        rtol=4e-6,
        atol=2e-12,
    )
    for operation in (current.getvrtdivspec, current.getpsichi):
        with pytest.raises(ValueError, match="nlat-2"):
            operation(ugrid, vgrid, ntrunc + 1)

    small_winds = (ugrid * np.float32(1e-3), vgrid * np.float32(1e-3))
    actual_psi, actual_chi = current.getpsichi(*small_winds)
    expected_psi, expected_chi = reference.getpsichi(*small_winds, ntrunc=ntrunc)
    for operation, actual, expected in zip(
        ("streamfunction", "velocity_potential"),
        (actual_psi, actual_chi),
        (expected_psi, expected_chi),
    ):
        parity.assert_close(
            f"regular_getpsichi_{operation}_default_float32",
            "float32",
            actual,
            expected,
            rtol=3e-5,
            atol=2e-2,
        )

    scalar_spectrum = _spectra(ntrunc, 1, scale=2e5, seed=605)[:, 0]
    actual_gradient = current.getgrad(scalar_spectrum)
    expected_gradient = reference.getgrad(scalar_spectrum)
    for component, (actual, expected) in enumerate(
        zip(actual_gradient, expected_gradient)
    ):
        parity.assert_close(
            f"regular_getgrad_component_{component}_float32",
            "float32",
            actual,
            expected,
            rtol=4e-6,
            atol=4e-6,
        )

    spectrum = _spectra(ntrunc, 2, scale=4e4, seed=602)
    band_limited = reference.spectogrd(spectrum)
    smoothing = np.linspace(0.3, 1.0, nlat, dtype=np.float32)
    expected_smoothed = reference.specsmooth(band_limited, smoothing)
    actual_smoothed = current.specsmooth(band_limited, smoothing)
    parity.assert_close(
        "regular_specsmooth_default_shared_band_float32",
        "float32",
        actual_smoothed,
        expected_smoothed,
        rtol=5e-6,
        atol=0.2,
    )

    expected_regridded = parity.REFERENCE.regrid(
        reference, target_reference, band_limited
    )
    actual_regridded = cpp.regrid(current, target_current, band_limited)
    parity.assert_close(
        "regular_regrid_default_shared_band_float32",
        "float32",
        actual_regridded,
        expected_regridded,
        rtol=5e-6,
        atol=0.2,
    )
    with pytest.raises(ValueError, match="nlat-2"):
        cpp.regrid(current, target_current, band_limited, ntrunc=ntrunc + 1)


@pytest.mark.contract("Spharmt.grdtospec")
@pytest.mark.contract("Spharmt.spectogrd")
def test_gaussian_maximum_degree_nlat_minus_one_is_full_rank():
    nlat, nlon, lmax = 8, 18, 7
    latitudes, _ = cpp.gaussian_lats_wts(nlat)
    for order in range(lmax + 1):
        matrix = _associated_legendre_matrix(latitudes, order, lmax)
        rank = int(np.linalg.matrix_rank(matrix, tol=1e-12))
        assert rank == lmax - order + 1
    reference, current = _pair("gaussian", nlat, nlon)
    coefficient_m, coefficient_l = _indices(lmax)
    for order in range(lmax + 1):
        index = int(
            np.flatnonzero((coefficient_m == order) & (coefficient_l == lmax))[0]
        )
        mode = np.zeros((lmax + 1) * (lmax + 2) // 2, dtype=np.complex64)
        mode[index] = 1
        parity.assert_close(
            f"gaussian_lmax_basis_m{order}_float32",
            "float32",
            current.spectogrd(mode),
            reference.spectogrd(mode),
            rtol=4e-6,
            atol=4e-6,
        )
    spectrum = _spectra(lmax, 1, scale=15.0, seed=992)[:, 0]
    expected_grid = reference.spectogrd(spectrum)
    actual_grid = current.spectogrd(spectrum)
    parity.assert_close(
        "gaussian_maxdegree_synthesis_float32",
        "float32",
        actual_grid,
        expected_grid,
        rtol=3e-6,
        atol=3e-5,
    )
    expected_coefficients = reference.grdtospec(expected_grid, lmax)
    actual_coefficients = current.grdtospec(expected_grid, lmax)
    parity.assert_close(
        "gaussian_maxdegree_analysis_float32",
        "float32",
        actual_coefficients,
        expected_coefficients,
        rtol=4e-6,
        atol=3e-6,
    )
    parity.assert_close(
        "gaussian_maxdegree_recovery_float32",
        "float32",
        actual_coefficients,
        spectrum,
        rtol=4e-6,
        atol=3e-6,
    )


@pytest.mark.contract("Spharmt.grdtospec")
@pytest.mark.contract("Spharmt.spectogrd")
@pytest.mark.parametrize("gridtype", ["regular", "gaussian"])
def test_even_longitude_scalar_nyquist_matches_spherepack_and_analytic_samples(
    gridtype,
):
    nlat, nlon, lmax = 10, 16, 8
    reference, current = _pair(gridtype, nlat, nlon)
    current64 = cpp.Spharmt(nlon, nlat, gridtype=gridtype, dtype=np.float64)
    order = nlon // 2
    coefficient_m, coefficient_l = _indices(lmax)
    index = int(np.flatnonzero((coefficient_m == order) & (coefficient_l == lmax))[0])
    latitudes = (
        np.linspace(90.0, -90.0, nlat)
        if gridtype == "regular"
        else cpp.gaussian_lats_wts(nlat)[0]
    )
    longitude_factor = (-1.0) ** np.arange(nlon)
    for imaginary in (False, True):
        spectrum = np.zeros((lmax + 1) * (lmax + 2) // 2, dtype=np.complex64)
        spectrum[index] = 1j if imaginary else 1
        spectrum64 = spectrum.astype(np.complex128)
        reference_grid = reference.spectogrd(spectrum)
        cpp_grid = current.spectogrd(spectrum)
        cpp64_grid = current64.spectogrd(spectrum64)
        expected_values = np.zeros((nlat, nlon), dtype=np.float64)
        if not imaginary:
            expected_values = np.stack(
                [
                    2.0 * _decimal_pbar(float(latitude), order, lmax) * longitude_factor
                    for latitude in latitudes
                ]
            )
        parity.assert_close(
            f"scalar_nyquist_synthesis_{gridtype}_float32",
            "float32",
            cpp_grid,
            reference_grid,
            rtol=3e-6,
            atol=2e-6,
        )
        parity.assert_close(
            f"scalar_nyquist_reference_analytic_{gridtype}_float32",
            "float32",
            reference_grid,
            expected_values,
            rtol=3e-6,
            atol=2e-6,
        )
        parity.assert_close(
            f"scalar_nyquist_cpp_analytic_{gridtype}_float32",
            "float32",
            cpp_grid,
            expected_values,
            rtol=3e-6,
            atol=2e-6,
        )
        expected_values64 = expected_values.astype(np.float64)
        parity.assert_close(
            f"scalar_nyquist_cpp_analytic_{gridtype}_float64",
            "float64",
            cpp64_grid,
            expected_values64,
            rtol=2e-13,
            atol=2e-14,
        )
        reference_from_reference_grid = reference.grdtospec(reference_grid, lmax)
        reference_from_cpp_grid = reference.grdtospec(cpp_grid, lmax)
        cpp_from_reference_grid = current.grdtospec(reference_grid, lmax)
        cpp_from_cpp_grid = current.grdtospec(cpp_grid, lmax)
        cpp64_from_cpp64_grid = current64.grdtospec(cpp64_grid, lmax)
        expected_spectrum = np.zeros_like(spectrum)
        if not imaginary:
            expected_spectrum[index] = 1.0
        for name, recovered in (
            ("reference_from_reference_grid", reference_from_reference_grid),
            ("reference_from_cpp_grid", reference_from_cpp_grid),
            ("cpp_from_reference_grid", cpp_from_reference_grid),
            ("cpp_from_cpp_grid", cpp_from_cpp_grid),
        ):
            parity.assert_close(
                f"scalar_nyquist_{gridtype}_{name}_{'imaginary' if imaginary else 'real'}_float32",
                "float32",
                recovered,
                expected_spectrum,
                rtol=4e-6,
                atol=4e-6,
            )
        expected_spectrum64 = expected_spectrum.astype(np.complex128)
        parity.assert_close(
            f"scalar_nyquist_{gridtype}_cpp64_recovery_{'imaginary' if imaginary else 'real'}",
            "float64",
            cpp64_from_cpp64_grid,
            expected_spectrum64,
            rtol=2e-13,
            atol=2e-14,
        )


@pytest.mark.contract("legendre_many")
@pytest.mark.contract("specintrp_many")
def test_additive_batch_apis_preserve_batch_and_coefficient_axes():
    ntrunc = 5
    latitudes = np.asarray([-53.0, 0.0, 62.0], dtype=np.float64)
    batched_pbar = cpp.legendre_many(latitudes, ntrunc)
    scalar_pbar = np.stack([cpp.legendre(float(lat), ntrunc) for lat in latitudes])
    np.testing.assert_allclose(batched_pbar, scalar_pbar, rtol=0, atol=2e-14)
    spectrum = np.zeros((ntrunc + 1) * (ntrunc + 2) // 2, dtype=np.complex128)
    spectrum[: ntrunc + 1] = np.arange(ntrunc + 1)
    longitudes = np.asarray([-20.0, 45.0, 365.0])
    values = cpp.specintrp_many(longitudes, spectrum, batched_pbar)
    repeated = np.asarray(
        [
            cpp.specintrp(lon, spectrum, row)
            for lon, row in zip(longitudes, batched_pbar)
        ]
    )
    np.testing.assert_allclose(values, repeated, rtol=0, atol=2e-13)


def pytest_collection_finish(session: pytest.Session) -> None:
    seen = {
        str(name)
        for item in session.items
        for marker in item.iter_markers("contract")
        for name in marker.args
    }
    difference_tests = {
        str(name)
        for item in session.items
        for marker in item.iter_markers("approved_difference")
        for name in marker.args
    }
    missing_contracts = set(PUBLIC_CONTRACTS) - seen
    missing_differences = set(APPROVED_DIFFERENCES) - difference_tests
    if missing_contracts or missing_differences:
        pytest.exit(
            "Parity manifest is incomplete. Missing public contracts: "
            f"{sorted(missing_contracts)}; missing approved-difference behavior "
            f"tests: {sorted(missing_differences)}",
            returncode=pytest.ExitCode.TESTS_FAILED,
        )
