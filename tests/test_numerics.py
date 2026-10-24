# SPDX-License-Identifier: GPL-2.0-or-later
import math
from decimal import Decimal, localcontext

import numpy as np
import pytest

import spharm_cpp as spharm


def test_gaussian_rule_and_spectral_indices():
    for nlat in range(1, 13):
        lat, weight = spharm.gaussian_lats_wts(nlat)
        assert lat.shape == weight.shape == (nlat,)
        assert np.all(np.diff(lat) < 0)
        assert np.all(weight > 0)
        assert weight.sum() == pytest.approx(2.0, abs=2e-14)
    m, n = spharm.getspecindx(3)
    assert m.dtype == n.dtype == np.dtype(np.intp)
    m0, n0 = spharm.getspecindx(0)
    assert m0.shape == n0.shape == (1,)
    assert m0.dtype == n0.dtype == np.dtype(np.intp)
    np.testing.assert_array_equal(m0, [0])
    np.testing.assert_array_equal(n0, [0])
    assert list(zip(m.tolist(), n.tolist())) == [
        (0, 0),
        (0, 1),
        (0, 2),
        (0, 3),
        (1, 1),
        (1, 2),
        (1, 3),
        (2, 2),
        (2, 3),
        (3, 3),
    ]


@pytest.mark.parametrize("nlat", [64, 128, 640, 2048])
def test_gaussian_rule_high_resolution(nlat):
    lat, weights = spharm.gaussian_lats_wts(nlat)
    assert lat.dtype == weights.dtype == np.dtype(np.float64)
    assert np.all(np.isfinite(lat)) and np.all(np.isfinite(weights))
    assert np.all(np.diff(lat) < 0)
    assert np.all(weights > 0)
    np.testing.assert_allclose(lat, -lat[::-1], rtol=0, atol=3e-13)
    np.testing.assert_allclose(weights, weights[::-1], rtol=0, atol=3e-15)

    # Independent polynomial moments, exact for Gauss-Legendre rules at nlat >= 3.
    x = np.sin(np.deg2rad(lat))
    for power in (0, 2, 4):
        np.testing.assert_allclose(
            np.dot(weights, x**power), 2.0 / (power + 1), rtol=0, atol=5e-13
        )

    if nlat <= 128:
        nodes, reference_weights = np.polynomial.legendre.leggauss(nlat)
        np.testing.assert_allclose(
            lat, np.rad2deg(np.arcsin(nodes[::-1])), rtol=0, atol=2e-11
        )
        np.testing.assert_allclose(
            weights, reference_weights[::-1], rtol=0, atol=4e-13
        )


def test_legendre_interpolation_and_poles():
    for dtype, cdtype, tol in [
        (np.float32, np.complex64, 5e-6),
        (np.float64, np.complex128, 2e-12),
    ]:
        lmax = 8
        p_north = spharm.legendre(90, lmax, dtype=dtype)
        p_south = spharm.legendre(-90, lmax, dtype=dtype)
        assert p_north.shape == p_south.shape == (45,)
        assert np.all(np.isfinite(p_north)) and np.all(np.isfinite(p_south))
        assert p_north[0] == pytest.approx(1 / np.sqrt(2), abs=tol)
        m, _ = spharm.getspecindx(lmax)
        np.testing.assert_allclose(p_north[m > 0], 0, atol=tol)

        t = spharm.Spharmt(18, 9, gridtype="gaussian", dtype=dtype)
        rng = np.random.default_rng(3)
        spec = (rng.normal(size=45) + 1j * rng.normal(size=45)).astype(cdtype)
        spec[0] = spec[0].real
        grid = t.spectogrd(spec)
        lats = spharm.gaussian_lats_wts(9)[0]
        for j, lat in enumerate(lats):
            pnm = spharm.legendre(float(lat), lmax, dtype=dtype)
            for i in (0, 4, 11):
                lon = i * 360.0 / t.nlon
                assert spharm.specintrp(lon, spec, pnm) == pytest.approx(
                    grid[j, i], abs=tol * 4, rel=tol * 4
                )
    high32 = spharm.legendre(89.999, 512, dtype=np.float32)
    high64 = spharm.legendre(89.999, 512, dtype=np.float64)
    assert np.all(np.isfinite(high32)) and np.all(np.isfinite(high64))
    np.testing.assert_array_equal(high32, high64.astype(np.float32))


@pytest.mark.parametrize(
    "dtype, rtol, atol, norm_rtol",
    [
        (np.float32, 3e-6, 2e-7, 5e-5),
        (np.float64, 5e-11, 5e-12, 5e-10),
    ],
)
def test_high_degree_legendre_against_independent_clenshaw(
    dtype, rtol, atol, norm_rtol
):
    lmax = 2048
    m, degree = spharm.getspecindx(lmax)
    reference_coefficients = np.zeros(lmax + 1)
    for latitude in (23.75, 89.999, -90.0):
        actual = spharm.legendre(latitude, lmax, dtype=dtype)
        assert np.all(np.isfinite(actual))
        x = math.sin(math.radians(latitude))
        for ell in (0, 511, 1024, 2048):
            reference_coefficients[ell] = 1.0
            # NumPy evaluates the Legendre series with Clenshaw's recurrence,
            # independent of the implementation's normalized forward recurrence.
            reference = math.sqrt((2 * ell + 1) / 2) * np.polynomial.legendre.legval(
                x, reference_coefficients
            )
            reference_coefficients[ell] = 0.0
            index = np.flatnonzero((m == 0) & (degree == ell))[0]
            assert actual[index] == pytest.approx(reference, rel=rtol, abs=atol)

        for ell in (511, 1024, 2048):
            row = actual[degree == ell]
            norm = row[0] ** 2 + 2 * np.dot(row[1:], row[1:])
            expected_norm = (2 * ell + 1) / 2
            assert norm == pytest.approx(
                expected_norm, rel=norm_rtol, abs=norm_rtol * expected_norm
            )


@pytest.mark.parametrize(
    "dtype, rtol, atol, norm_rtol",
    [
        (np.float32, 3e-6, 2e-7, 5e-5),
        (np.float64, 5e-11, 5e-12, 5e-10),
    ],
)
def test_high_degree_batched_legendre_simd_lanes(dtype, rtol, atol, norm_rtol):
    lmax = 1024
    latitudes = np.array([-90.0, -89.999, -45.0, 0.0, 22.5, 35.0, 89.999, 90.0])
    actual = spharm.legendre_many(latitudes, lmax, dtype=dtype)
    scalar = np.stack(
        [spharm.legendre(float(latitude), lmax, dtype=dtype) for latitude in latitudes]
    )
    np.testing.assert_allclose(actual, scalar, rtol=rtol, atol=atol)

    m, degree = spharm.getspecindx(lmax)
    reference_coefficients = np.zeros(lmax + 1)
    for row, latitude in enumerate(latitudes):
        x = math.sin(math.radians(float(latitude)))
        for ell in (511, 1024):
            reference_coefficients[ell] = 1.0
            reference = math.sqrt((2 * ell + 1) / 2) * np.polynomial.legendre.legval(
                x, reference_coefficients
            )
            reference_coefficients[ell] = 0.0
            index = np.flatnonzero((m == 0) & (degree == ell))[0]
            assert actual[row, index] == pytest.approx(reference, rel=rtol, abs=atol)

        coefficients = actual[row, degree == lmax]
        norm = coefficients[0] ** 2 + 2 * np.dot(coefficients[1:], coefficients[1:])
        expected_norm = (2 * lmax + 1) / 2
        assert norm == pytest.approx(
            expected_norm, rel=norm_rtol, abs=norm_rtol * expected_norm
        )


def _decimal_legendre_reference(latitude, order, degree):
    """High precision normalized Pbar_lm, without the Condon-Shortley phase."""
    phi = float(latitude) * math.pi / 180.0
    x = Decimal.from_float(math.sin(phi))
    sin_theta = Decimal.from_float(math.cos(phi))
    with localcontext() as context:
        context.prec = 100
        pmm = Decimal(1) / Decimal(2).sqrt()
        for m in range(1, order + 1):
            pmm *= sin_theta * (Decimal(2 * m + 1) / Decimal(2 * m)).sqrt()
        if order == degree:
            return float(pmm)

        previous = pmm
        current = x * Decimal(2 * order + 3).sqrt() * pmm
        if order + 1 == degree:
            return float(current)
        for ell in range(order + 2, degree + 1):
            ld = Decimal(ell)
            md = Decimal(order)
            lp, lm = ld + md, ld - md
            a = (((2 * ld - 1) * (2 * ld + 1)) / (lp * lm)).sqrt()
            b = (((2 * ld + 1) * (lp - 1) * (lm - 1)) / ((2 * ld - 3) * lp * lm)).sqrt()
            previous, current = current, a * x * current - b * previous
        return float(current)


@pytest.mark.parametrize(
    "dtype, norm_rtol, coefficient_rtol",
    [
        (np.float32, 3e-6, 5e-6),
        (np.float64, 2e-11, 2e-9),
    ],
)
def test_degree_2048_legendre_norm_across_latitudes(dtype, norm_rtol, coefficient_rtol):
    lmax = 2048
    latitudes = np.array([45, 55, 60, 65, 68, 70, 75, 80], dtype=np.float64)
    orders, degrees = spharm.getspecindx(lmax)
    batch = spharm.legendre_many(latitudes, lmax, dtype=dtype)
    expected_norm = (2 * lmax + 1) / 2

    for row_index, latitude in enumerate(latitudes):
        scalar = spharm.legendre(float(latitude), lmax, dtype=dtype)
        np.testing.assert_allclose(
            batch[row_index],
            scalar,
            rtol=coefficient_rtol,
            atol=2e-7 if dtype is np.float32 else 2e-12,
        )
        for values in (scalar, batch[row_index]):
            coefficients = values[degrees == lmax].astype(np.float64)
            norm = coefficients[0] ** 2 + 2 * np.dot(coefficients[1:], coefficients[1:])
            assert norm == pytest.approx(
                expected_norm, rel=norm_rtol, abs=norm_rtol * expected_norm
            ), f"degree {lmax} norm failed at latitude {latitude}"


@pytest.mark.parametrize(
    "dtype, rtol",
    [
        (np.float32, 5e-6),
        (np.float64, 2e-9),
    ],
)
def test_high_order_legendre_coefficients_against_decimal_reference(dtype, rtol):
    lmax = 2048
    cases = [
        (45.0, 1441),
        (55.0, 1167),
        (60.0, 1100),
        (60.0, 1200),
        (65.0, 870),
        (68.0, 755),
        (70.0, 695),
        (89.9, 16),
        (89.999, 2),
    ]
    latitudes = np.array(sorted({latitude for latitude, _ in cases}))
    batch = spharm.legendre_many(latitudes, lmax, dtype=dtype)
    scalar_rows = {
        latitude: spharm.legendre(latitude, lmax, dtype=dtype) for latitude in latitudes
    }
    orders, degrees = spharm.getspecindx(lmax)
    latitude_rows = {float(latitude): i for i, latitude in enumerate(latitudes)}

    for latitude, order in cases:
        index = np.flatnonzero((orders == order) & (degrees == lmax))[0]
        expected = _decimal_legendre_reference(latitude, order, lmax)
        for values in (scalar_rows[latitude], batch[latitude_rows[latitude]]):
            assert float(values[index]) == pytest.approx(expected, rel=rtol, abs=0.0), (
                f"Pbar({lmax},{order}) failed at latitude {latitude}"
            )


@pytest.mark.parametrize(
    "dtype, cdtype, tol",
    [
        (np.float32, np.complex64, 8e-6),
        (np.float64, np.complex128, 3e-12),
    ],
)
def test_batched_legendre_and_spectral_interpolation(dtype, cdtype, tol):
    lmax = 8
    latitude_storage = np.array([-90, -65, -20, 0, 35, 71, 90], dtype=np.float64)
    latitudes = latitude_storage[::2]
    pnm = spharm.legendre_many(latitudes, lmax, dtype=dtype)
    assert pnm.shape == (4, (lmax + 1) * (lmax + 2) // 2)
    assert pnm.dtype == np.dtype(dtype)
    expected_pnm = np.stack(
        [spharm.legendre(float(latitude), lmax, dtype=dtype) for latitude in latitudes]
    )
    np.testing.assert_allclose(pnm, expected_pnm, rtol=tol, atol=tol)

    rng = np.random.default_rng(15)
    spec = (rng.normal(size=pnm.shape[1]) + 1j * rng.normal(size=pnm.shape[1])).astype(
        cdtype
    )
    spec[0] = spec[0].real
    longitude_storage = np.array(
        [-30, 0, 45, 120, 211, 360, 405, 720], dtype=np.float64
    )
    longitudes = longitude_storage[::2]
    values = spharm.specintrp_many(longitudes, spec, pnm)
    expected = np.array(
        [
            spharm.specintrp(float(lon), spec, pbar)
            for lon, pbar in zip(longitudes, expected_pnm)
        ]
    )
    assert values.dtype == np.float64
    np.testing.assert_allclose(values, expected, rtol=tol, atol=tol)

    with pytest.raises(ValueError, match="rank one"):
        spharm.legendre_many(np.zeros((2, 2)), lmax)
    with pytest.raises(ValueError, match="shape"):
        spharm.specintrp_many(longitudes, spec, pnm[:-1])


@pytest.mark.parametrize(
    "dtype, cdtype, tol",
    [
        (np.float32, np.complex64, 1e-5),
        (np.float64, np.complex128, 5e-12),
    ],
)
def test_batched_interpolation_complex_coefficient_all_eight_lanes(dtype, cdtype, tol):
    lmax = 8
    order = 3
    degree = 8
    index = order * (2 * lmax - order + 3) // 2 + (degree - order)
    latitudes = np.array([-83.0, -61.0, -38.0, -11.0, 17.0, 42.0, 66.0, 88.0])
    longitudes = np.array([-30.0, 0.0, 45.0, 120.0, 211.0, 360.0, 405.0, 720.0])
    pbar = spharm.legendre_many(latitudes, lmax, dtype=dtype)
    coefficient = 0.75 - 0.25j
    spec = np.zeros(pbar.shape[1], dtype=cdtype)
    spec[index] = coefficient

    values = spharm.specintrp_many(longitudes, spec, pbar)
    phase = np.deg2rad(order * longitudes)
    expected = (
        2.0
        * pbar[:, index]
        * (coefficient.real * np.cos(phase) - coefficient.imag * np.sin(phase))
    )
    np.testing.assert_allclose(values, expected, rtol=tol, atol=tol)


@pytest.mark.parametrize(
    "dtype, cdtype, rtol",
    [
        (np.float32, np.complex64, 1e-5),
        (np.float64, np.complex128, 5e-12),
    ],
)
def test_batched_interpolation_high_order_phase_and_periodicity(dtype, cdtype, rtol):
    lmax = 2048
    nalm = (lmax + 1) * (lmax + 2) // 2
    high_m = 1100
    high_index = high_m * (2 * lmax - high_m + 3) // 2 + (lmax - high_m)
    longitudes = np.array(
        [
            0.125,
            360.0 * (2**20) + 0.125,
            90.125,
            90.125 - 360.0 * (2**20),
            180.125,
            180.125 + 360.0 * (2**20),
            -90.125,
            -90.125 - 360.0 * (2**20),
        ]
    )
    latitudes = np.full(longitudes.shape, 60.0)
    pbar = spharm.legendre_many(latitudes, lmax, dtype=dtype)
    spec = np.zeros(nalm, dtype=cdtype)

    def exact_phase_radians(longitude):
        # Reduce the exact binary64 input after multiplying by m. This keeps
        # the independent reference out of libm's large-argument reduction.
        with localcontext() as ctx:
            ctx.prec = 80
            degrees = (Decimal.from_float(float(longitude)) * high_m) % 360
            if degrees < 0:
                degrees += 360
            return float(degrees) * math.pi / 180.0

    phases = list(map(exact_phase_radians, longitudes))
    for test_coefficient in (0.75 - 0.25j, 0.0 + 0.25j, 0.75 + 0.0j):
        spec[high_index] = test_coefficient
        values = spharm.specintrp_many(longitudes, spec, pbar)
        coefficient = complex(spec[high_index])
        expected = np.array(
            [
                2.0
                * float(row[high_index])
                * (
                    coefficient.real * math.cos(phase)
                    - coefficient.imag * math.sin(phase)
                )
                for row, phase in zip(pbar, phases)
            ]
        )
        np.testing.assert_allclose(values, expected, rtol=rtol, atol=1e-20)

        if test_coefficient == 0.75 - 0.25j:
            mixed_values = values.copy()
            mixed_expected = expected

    values = mixed_values
    expected = mixed_expected
    spec[high_index] = 0.75 - 0.25j
    tail_indices = np.array([1, 2, 6])
    tail_values = spharm.specintrp_many(
        longitudes[tail_indices], spec, pbar[tail_indices]
    )
    np.testing.assert_allclose(
        tail_values, expected[tail_indices], rtol=rtol, atol=1e-20
    )
    # Compare individual lanes against the scalar API at reduced longitudes;
    # large-turn cases above use the independent exact-reduction reference.
    for index in (0, 2, 4, 6):
        reduced_longitude = float(np.remainder(longitudes[index], 360.0))
        scalar = spharm.specintrp(reduced_longitude, spec, pbar[index])
        np.testing.assert_allclose(values[index], scalar, rtol=rtol, atol=1e-20)
    # Inputs separated by whole turns produce identical reduced phases,
    # including for high order m where direct m*longitude loses precision.
    np.testing.assert_array_equal(values[[0, 2, 4, 6]], values[[1, 3, 5, 7]])

    del pbar
    latitudes = np.array([45.0, 65.0, 68.0, 70.0, 75.0, 80.0, -60.0, 89.999])
    longitudes = np.array([0.0, 90.0, 180.0, 360.0, -90.0, 720.0, -360.0, 1080.0])
    pbar = spharm.legendre_many(latitudes, lmax, dtype=dtype)
    near_pole_m = 17
    near_pole_index = near_pole_m * (2 * lmax - near_pole_m + 3) // 2 + (
        lmax - near_pole_m
    )
    spec.fill(0)
    spec[near_pole_index] = -0.5 + 0.375j
    values = spharm.specintrp_many(longitudes, spec, pbar)
    reduced_phase = np.remainder(longitudes, 360.0) * math.pi / 180.0
    coefficient = complex(spec[near_pole_index])
    expected = np.array(
        [
            2.0
            * float(row[near_pole_index])
            * (
                coefficient.real * math.cos(near_pole_m * phase)
                - coefficient.imag * math.sin(near_pole_m * phase)
            )
            for row, phase in zip(pbar, reduced_phase)
        ]
    )
    np.testing.assert_allclose(values, expected, rtol=rtol, atol=1e-14)


def test_scalar_interpolation_reduces_large_finite_longitudes():
    lmax = 1
    latitude = 23.5
    order = 1
    degree = 1
    index = order * (2 * lmax - order + 3) // 2 + (degree - order)
    pbar = spharm.legendre(latitude, lmax)
    spec = np.zeros((lmax + 1) * (lmax + 2) // 2, dtype=np.complex128)
    spec[index] = 0.75 - 0.25j
    longitudes = (1e20, -1e20, 1e308, -1e308)

    for longitude in longitudes:
        # Reduce the exact binary64 input in Decimal before converting the
        # small phase to double. This avoids using libm's huge-angle trig
        # reduction as the numerical reference.
        with localcontext() as ctx:
            ctx.prec = 400
            reduced = float(Decimal.from_float(longitude) % Decimal(360))
        phase = math.radians(reduced)
        expected = (
            2.0
            * pbar[index]
            * (spec[index].real * math.cos(phase) - spec[index].imag * math.sin(phase))
        )
        actual = spharm.specintrp(longitude, spec, pbar)
        assert actual == pytest.approx(expected, rel=2e-14, abs=2e-14)
