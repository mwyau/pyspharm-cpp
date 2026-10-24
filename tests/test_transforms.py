# SPDX-FileCopyrightText: 2026 Albert M. W. Yau
# SPDX-License-Identifier: GPL-2.0-or-later
import math

import numpy as np
import pytest

import spharm_cpp as spharm


@pytest.mark.parametrize(
    "dtype, cdtype, atol",
    [
        (np.float32, np.complex64, 3e-5),
        (np.float64, np.complex128, 2e-12),
    ],
)
@pytest.mark.parametrize("gridtype", ["regular", "gaussian"])
def test_scalar_roundtrip(dtype, cdtype, atol, gridtype):
    nlon, nlat = 18, 9
    lmax = nlat - 2 if gridtype == "regular" else nlat - 1
    transform = spharm.Spharmt(nlon, nlat, 2.0, gridtype, "computed", dtype=dtype)
    rng = np.random.default_rng(1)
    nalm = (lmax + 1) * (lmax + 2) // 2
    spec = (rng.normal(size=nalm) + 1j * rng.normal(size=nalm)).astype(cdtype)
    spec[: lmax + 1] = spec[: lmax + 1].real

    grid = transform.spectogrd(spec)
    recovered = transform.grdtospec(grid, lmax)
    np.testing.assert_allclose(recovered, spec, rtol=atol, atol=atol)
    assert grid.shape == (nlat, nlon)
    assert grid.dtype == np.dtype(dtype)


@pytest.mark.parametrize(
    "dtype, cdtype, atol",
    [
        (np.float32, np.complex64, 6e-5),
        (np.float64, np.complex128, 3e-12),
    ],
)
@pytest.mark.parametrize("nlon, nlat", [(16, 10), (17, 9), (16, 9), (17, 10)])
@pytest.mark.parametrize("gridtype", ["regular", "gaussian"])
def test_odd_even_grid_roundtrips(dtype, cdtype, atol, nlon, nlat, gridtype):
    t = spharm.Spharmt(nlon, nlat, gridtype=gridtype, dtype=dtype)
    lmax = min(5, nlat - 2 if gridtype == "regular" else nlat - 1)
    rng = np.random.default_rng(nlon * 31 + nlat)
    nalm = (lmax + 1) * (lmax + 2) // 2
    spec = (rng.normal(size=nalm) + 1j * rng.normal(size=nalm)).astype(cdtype)
    spec[: lmax + 1] = spec[: lmax + 1].real

    grid = t.spectogrd(spec)
    recovered = t.grdtospec(grid, lmax)
    np.testing.assert_allclose(recovered, spec, rtol=atol, atol=atol)


@pytest.mark.parametrize(
    "dtype, cdtype, atol",
    [
        (np.float32, np.complex64, 3e-5),
        (np.float64, np.complex128, 2e-12),
    ],
)
@pytest.mark.parametrize("gridtype", ["regular", "gaussian"])
def test_even_longitude_nyquist_mode(dtype, cdtype, atol, gridtype):
    nlon, nlat, lmax = 16, 10, 8
    t = spharm.Spharmt(nlon, nlat, gridtype=gridtype, dtype=dtype)
    m, n = spharm.getspecindx(lmax)
    nyquist_index = np.flatnonzero((m == nlon // 2) & (n == nlon // 2))[0]
    spec = np.zeros((lmax + 1) * (lmax + 2) // 2, dtype=cdtype)
    spec[nyquist_index] = 1

    grid = t.spectogrd(spec)
    latitudes = (
        np.linspace(90.0, -90.0, nlat)
        if gridtype == "regular"
        else spharm.gaussian_lats_wts(nlat)[0]
    )
    expected = np.stack(
        [
            2
            * spharm.legendre(float(latitude), lmax, dtype=dtype)[nyquist_index]
            * (-1.0) ** np.arange(nlon)
            for latitude in latitudes
        ]
    ).astype(dtype)
    np.testing.assert_allclose(grid, expected, rtol=atol, atol=atol)
    np.testing.assert_allclose(
        grid[:, ::2],
        np.broadcast_to(grid[:, :1], grid[:, ::2].shape),
        rtol=atol,
        atol=atol,
    )
    np.testing.assert_allclose(
        grid[:, 1::2],
        np.broadcast_to(-grid[:, :1], grid[:, 1::2].shape),
        rtol=atol,
        atol=atol,
    )
    analyzed = t.grdtospec(grid, lmax)
    assert analyzed[nyquist_index].real == pytest.approx(1, rel=atol, abs=atol)
    assert np.max(np.abs(np.delete(analyzed, nyquist_index))) < atol


@pytest.mark.parametrize("gridtype", ["regular", "gaussian"])
@pytest.mark.parametrize("nlon", [8, 9])
@pytest.mark.parametrize(
    "dtype, cdtype",
    [
        (np.float32, np.complex64),
        (np.float64, np.complex128),
    ],
)
def test_analysis_truncation_respects_longitude_resolution(
    gridtype, nlon, dtype, cdtype
):
    nlat = 12
    transform = spharm.Spharmt(nlon, nlat, gridtype=gridtype, dtype=dtype)
    longitude_limit = nlon // 2
    count = (longitude_limit + 1) * (longitude_limit + 2) // 2
    vector_limit = (nlon - 1) // 2
    vector_count = (vector_limit + 1) * (vector_limit + 2) // 2
    grid = np.zeros((nlat, nlon), dtype=dtype)

    default_spectrum = transform.grdtospec(grid)
    assert default_spectrum.shape == (count,)
    assert default_spectrum.dtype == np.dtype(cdtype)
    assert transform.grdtospec(grid, longitude_limit).shape == (count,)
    default_vorticity, default_divergence = transform.getvrtdivspec(grid, grid)
    assert default_vorticity.shape == default_divergence.shape == (vector_count,)
    vorticity, divergence = transform.getvrtdivspec(grid, grid, vector_limit)
    assert vorticity.shape == divergence.shape == (vector_count,)
    streamfunction, potential = transform.getpsichi(grid, grid, vector_limit)
    assert streamfunction.shape == potential.shape == (nlat, nlon)
    target = spharm.Spharmt(16, 8, gridtype=gridtype, dtype=dtype)
    assert spharm.regrid(transform, target, grid).shape == (8, 16)

    with pytest.raises(ValueError, match="nlon/2"):
        transform.grdtospec(grid, longitude_limit + 1)
    with pytest.raises(ValueError, match=r"floor\(\(nlon-1\)/2\)"):
        transform.getvrtdivspec(grid, grid, vector_limit + 1)
    with pytest.raises(ValueError, match=r"floor\(\(nlon-1\)/2\)"):
        transform.getpsichi(grid, grid, vector_limit + 1)
    with pytest.raises(ValueError, match="nlon/2"):
        spharm.regrid(transform, target, grid, longitude_limit + 1)


def _diagonal_harmonic_gradient(
    order: int,
    coefficient: complex,
    latitudes: np.ndarray,
    longitudes: np.ndarray,
    radius: float,
) -> tuple[np.ndarray, np.ndarray]:
    """Analytic gradient for the normalized no-Condon–Shortley mode l=m."""
    normalization = math.sqrt(
        math.factorial(2 * order)
        * (2 * order + 1)
        / (2 ** (2 * order + 1) * math.factorial(order) ** 2)
    )
    phi = np.deg2rad(latitudes)[:, None]
    lam = np.deg2rad(longitudes)[None, :]
    cos_power = np.cos(phi) ** (order - 1)
    sin_mlam = np.sin(order * lam)
    cos_mlam = np.cos(order * lam)
    east = (
        -2
        * order
        * normalization
        * cos_power
        / radius
        * (coefficient.real * sin_mlam + coefficient.imag * cos_mlam)
    )
    north = (
        -2
        * order
        * normalization
        * np.sin(phi)
        * cos_power
        / radius
        * (coefficient.real * cos_mlam - coefficient.imag * sin_mlam)
    )
    return east, north


@pytest.mark.parametrize("gridtype", ["regular", "gaussian"])
@pytest.mark.parametrize("order", [7, 8])
@pytest.mark.parametrize("coefficient", [1.0 + 0.0j, 0.0 + 1.0j])
def test_spin_one_synthesis_matches_analytic_modes_at_and_below_nyquist(
    gridtype, order, coefficient
):
    nlat, nlon = 10, 16
    transform = spharm.Spharmt(nlon, nlat, gridtype=gridtype, dtype=np.float64)
    orders, degrees = spharm.getspecindx(order)
    index = int(np.flatnonzero((orders == order) & (degrees == order))[0])
    spectrum = np.zeros((order + 1) * (order + 2) // 2, dtype=np.complex128)
    spectrum[index] = coefficient
    zero = np.zeros_like(spectrum)
    latitudes = (
        np.linspace(90.0, -90.0, nlat)
        if gridtype == "regular"
        else spharm.gaussian_lats_wts(nlat)[0]
    )
    longitudes = np.arange(nlon, dtype=np.float64) * (360.0 / nlon)

    expected_east, expected_north = _diagonal_harmonic_gradient(
        order, coefficient, latitudes, longitudes, transform.rsphere
    )
    east, north = transform.getgrad(spectrum)
    np.testing.assert_allclose(east, expected_east, rtol=3e-13, atol=2e-19)
    np.testing.assert_allclose(north, expected_north, rtol=3e-13, atol=2e-19)

    eigenvalue = order * (order + 1)
    potential_coefficient = -(transform.rsphere**2) * coefficient / eigenvalue
    potential_east, potential_north = _diagonal_harmonic_gradient(
        order,
        potential_coefficient,
        latitudes,
        longitudes,
        transform.rsphere,
    )
    vorticity_east, vorticity_north = transform.getuv(spectrum, zero)
    divergence_east, divergence_north = transform.getuv(zero, spectrum)
    np.testing.assert_allclose(divergence_east, potential_east, rtol=3e-13, atol=1e-8)
    np.testing.assert_allclose(divergence_north, potential_north, rtol=3e-13, atol=1e-8)
    np.testing.assert_allclose(vorticity_east, -potential_north, rtol=3e-13, atol=1e-8)
    np.testing.assert_allclose(vorticity_north, potential_east, rtol=3e-13, atol=1e-8)

    if order == nlon // 2 and coefficient.imag:
        scalar_samples = transform.spectogrd(spectrum)
        assert np.max(np.abs(scalar_samples)) == 0
        assert np.max(np.abs(east)) > 0


@pytest.mark.parametrize(
    "dtype, cdtype, atol",
    [
        (np.float32, np.complex64, 2e-4),
        (np.float64, np.complex128, 2e-12),
    ],
)
def test_batched_noncontiguous_scalar_and_wind_transforms(dtype, cdtype, atol):
    t = spharm.Spharmt(16, 9, 3.0, "gaussian", dtype=dtype)
    rng = np.random.default_rng(2)
    raw = rng.normal(size=(9, 16, 4)).astype(dtype)
    maps = raw[:, :, ::-2]
    assert not maps.flags.c_contiguous
    specs = t.grdtospec(maps, 5)
    assert specs.shape == (21, 2)
    recovered = t.spectogrd(specs)
    expected = t.spectogrd(t.grdtospec(np.ascontiguousarray(maps), 5))
    np.testing.assert_allclose(recovered, expected, rtol=atol, atol=atol)

    nalm = 21
    vrt = (rng.normal(size=(nalm, 2)) + 1j * rng.normal(size=(nalm, 2))).astype(cdtype)
    div = (rng.normal(size=(nalm, 2)) + 1j * rng.normal(size=(nalm, 2))).astype(cdtype)
    vrt[:6] = vrt[:6].real
    div[:6] = div[:6].real
    vrt[0] = 0
    div[0] = 0
    vrt = vrt[:, ::-1]
    u, v = t.getuv(vrt, div)
    vrt2, div2 = t.getvrtdivspec(u, v, 5)
    np.testing.assert_allclose(vrt2, vrt, rtol=atol, atol=atol)
    np.testing.assert_allclose(div2, div, rtol=atol, atol=atol)


@pytest.mark.parametrize(
    "dtype, cdtype, atol",
    [
        (np.float32, np.complex64, 4e-5),
        (np.float64, np.complex128, 2e-12),
    ],
)
def test_gradient_vorticity_potential_and_smoothing(dtype, cdtype, atol):
    t = spharm.Spharmt(18, 9, 2.0, "gaussian", dtype=dtype)
    spec = np.zeros(21, dtype=cdtype)
    spec[1] = 1
    chi_grid = t.spectogrd(spec)
    u, v = t.getgrad(spec)
    vrt, div = t.getvrtdivspec(u, v, 5)
    assert np.max(np.abs(vrt)) < atol
    np.testing.assert_allclose(div[1], -2.0 / (t.rsphere**2), rtol=atol, atol=atol)
    psi, chi = t.getpsichi(u, v, 5)
    assert np.max(np.abs(psi)) < atol
    np.testing.assert_allclose(chi, chi_grid, rtol=atol, atol=atol)
    smooth = np.ones(t.nlat, dtype=dtype)
    np.testing.assert_allclose(
        t.specsmooth(chi_grid, smooth), chi_grid, rtol=atol, atol=atol
    )

    # Independent non-zonal l=1,m=1 gradient harmonic. With the API's
    # no-Condon-Shortley normalization, f=sqrt(3) cos(phi) cos(lambda).
    nonzonal = np.zeros(3, dtype=cdtype)
    nonzonal[2] = 1
    eastward, northward = t.getgrad(nonzonal)
    latitudes, _ = spharm.gaussian_lats_wts(t.nlat)
    longitudes = np.arange(t.nlon, dtype=np.float64) * (360.0 / t.nlon)
    phi = np.deg2rad(latitudes)[:, None]
    lam = np.deg2rad(longitudes)[None, :]
    expected_eastward = np.broadcast_to(
        -np.sqrt(3.0) / t.rsphere * np.sin(lam), eastward.shape
    )
    expected_northward = np.broadcast_to(
        -np.sqrt(3.0) / t.rsphere * np.sin(phi) * np.cos(lam),
        northward.shape,
    )
    np.testing.assert_allclose(eastward, expected_eastward, rtol=atol, atol=atol)
    np.testing.assert_allclose(northward, expected_northward, rtol=atol, atol=atol)
    harmonic_vorticity, harmonic_divergence = t.getvrtdivspec(eastward, northward, 1)
    assert np.max(np.abs(harmonic_vorticity)) < atol
    np.testing.assert_allclose(
        harmonic_divergence[2], -2.0 / t.rsphere**2, rtol=atol, atol=atol
    )


def test_geodesic_points_and_regrid_precision_conversion():
    for m in range(1, 20):
        lat, lon = spharm.getgeodesicpts(m)
        lat32, lon32 = spharm.getgeodesicpts(m, dtype=np.float32)
        assert lat.dtype == lon.dtype == np.float64
        assert lat32.dtype == lon32.dtype == np.float32
        assert lat.shape == lon.shape == (10 * (m - 1) ** 2 + 2,)
        assert (lat[0], lat[1], lon[0], lon[1]) == (90, -90, 0, 0)
        assert np.all(np.isfinite(lat)) and np.all(np.isfinite(lon))
        assert np.all((-90 <= lat) & (lat <= 90))
        assert np.all((0 <= lon) & (lon <= 360))
        np.testing.assert_array_equal(lat, spharm.getgeodesicpts(m)[0])
        np.testing.assert_array_equal(
            lat32, spharm.getgeodesicpts(m, dtype=np.float32)[0]
        )
        np.testing.assert_allclose(lat32, lat, rtol=0, atol=8e-5)
        longitude_difference = (lon32.astype(np.float64) - lon + 180.0) % 360.0 - 180.0
        np.testing.assert_allclose(longitude_difference, 0.0, rtol=0, atol=8e-5)

    source = spharm.Spharmt(16, 9, gridtype="gaussian", dtype=np.float32)
    target = spharm.Spharmt(12, 7, gridtype="regular", dtype=np.float64)
    field = np.ones((9, 16), dtype=np.float32)
    result = spharm.regrid(source, target, field)
    assert result.dtype == np.float64
    np.testing.assert_allclose(result, 1, rtol=2e-7, atol=2e-7)
    smoothed = spharm.regrid(source, target, field, smooth=np.ones(7, dtype=np.float64))
    np.testing.assert_allclose(smoothed, 1, rtol=2e-7, atol=2e-7)


@pytest.mark.parametrize("dtype", [np.float32, np.float64])
def test_strided_readonly_smoothing_factors(dtype):
    transform = spharm.Spharmt(18, 9, gridtype="gaussian", dtype=dtype)
    rng = np.random.default_rng(73)
    grid = rng.normal(size=(9, 18)).astype(dtype)
    smooth = np.linspace(0.3, 1.1, 9, dtype=dtype)
    backing = np.empty(18, dtype=dtype)
    backing[::2] = smooth
    backing[1::2] = -17
    strided = backing[::2]
    strided.flags.writeable = False

    tolerance = 2e-6 if dtype == np.float32 else 2e-13
    expected = transform.specsmooth(grid, smooth)
    np.testing.assert_allclose(
        transform.specsmooth(grid, strided), expected, rtol=tolerance, atol=tolerance
    )

    target = spharm.Spharmt(18, 9, gridtype="gaussian", dtype=dtype)
    expected_regrid = spharm.regrid(transform, target, grid, smooth=smooth)
    actual_regrid = spharm.regrid(transform, target, grid, smooth=strided)
    np.testing.assert_allclose(
        actual_regrid, expected_regrid, rtol=tolerance, atol=tolerance
    )


@pytest.mark.parametrize("gridtype", ["regular", "gaussian"])
def test_float32_radius_scaling_without_premature_casts(gridtype):
    baseline = spharm.Spharmt(16, 9, 1.0, gridtype, dtype=np.float32)
    large = spharm.Spharmt(16, 9, 1e40, gridtype, dtype=np.float32)
    spec = np.zeros(3, dtype=np.complex64)
    spec[1] = 1e20
    zero = np.zeros_like(spec)

    reference_grad = baseline.getgrad(spec)
    large_grad = large.getgrad(spec)
    for actual, reference in zip(large_grad, reference_grad):
        np.testing.assert_allclose(
            actual.astype(np.float64),
            reference.astype(np.float64) / large.rsphere,
            rtol=5e-5,
            atol=1e-26,
        )
    assert np.any(large_grad[1] != 0)

    small_spec = np.zeros_like(spec)
    small_spec[1] = 1e-20
    for actual, reference in zip(
        large.getuv(small_spec, zero), baseline.getuv(spec, zero)
    ):
        np.testing.assert_allclose(
            actual / 1e20, reference / 1e20, rtol=5e-5, atol=5e-6
        )

    huge = spharm.Spharmt(16, 9, 1e50, gridtype, dtype=np.float32)
    _, divergence = huge.getvrtdivspec(*reference_grad, 1)
    _, reference_divergence = baseline.getvrtdivspec(*reference_grad, 1)
    np.testing.assert_allclose(
        divergence.astype(np.complex128),
        reference_divergence.astype(np.complex128) / huge.rsphere,
        rtol=5e-5,
        atol=1e-36,
    )
    assert np.any(divergence != 0)

    medium = spharm.Spharmt(16, 9, 1e23, gridtype, dtype=np.float32)
    winds = tuple(component * np.float32(1e-9) for component in reference_grad)
    reference_potentials = baseline.getpsichi(*winds, 1)
    medium_potentials = medium.getpsichi(*winds, 1)
    for actual, reference in zip(medium_potentials, reference_potentials):
        np.testing.assert_allclose(
            actual.astype(np.float64) / medium.rsphere,
            reference.astype(np.float64),
            rtol=5e-5,
            atol=1e6,
        )
