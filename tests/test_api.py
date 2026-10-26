# SPDX-FileCopyrightText: 2026 Albert M. W. Yau
# SPDX-License-Identifier: GPL-2.0-or-later
import gc
from importlib.metadata import distribution, version

import numpy as np
import pytest

import spharm_cpp as spharm


def test_installed_distribution_metadata():
    dist = distribution("pyspharm-cpp")
    assert spharm.__version__ == version("pyspharm-cpp")
    assert dist.metadata["Name"] == "pyspharm-cpp"
    assert dist.metadata["Version"] == "0.1.0.dev0"
    assert dist.metadata["License-Expression"] == "GPL-2.0-or-later"
    assert len(dist.requires or []) == 1
    assert dist.requires[0].lower().startswith("numpy")


def test_auxiliary_overloads_replace_many_exports():
    assert "legendre" in spharm.__all__
    assert "specintrp" in spharm.__all__
    assert not hasattr(spharm, "legendre_many")
    assert not hasattr(spharm, "specintrp_many")
    assert not hasattr(spharm._spharm_cpp, "legendre_many")
    assert not hasattr(spharm._spharm_cpp, "specintrp_many")


@pytest.mark.parametrize("m", [1, 2, 4, 10, 20])
def test_geodesic_precision_and_ordering(m):
    latitude64, longitude64 = spharm.getgeodesicpts(m)
    latitude32, longitude32 = spharm.getgeodesicpts(m, dtype=np.float32)
    expected_count = 10 * (m - 1) ** 2 + 2
    assert latitude64.shape == longitude64.shape == (expected_count,)
    assert latitude32.shape == longitude32.shape == (expected_count,)
    assert latitude64.dtype == longitude64.dtype == np.float64
    assert latitude32.dtype == longitude32.dtype == np.float32
    for latitude, longitude, tolerance in (
        (latitude64, longitude64, 4e-16),
        (latitude32, longitude32, 2e-7),
    ):
        assert (latitude[0], latitude[1], longitude[0], longitude[1]) == (
            90,
            -90,
            0,
            0,
        )
        phi = np.deg2rad(latitude.astype(np.float64))
        lam = np.deg2rad(longitude.astype(np.float64))
        xyz_norm = np.sqrt(
            (np.cos(phi) * np.cos(lam)) ** 2
            + (np.cos(phi) * np.sin(lam)) ** 2
            + np.sin(phi) ** 2
        )
        np.testing.assert_allclose(xyz_norm, 1.0, rtol=0, atol=tolerance)


@pytest.mark.parametrize("dtype", [np.float32, np.float64])
@pytest.mark.parametrize("m", [1, 2, 32, 128])
def test_geodesic_output_buffers_are_independent(dtype, m):
    latitudes, longitudes = spharm.getgeodesicpts(m, dtype=dtype)
    assert latitudes.flags.c_contiguous and longitudes.flags.c_contiguous
    assert latitudes.flags.writeable and longitudes.flags.writeable
    assert not np.shares_memory(latitudes, longitudes)

    # Ensure all entries were written deterministically, including poles.
    latitudes_again, longitudes_again = spharm.getgeodesicpts(m, dtype=dtype)
    bits = np.uint32 if dtype == np.float32 else np.uint64
    np.testing.assert_array_equal(latitudes.view(bits), latitudes_again.view(bits))
    np.testing.assert_array_equal(longitudes.view(bits), longitudes_again.view(bits))

    # Modifying one returned array must not affect the other.
    expected_longitudes = longitudes.copy()
    latitudes[:] = 0
    np.testing.assert_array_equal(longitudes, expected_longitudes)


def test_constructor_and_regular_analysis_limit():
    t = spharm.Spharmt(16, 9, 2.0, "regular", "computed")
    assert t.dtype == "float64"
    assert (t.nlon, t.nlat, t.rsphere, t.gridtype, t.legfunc) == (
        16,
        9,
        2.0,
        "regular",
        "computed",
    )
    assert "Spharmt" in repr(t)
    t.extra = "allowed"
    assert t.extra == "allowed"
    with pytest.raises((AttributeError, TypeError)):
        t.nlon = 12
    with pytest.raises((AttributeError, TypeError)):
        del t.nlat
    with pytest.raises(ValueError, match="computed"):
        spharm.Spharmt(16, 9, legfunc="stored")
    with pytest.raises(ValueError, match="nlat-2"):
        t.grdtospec(np.zeros((9, 16)), 8)
    # DUCC synthesis supports the highest regular-grid degree.
    assert t.spectogrd(np.zeros(45, dtype=np.complex128)).shape == (9, 16)
    with pytest.raises(ValueError, match="nlon"):
        spharm.Spharmt(3, 9)
    with pytest.raises(ValueError, match="nlat"):
        spharm.Spharmt(16, 2)


def test_transform_and_utility_results_own_their_storage():
    transform = spharm.Spharmt(18, 10, gridtype="gaussian")
    field = np.arange(180, dtype=np.float64).reshape(10, 18)
    spectrum = transform.grdtospec(field, 8)
    spectrum_reference = spectrum.copy()
    assert not np.shares_memory(spectrum, field)
    field[:] = -1
    np.testing.assert_array_equal(spectrum, spectrum_reference)

    grid = transform.spectogrd(spectrum)
    grid_reference = grid.copy()
    assert not np.shares_memory(grid, spectrum)
    spectrum[:] = 0
    np.testing.assert_array_equal(grid, grid_reference)

    latitudes = np.array([-60.0, 0.0, 60.0])
    pnm = spharm.legendre(latitudes, 8)
    pnm_reference = pnm.copy()
    assert not np.shares_memory(pnm, latitudes)
    latitudes[:] = 90
    np.testing.assert_array_equal(pnm, pnm_reference)

    coefficients = np.arange(pnm.shape[1], dtype=np.float64).astype(np.complex128)
    longitudes = np.array([10.0, 90.0, 230.0])
    values = spharm.specintrp(longitudes, coefficients, pnm)
    values_reference = values.copy()
    assert not np.shares_memory(values, longitudes)
    assert not np.shares_memory(values, coefficients)
    assert not np.shares_memory(values, pnm)
    coefficients[:] = 0
    pnm[:] = 0
    longitudes[:] = 0
    del coefficients, pnm, longitudes
    gc.collect()
    np.testing.assert_array_equal(values, values_reference)


def test_validation_readonly_and_shape_errors():
    t = spharm.Spharmt(16, 9)
    readonly = np.ones((9, 16))
    readonly.flags.writeable = False
    assert t.grdtospec(readonly).shape == (36,)
    with pytest.raises(ValueError, match="shape"):
        t.grdtospec(np.ones((8, 16)))
    with pytest.raises(TypeError, match="float64"):
        t.grdtospec(np.ones((9, 16), dtype=np.float32))
    with pytest.raises(ValueError, match="triangular"):
        t.spectogrd(np.zeros(5, dtype=np.complex128))
    with pytest.raises(ValueError, match="inconsistent"):
        spharm.specintrp(0, np.zeros(6, dtype=np.complex128), np.zeros(3))
    with pytest.raises(ValueError, match="rank-one"):
        spharm.specintrp(0, np.zeros((6, 1), dtype=np.complex128), np.zeros(6))
    with pytest.raises(ValueError, match="finite"):
        spharm.legendre(float("nan"), 3)
    with pytest.raises(ValueError, match="at least 1"):
        spharm.getgeodesicpts(0)
    with pytest.raises(TypeError, match="incompatible function arguments"):
        spharm.getgeodesicpts(2, np.float32)
    with pytest.raises(TypeError, match="float32 or numpy.float64"):
        spharm.getgeodesicpts(2, dtype=np.int32)


@pytest.mark.parametrize(
    "kwargs, error, message",
    [
        ({"rsphere": 0.0}, ValueError, "rsphere"),
        ({"rsphere": float("inf")}, ValueError, "rsphere"),
        ({"gridtype": "invalid"}, ValueError, "gridtype"),
        ({"dtype": np.int32}, TypeError, "dtype"),
    ],
)
def test_constructor_rejects_invalid_options(kwargs, error, message):
    with pytest.raises(error, match=message):
        spharm.Spharmt(16, 9, **kwargs)


def test_transform_rejects_invalid_array_shapes():
    transform = spharm.Spharmt(16, 9)
    field = np.zeros((9, 16))
    spectrum = np.zeros(6, dtype=np.complex128)
    empty_batch = np.empty((6, 0), dtype=np.complex128)
    invalid_calls = [
        (lambda: transform.grdtospec(np.empty((9, 16, 0))), "fields"),
        (lambda: transform.spectogrd(empty_batch), "fields"),
        (lambda: transform.getvrtdivspec(field, field[:, :, None]), "same shape"),
        (lambda: transform.getvrtdivspec(field, field[:, :-1]), "same shape"),
        (lambda: transform.getuv(spectrum, spectrum[:, None]), "matching shape"),
        (lambda: transform.getuv(empty_batch, empty_batch), "fields"),
        (lambda: transform.getgrad(empty_batch), "fields"),
        (lambda: transform.specsmooth(field, np.ones(8)), "smooth"),
        (lambda: transform.specsmooth(field, np.ones((9, 1))), "smooth"),
    ]
    for call, message in invalid_calls:
        with pytest.raises(ValueError, match=message):
            call()
