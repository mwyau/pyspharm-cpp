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


@pytest.mark.parametrize("m", [1, 2, 4, 10])
def test_geodesic_float32_is_rounded_float64(m):
    latitude64, longitude64 = spharm.getgeodesicpts(m)
    latitude32, longitude32 = spharm.getgeodesicpts(m, dtype=np.float32)
    assert latitude64.dtype == longitude64.dtype == np.float64
    assert latitude32.dtype == longitude32.dtype == np.float32
    np.testing.assert_array_equal(latitude32, latitude64.astype(np.float32))
    np.testing.assert_array_equal(longitude32, longitude64.astype(np.float32))


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
    pnm = spharm.legendre_many(latitudes, 8)
    pnm_reference = pnm.copy()
    assert not np.shares_memory(pnm, latitudes)
    latitudes[:] = 90
    np.testing.assert_array_equal(pnm, pnm_reference)

    coefficients = np.arange(pnm.shape[1], dtype=np.float64).astype(np.complex128)
    longitudes = np.array([10.0, 90.0, 230.0])
    values = spharm.specintrp_many(longitudes, coefficients, pnm)
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
