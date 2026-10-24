# SPDX-License-Identifier: GPL-2.0-or-later
import sys
import sysconfig
from concurrent.futures import ThreadPoolExecutor
from threading import Barrier

import numpy as np

import spharm_cpp as spharm


def test_free_threaded_runtime_remains_disabled_after_import():
    if sysconfig.get_config_var("Py_GIL_DISABLED"):
        assert hasattr(sys, "_is_gil_enabled")
        assert sys._is_gil_enabled() is False


def test_cpu_info_uses_ducc_metadata_schema():
    info = spharm.cpu_info()
    assert isinstance(info["architecture"], str) and info["architecture"]
    assert isinstance(info["features"], list)
    assert all(isinstance(feature, str) for feature in info["features"])
    assert isinstance(info["multiarch"], bool)
    assert set(info) == (
        {
            "architecture",
            "features",
            "multiarch",
            "compiled_profiles",
            "available_profiles",
            "configured_limit",
            "active_profile",
        }
        if info["multiarch"]
        else {"architecture", "features", "multiarch"}
    )
    if info["multiarch"]:
        assert info["active_profile"] in info["compiled_profiles"]
        assert info["active_profile"] in info["available_profiles"]


def test_concurrent_transform_instances_are_independent():
    seeds = range(4)

    def roundtrip(seed):
        transform = spharm.Spharmt(18, 10, gridtype="gaussian")
        rng = np.random.default_rng(seed)
        field = rng.normal(size=(10, 18))
        result = transform.spectogrd(transform.grdtospec(field, 8))
        return field, result

    sequential = [roundtrip(seed) for seed in seeds]
    with ThreadPoolExecutor(max_workers=4) as pool:
        concurrent = list(pool.map(roundtrip, seeds))
    for (expected_field, expected), (field, result) in zip(sequential, concurrent):
        np.testing.assert_array_equal(field, expected_field)
        np.testing.assert_allclose(result, expected, rtol=2e-12, atol=2e-12)


def test_concurrent_calls_on_shared_transform_are_independent():
    transform = spharm.Spharmt(96, 48, gridtype="gaussian")
    rng = np.random.default_rng(31)
    fields = [rng.normal(size=(48, 96)) for _ in range(4)]

    def roundtrip(field, start_gate=None):
        if start_gate is not None:
            start_gate.wait(timeout=30)
        return transform.spectogrd(transform.grdtospec(field, 47))

    expected = [roundtrip(field) for field in fields]
    start_gate = Barrier(len(fields))
    with ThreadPoolExecutor(max_workers=4) as pool:
        concurrent = list(pool.map(lambda field: roundtrip(field, start_gate), fields))
    for result, reference in zip(concurrent, expected):
        np.testing.assert_allclose(result, reference, rtol=2e-12, atol=2e-12)


def test_concurrent_public_operations_match_sequential_references():
    transform = spharm.Spharmt(18, 9, gridtype="gaussian")
    target = spharm.Spharmt(16, 8, gridtype="gaussian")
    rng = np.random.default_rng(47)
    field = rng.normal(size=(9, 18))
    u = rng.normal(size=(9, 18))
    v = rng.normal(size=(9, 18))
    smoothing = np.linspace(0.4, 1.0, 9)
    latitudes = np.array([-71.0, -12.0, 35.0, 82.0])
    longitudes = np.array([-30.0, 0.0, 45.0, 120.0])
    pbar = spharm.legendre_many(latitudes, 6)
    coefficients = transform.grdtospec(field, 6)

    def operations():
        spectrum = transform.grdtospec(field, 6)
        synthesized = transform.spectogrd(spectrum)
        vorticity, divergence = transform.getvrtdivspec(u, v, 6)
        winds = transform.getuv(vorticity, divergence)
        gradient = transform.getgrad(spectrum)
        potentials = transform.getpsichi(u, v, 6)
        smoothed = transform.specsmooth(field, smoothing)
        regridded = spharm.regrid(transform, target, field, ntrunc=6)
        gaussian_lats, gaussian_weights = spharm.gaussian_lats_wts(9)
        order_indices, degree_indices = spharm.getspecindx(6)
        geodesic_lats, geodesic_lons = spharm.getgeodesicpts(3)
        return (
            spectrum,
            synthesized,
            vorticity,
            divergence,
            *winds,
            *gradient,
            *potentials,
            smoothed,
            regridded,
            gaussian_lats,
            gaussian_weights,
            order_indices,
            degree_indices,
            spharm.legendre(23.0, 6),
            pbar,
            spharm.specintrp(45.0, coefficients, pbar[2]),
            spharm.specintrp_many(longitudes, coefficients, pbar),
            geodesic_lats,
            geodesic_lons,
            spharm.cpu_info(),
        )

    expected = operations()
    start_gate = Barrier(4)

    def concurrent_operations(_):
        start_gate.wait(timeout=30)
        return operations()

    with ThreadPoolExecutor(max_workers=4) as pool:
        concurrent = list(pool.map(concurrent_operations, range(4)))

    for result in concurrent:
        assert len(result) == len(expected)
        for actual, reference in zip(result, expected):
            if isinstance(reference, np.ndarray):
                if np.issubdtype(reference.dtype, np.integer):
                    np.testing.assert_array_equal(actual, reference)
                else:
                    np.testing.assert_allclose(
                        actual, reference, rtol=2e-12, atol=2e-12
                    )
            else:
                assert actual == reference
