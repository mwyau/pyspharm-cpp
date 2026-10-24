# SPDX-License-Identifier: GPL-2.0-or-later
import os
import subprocess
import sys
import sysconfig
import textwrap
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


def test_threaded_numeric_kernels_and_transforms_match_serial(tmp_path):
    script = textwrap.dedent(
        """
        import sys
        import numpy as np
        import spharm_cpp as spharm

        rng = np.random.default_rng(1701)
        ntrunc = 128
        npoints = 65
        latitudes = np.linspace(-89.0, 89.0, npoints)
        longitudes = np.linspace(-1440.0, 1440.0, npoints)
        output = {}
        for suffix, dtype, cdtype in (
            ("32", np.float32, np.complex64),
            ("64", np.float64, np.complex128),
        ):
            legendre = spharm.legendre_many(latitudes, ntrunc, dtype=dtype)
            ncoeff = legendre.shape[1]
            spectrum = (
                rng.standard_normal(ncoeff) + 1j * rng.standard_normal(ncoeff)
            ).astype(cdtype)
            output[f"legendre_{suffix}"] = legendre
            output[f"interpolation_{suffix}"] = spharm.specintrp_many(
                longitudes, spectrum, legendre
            )

        lmax = 60
        for suffix, dtype in (("32", np.float32), ("64", np.float64)):
            transform = spharm.Spharmt(128, 64, gridtype="gaussian", dtype=dtype)
            grid = rng.standard_normal((64, 128, 4)).astype(dtype)
            east = rng.standard_normal((64, 128, 4)).astype(dtype)
            north = rng.standard_normal((64, 128, 4)).astype(dtype)
            spectrum = transform.grdtospec(grid, lmax)
            vort, div = transform.getvrtdivspec(east, north, lmax)
            output[f"sht_analysis_{suffix}"] = spectrum
            output[f"sht_synthesis_{suffix}"] = transform.spectogrd(spectrum)
            output[f"sht_vorticity_{suffix}"] = vort
            output[f"sht_divergence_{suffix}"] = div
            output[f"sht_winds_{suffix}"] = np.stack(transform.getuv(vort, div))
        np.savez(sys.argv[1], **output)
        """
    )
    configurations = (
        ("default", None, None),
        ("ducc-four", "4", "1"),
        ("omp-four", None, "4"),
        ("all-threads", "0", "1"),
    )
    output_paths = {}
    for name, ducc_threads, omp_threads in configurations:
        env = os.environ.copy()
        if ducc_threads is None:
            env.pop("DUCC0_NUM_THREADS", None)
        else:
            env["DUCC0_NUM_THREADS"] = ducc_threads
        if omp_threads is None:
            env.pop("OMP_NUM_THREADS", None)
        else:
            env["OMP_NUM_THREADS"] = omp_threads
        output_path = tmp_path / f"threaded-{name}.npz"
        completed = subprocess.run(
            [sys.executable, "-c", script, str(output_path)],
            check=False,
            capture_output=True,
            text=True,
            timeout=120,
            env=env,
        )
        assert completed.returncode == 0, (
            f"{name} subprocess failed:\n{completed.stdout}\n{completed.stderr}"
        )
        output_paths[name] = output_path

    with np.load(output_paths["default"]) as reference:
        for name in ("ducc-four", "omp-four", "all-threads"):
            with np.load(output_paths[name]) as actual:
                for key in reference.files:
                    if key.startswith("sht_"):
                        tolerance = 3e-6 if key.endswith("32") else 3e-12
                        np.testing.assert_allclose(
                            actual[key],
                            reference[key],
                            rtol=tolerance,
                            atol=tolerance,
                        )
                    else:
                        np.testing.assert_array_equal(actual[key], reference[key])
