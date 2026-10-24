# SPDX-License-Identifier: GPL-2.0-or-later
"""Reproducible warmed transform and utility benchmarks.

The optional reference packages must already be installed in the environment.
The benchmark does not install dependencies. Run it with one backend process per
thread count so DUCC reads the requested worker limit before its first call.
"""

from __future__ import annotations

import argparse
import csv
import gc
import importlib.util
import os
import platform
import re
import statistics
import sys
import time
import tracemalloc
from pathlib import Path


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--sizes", nargs="+", default=["256x128", "512x256"],
                        help="nlon x nlat resolutions")
    parser.add_argument("--fields", type=int, default=8)
    parser.add_argument("--repeats", type=int, default=10)
    parser.add_argument("--warmups", type=int, default=3)
    parser.add_argument("--threads", type=int, default=1)
    parser.add_argument("--precision", choices=("both", "float32", "float64"),
                        default="both")
    parser.add_argument("--output", type=Path, default=None,
                        help="optional CSV output path")
    return parser.parse_args()


def result_nbytes(value):
    if isinstance(value, (tuple, list)):
        return sum(result_nbytes(item) for item in value)
    return int(getattr(value, "nbytes", 0))


def measure(function, *, repeats, warmups):
    for _ in range(warmups):
        result = function()
    del result
    gc.collect()
    samples = []
    for _ in range(repeats):
        start = time.perf_counter_ns()
        result = function()
        samples.append((time.perf_counter_ns() - start) * 1e-6)
    output_bytes = result_nbytes(result)
    del result
    gc.collect()
    tracemalloc.start()
    result = function()
    _, peak = tracemalloc.get_traced_memory()
    tracemalloc.stop()
    output_bytes = max(output_bytes, result_nbytes(result))
    del result
    ordered = sorted(samples)
    p90 = ordered[min(len(ordered) - 1, int(0.9 * len(ordered)))]
    return statistics.median(samples), p90, output_bytes, peak


def ducc_provenance():
    cmake_file = Path(__file__).resolve().parents[1] / "CMakeLists.txt"
    if not cmake_file.exists():
        return "", ""
    source = cmake_file.read_text(encoding="utf-8")

    def value(name):
        match = re.search(rf'^set\({name} "([^"]+)"\)$', source, re.MULTILINE)
        return match.group(1) if match else ""

    return value("SPHARM_CPP_DUCC_TAG"), value("SPHARM_CPP_DUCC_REVISION")


def main():
    args = parse_args()
    if args.repeats < 1 or args.warmups < 0 or args.fields < 1 or args.threads < 0:
        raise SystemExit("repeats and fields must be positive; warmups and threads nonnegative")
    os.environ["DUCC0_NUM_THREADS"] = str(args.threads)

    import numpy as np
    import spharm_cpp

    try:
        import ducc0
    except ImportError:
        ducc0 = None
    try:
        import spharm
    except ImportError:
        spharm = None

    dtypes = []
    if args.precision in ("both", "float32"):
        dtypes.append(("float32", np.float32, np.complex64))
    if args.precision in ("both", "float64"):
        dtypes.append(("float64", np.float64, np.complex128))

    fields = ["backend", "version", "active_profile", "compiled_profiles",
              "cpu_architecture", "cpu_features", "ducc_tag", "ducc_revision",
              "dtype", "nlon", "nlat", "nfields",
              "batch_size", "operation", "threads", "median_ms", "p90_ms",
              "gridpoints_per_s", "result_bytes", "python_traced_peak_bytes"]
    output = args.output.open("w", newline="") if args.output else sys.stdout
    writer = csv.DictWriter(output, fieldnames=fields)
    writer.writeheader()

    profiles = spharm_cpp.cpu_info()
    ducc_tag, ducc_revision = ducc_provenance()
    print(f"# platform={platform.platform()} python={platform.python_version()} "
          f"cpu_info={profiles} DUCC0_NUM_THREADS={args.threads}", file=sys.stderr)
    if spharm is None:
        print("# PySpharm reference not installed; skipping it", file=sys.stderr)
    if ducc0 is None:
        print("# ducc0 reference not installed; skipping it", file=sys.stderr)
    if importlib.util.find_spec("pyspharm_rs") is None:
        print("# pyspharm-rs not installed; skipping it", file=sys.stderr)

    backends = [("pyspharm-cpp", spharm_cpp.__version__, spharm_cpp)]
    if spharm is not None:
        backends.append(("PySpharm", getattr(spharm, "__version__", "unknown"), spharm))
    if ducc0 is not None:
        backends.append(("ducc0", ducc0.__version__, ducc0))

    for size in args.sizes:
        try:
            nlon, nlat = (int(part) for part in size.lower().split("x", 1))
        except (ValueError, TypeError):
            raise SystemExit(f"invalid resolution {size!r}; expected nlonxnlat")
        if nlon < 4 or nlat < 3:
            raise SystemExit(f"resolution {size!r} is too small")
        lmax = nlat - 1
        geometry = "GL"
        rng = np.random.default_rng(20261023)

        for label, dtype, cdtype in dtypes:
            scalar = rng.normal(size=(nlat, nlon)).astype(dtype)
            batch = rng.normal(size=(nlat, nlon, args.fields)).astype(dtype)
            u = rng.normal(size=(nlat, nlon)).astype(dtype)
            v = rng.normal(size=(nlat, nlon)).astype(dtype)

            if spharm is not None and label == "float32":
                transform = spharm.Spharmt(nlon, nlat, gridtype="gaussian")
                legacy_backends = [("PySpharm", transform)]
            else:
                legacy_backends = []

            cpp_transform = spharm_cpp.Spharmt(
                nlon, nlat, gridtype="gaussian", dtype=dtype
            )
            package_backends = [("pyspharm-cpp", cpp_transform)]
            package_backends.extend(legacy_backends)

            for backend_name, transform in package_backends:
                if backend_name == "PySpharm":
                    truncation = nlat - 1
                else:
                    truncation = lmax
                coefficients = transform.grdtospec(scalar, truncation)
                batch_coefficients = (
                    transform.grdtospec(batch, truncation)
                    if backend_name == "pyspharm-cpp" else None
                )
                wind_specs = transform.getvrtdivspec(u, v, truncation)
                smooth = np.ones(nlat, dtype=dtype)
                target_nlat = max(3, nlat // 2)
                target_nlon = max(4, nlon // 2)
                if backend_name == "PySpharm":
                    target_transform = spharm.Spharmt(
                        target_nlon, target_nlat, gridtype="regular"
                    )
                    regrid = lambda: spharm.regrid(
                        transform, target_transform, scalar,
                        min(truncation, target_nlat - 2),
                    )
                    legendre = lambda: spharm.legendre(23.5, truncation)
                    geodesic = lambda: spharm.getgeodesicpts(32)
                    interpolate = lambda: spharm.specintrp(
                        41.0, coefficients if coefficients.ndim == 1 else coefficients[:, 0],
                        spharm.legendre(23.5, truncation),
                    )
                    get_indices = lambda: spharm.getspecindx(truncation)
                    gaussian_rule = lambda: spharm.gaussian_lats_wts(nlat)
                else:
                    target_transform = spharm_cpp.Spharmt(
                        target_nlon, target_nlat, gridtype="regular", dtype=dtype
                    )
                    regrid = lambda: spharm_cpp.regrid(
                        transform, target_transform, scalar,
                        min(truncation, target_nlat - 2),
                    )
                    legendre = lambda: spharm_cpp.legendre(
                        23.5, truncation, dtype=dtype
                    )
                    geodesic = lambda: spharm_cpp.getgeodesicpts(32)
                    interpolate = lambda: spharm_cpp.specintrp(
                        41.0, coefficients if coefficients.ndim == 1 else coefficients[:, 0],
                        spharm_cpp.legendre(23.5, truncation, dtype=dtype),
                    )
                    get_indices = lambda: spharm_cpp.getspecindx(truncation)
                    gaussian_rule = lambda: spharm_cpp.gaussian_lats_wts(nlat)

                operations = [
                    ("scalar_analysis", lambda t=transform, x=scalar, L=truncation:
                     t.grdtospec(x, L), 1),
                    ("scalar_synthesis", lambda t=transform, x=coefficients:
                     t.spectogrd(x), 1),
                    ("wind_analysis", lambda t=transform, x=u, y=v, L=truncation:
                     t.getvrtdivspec(x, y, L), 1),
                    ("wind_synthesis", lambda t=transform, x=wind_specs:
                     t.getuv(*x), 1),
                    ("multi_field_analysis", lambda t=transform, x=batch, L=truncation:
                     t.grdtospec(x, L), args.fields),
                    ("specsmooth", lambda t=transform, x=scalar, s=smooth:
                     t.specsmooth(x, s), 1),
                    ("regrid", regrid, 1),
                    ("legendre", legendre, 0),
                    ("geodesic_points", geodesic, 0),
                    ("specintrp", interpolate, 0),
                    ("getspecindx", get_indices, 0),
                    ("gaussian_lats_wts", gaussian_rule, 0),
                ]
                if backend_name == "pyspharm-cpp":
                    # Compare the public batched API with repeated scalar
                    # calls on field views; each uses the same DUCC kernel.
                    scalar_batch_coefficients = np.stack([
                        transform.grdtospec(batch[:, :, field], truncation)
                        for field in range(args.fields)
                    ], axis=1)
                    np.testing.assert_allclose(batch_coefficients,
                                               scalar_batch_coefficients)
                    batch_grid = transform.spectogrd(batch_coefficients)
                    scalar_batch_grid = np.stack([
                        transform.spectogrd(batch_coefficients[:, field])
                        for field in range(args.fields)
                    ], axis=2)
                    np.testing.assert_allclose(batch_grid, scalar_batch_grid)
                    operations.extend([
                        ("multi_field_analysis_scalar_loop", lambda t=transform,
                         x=batch, L=truncation: np.stack([
                             t.grdtospec(x[:, :, field], L)
                             for field in range(args.fields)
                         ], axis=1), args.fields),
                        ("multi_field_synthesis", lambda t=transform,
                         x=batch_coefficients: t.spectogrd(x), args.fields),
                        ("multi_field_synthesis_scalar_loop", lambda t=transform,
                         x=batch_coefficients: np.stack([
                             t.spectogrd(x[:, field])
                             for field in range(args.fields)
                         ], axis=2), args.fields),
                    ])
                    batch_latitudes = np.linspace(-85.0, 85.0, args.fields)
                    batch_longitudes = np.linspace(0.0, 360.0, args.fields, endpoint=False)
                    batch_legendre = spharm_cpp.legendre_many(
                        batch_latitudes, truncation, dtype=dtype
                    )
                    scalar_coefficients = coefficients[:, 0] if coefficients.ndim == 2 else coefficients
                    utility_tolerance = 8e-6 if label == "float32" else 3e-12
                    scalar_legendre_batch = np.stack([
                        spharm_cpp.legendre(float(latitude), truncation,
                                            dtype=dtype)
                        for latitude in batch_latitudes
                    ])
                    np.testing.assert_allclose(
                        batch_legendre, scalar_legendre_batch,
                        rtol=utility_tolerance, atol=utility_tolerance,
                    )
                    scalar_interpolations = np.array([
                        spharm_cpp.specintrp(float(longitude),
                                             scalar_coefficients, pbar)
                        for longitude, pbar in zip(batch_longitudes,
                                                   batch_legendre)
                    ])
                    batch_interpolations = spharm_cpp.specintrp_many(
                        batch_longitudes, scalar_coefficients, batch_legendre
                    )
                    np.testing.assert_allclose(
                        batch_interpolations, scalar_interpolations,
                        rtol=utility_tolerance, atol=utility_tolerance,
                    )
                    operations.extend([
                        ("legendre_scalar_batch", lambda: np.stack([
                            spharm_cpp.legendre(float(latitude), truncation,
                                                dtype=dtype)
                            for latitude in batch_latitudes
                        ]), 0),
                        ("specintrp_scalar_batch", lambda: np.fromiter((
                            spharm_cpp.specintrp(float(longitude),
                                                  scalar_coefficients, pbar)
                            for longitude, pbar in zip(batch_longitudes,
                                                       batch_legendre)
                        ), dtype=np.float64, count=args.fields), 0),
                        ("legendre_many", lambda: spharm_cpp.legendre_many(
                            batch_latitudes, truncation, dtype=dtype), 0),
                        ("specintrp_many", lambda: spharm_cpp.specintrp_many(
                            batch_longitudes, scalar_coefficients, batch_legendre), 0),
                    ])
                for name, function, nfields in operations:
                    median, p90, nbytes, peak = measure(
                        function, repeats=args.repeats, warmups=args.warmups
                    )
                    points = nlat * nlon * max(1, nfields)
                    writer.writerow({
                        "backend": backend_name,
                        "version": (spharm_cpp.__version__ if backend_name == "pyspharm-cpp"
                                    else spharm.__version__),
                        "active_profile": profiles.get("active_profile", ""),
                        "compiled_profiles": ";".join(
                            profiles.get("compiled_profiles", [])
                        ),
                        "cpu_architecture": profiles["architecture"],
                        "cpu_features": ";".join(profiles["features"]),
                        "ducc_tag": ducc_tag,
                        "ducc_revision": ducc_revision,
                        "dtype": label,
                        "nlon": nlon,
                        "nlat": nlat,
                        "nfields": nfields,
                        "batch_size": (args.fields if name.endswith("_many")
                                       or name.endswith("_scalar_batch") else nfields),
                        "operation": name,
                        "threads": args.threads,
                        "median_ms": f"{median:.6f}",
                        "p90_ms": f"{p90:.6f}",
                        "gridpoints_per_s": "" if not nfields else f"{points / (median * 1e-3):.3f}",
                        "result_bytes": nbytes,
                        "python_traced_peak_bytes": peak,
                    })

            if ducc0 is not None:
                alm = ducc0.sht.analysis_2d(
                    map=scalar[None, :, :], spin=0, lmax=lmax,
                    geometry=geometry, nthreads=args.threads,
                )
                winds = np.stack((u, v))
                wind_alms = ducc0.sht.analysis_2d(
                    map=winds, spin=1, lmax=lmax,
                    geometry=geometry, nthreads=args.threads,
                )
                operations = [
                    ("scalar_analysis", lambda: ducc0.sht.analysis_2d(
                        map=scalar[None, :, :], spin=0, lmax=lmax,
                        geometry=geometry, nthreads=args.threads), 1),
                    ("scalar_synthesis", lambda: ducc0.sht.synthesis_2d(
                        alm=alm, spin=0, lmax=lmax, geometry=geometry,
                        ntheta=nlat, nphi=nlon, nthreads=args.threads), 1),
                    ("wind_analysis", lambda: ducc0.sht.analysis_2d(
                        map=winds, spin=1, lmax=lmax, geometry=geometry,
                        nthreads=args.threads), 1),
                    ("wind_synthesis", lambda: ducc0.sht.synthesis_2d(
                        alm=wind_alms, spin=1, lmax=lmax, geometry=geometry,
                        ntheta=nlat, nphi=nlon, nthreads=args.threads), 1),
                    ("multi_field_analysis", lambda: [ducc0.sht.analysis_2d(
                        map=batch[:, :, i][None, :, :], spin=0, lmax=lmax,
                        geometry=geometry, nthreads=args.threads)
                        for i in range(args.fields)], args.fields),
                ]
                for name, function, nfields in operations:
                    median, p90, nbytes, peak = measure(
                        function, repeats=args.repeats, warmups=args.warmups
                    )
                    points = nlat * nlon * max(1, nfields)
                    writer.writerow({
                        "backend": "ducc0", "version": ducc0.__version__,
                        "dtype": label, "nlon": nlon, "nlat": nlat,
                        "nfields": nfields,
                        "batch_size": nfields if name == "multi_field_analysis" else 0,
                        "operation": name,
                        "threads": args.threads,
                        "median_ms": f"{median:.6f}", "p90_ms": f"{p90:.6f}",
                        "gridpoints_per_s": f"{points / (median * 1e-3):.3f}",
                        "result_bytes": nbytes,
                        "python_traced_peak_bytes": peak,
                    })

    if args.output:
        output.close()


if __name__ == "__main__":
    main()
