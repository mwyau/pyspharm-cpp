# SPDX-FileCopyrightText: 2026 Albert M. W. Yau
# SPDX-License-Identifier: GPL-2.0-or-later
"""Measure scalar and wind SHT analysis/synthesis on five grid cases."""

from __future__ import annotations

import argparse
import csv
import gc
import json
import math
import os
import platform
import statistics
import time
from importlib.metadata import version as package_version
from pathlib import Path
from typing import TypedDict

# Case name, grid type, longitudes, latitudes, chart label.
CASES = (
    ("cc-2p5", "regular", 144, 73, "Clenshaw–Curtis · 2.5° · 144 × 73"),
    ("cc-2p5-1t", "regular", 144, 73, "Clenshaw–Curtis · 2.5° · 144 × 73"),
    ("cc-0p25", "regular", 1440, 721, "Clenshaw–Curtis · 0.25° · 1440 × 721"),
    ("gl-128", "gaussian", 256, 128, "Gaussian · 128 latitude rings · 256 × 128"),
    ("gl-640", "gaussian", 1280, 640, "Gaussian · 640 latitude rings · 1280 × 640"),
)
OPERATIONS = (
    "scalar_analysis",
    "scalar_synthesis",
    "wind_analysis",
    "wind_synthesis",
)


class CpuResources(TypedDict):
    logical_cpus: int
    affinity_cpus: int
    quota_cpus: float | None
    available_cpus: int


def cpu_resources() -> CpuResources:
    logical_cpus = os.cpu_count() or 1
    if hasattr(os, "sched_getaffinity"):
        affinity_cpus = len(os.sched_getaffinity(0))
    else:
        affinity_cpus = logical_cpus

    quota_cpus = None
    try:
        value, period = Path("/sys/fs/cgroup/cpu.max").read_text().split()
        if value != "max":
            quota_cpus = int(value) / int(period)
    except (OSError, ValueError):
        try:
            quota = int(Path("/sys/fs/cgroup/cpu/cpu.cfs_quota_us").read_text())
            period = int(Path("/sys/fs/cgroup/cpu/cpu.cfs_period_us").read_text())
            if quota > 0:
                quota_cpus = quota / period
        except (OSError, ValueError):
            pass

    available_cpus = affinity_cpus
    if quota_cpus is not None:
        available_cpus = min(available_cpus, math.floor(quota_cpus))

    return {
        "logical_cpus": logical_cpus,
        "affinity_cpus": affinity_cpus,
        "quota_cpus": quota_cpus,
        "available_cpus": available_cpus,
    }


def cpu_model() -> str:
    try:
        for line in Path("/proc/cpuinfo").read_text().splitlines():
            if line.lower().startswith("model name"):
                return line.split(":", 1)[1].strip()
    except OSError:
        pass
    return platform.processor() or "unknown"


def measure(callback, warmups, repeats):
    for _ in range(warmups):
        callback()
    gc.collect()
    samples = []
    for _ in range(repeats):
        start = time.perf_counter_ns()
        callback()
        samples.append((time.perf_counter_ns() - start) * 1e-6)
    samples.sort()
    return statistics.median(samples), samples[math.ceil(0.9 * repeats) - 1], samples


def pyspharm_operations(transform, grid, east, north, lmax):
    # Spectra for synthesis are prepared outside the timed region.
    scalar = transform.grdtospec(grid, lmax)
    vorticity, divergence = transform.getvrtdivspec(east, north, lmax)
    return {
        "scalar_analysis": lambda: transform.grdtospec(grid, lmax),
        "scalar_synthesis": lambda: transform.spectogrd(scalar),
        "wind_analysis": lambda: transform.getvrtdivspec(east, north, lmax),
        "wind_synthesis": lambda: transform.getuv(vorticity, divergence),
    }


def ducc_operations(np, ducc0, grid, east, north, gridtype, lmax, threads):
    geometry = "CC" if gridtype == "regular" else "GL"
    nlat, nlon = grid.shape
    scalar_map = grid[None, :, :]
    wind_map = np.stack((east, north))
    scalar = ducc0.sht.analysis_2d(
        map=scalar_map, spin=0, lmax=lmax, geometry=geometry, nthreads=threads
    )
    wind = ducc0.sht.analysis_2d(
        map=wind_map, spin=1, lmax=lmax, geometry=geometry, nthreads=threads
    )

    def analyze(spin, data):
        return ducc0.sht.analysis_2d(
            map=data, spin=spin, lmax=lmax, geometry=geometry, nthreads=threads
        )

    def synthesize(spin, alms):
        return ducc0.sht.synthesis_2d(
            alm=alms, spin=spin, lmax=lmax, geometry=geometry,
            ntheta=nlat, nphi=nlon, nthreads=threads,
        )

    return {
        "scalar_analysis": lambda: analyze(0, scalar_map),
        "scalar_synthesis": lambda: synthesize(0, scalar),
        "wind_analysis": lambda: analyze(1, wind_map),
        "wind_synthesis": lambda: synthesize(1, wind),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--threads", type=int, default=0)
    parser.add_argument(
        "--cases",
        nargs="+",
        choices=[case[0] for case in CASES],
        default=["cc-2p5", "cc-0p25", "gl-128", "gl-640"],
    )
    parser.add_argument("--warmups", type=int, default=2)
    parser.add_argument("--repeats", type=int, default=7)
    parser.add_argument(
        "--expected-profile",
        choices=("x86-64-v1", "x86-64-v3", "x86-64-v4"),
        help="verify the requested SIMD profile before measurement",
    )
    parser.add_argument(
        "--print-samples",
        action="store_true",
        help="print individual timing samples for diagnosing variability",
    )
    args = parser.parse_args()
    if args.threads < 0 or args.warmups < 1 or args.repeats < 3:
        parser.error("threads >= 0; warmups >= 1; repeats >= 3")
    resources = cpu_resources()
    threads = args.threads or resources["available_cpus"]
    if threads < 1:
        parser.error("no CPUs are available to run the benchmark")
    if threads > resources["available_cpus"]:
        parser.error(
            f"requested {threads} threads but only "
            f"{resources['available_cpus']} CPUs are available from affinity and quota"
        )
    if "cc-2p5-1t" in args.cases and threads != 1:
        parser.error("cc-2p5-1t requires --threads 1")

    # pyspharm-cpp reads DUCC0_NUM_THREADS once, at its first transform.
    for variable in ("DUCC0_NUM_THREADS", "OMP_NUM_THREADS"):
        os.environ[variable] = str(threads)
    os.environ["OPENBLAS_NUM_THREADS"] = "1"

    import numpy as np
    import ducc0
    import spharm as pyspharm
    import spharm_cpp

    active_profile = spharm_cpp.cpu_info()["active_profile"]
    if args.expected_profile and active_profile != args.expected_profile:
        parser.error(f"expected SIMD {args.expected_profile}, got {active_profile}")

    metadata = {
        "threads": threads,
        **resources,
        "cpu_model": cpu_model(),
        "platform": platform.platform(),
        "python": platform.python_version(),
        "cpu_info": spharm_cpp.cpu_info(),
        "pyspharm_version": package_version("pyspharm"),
        "spharm_module_version": getattr(pyspharm, "__version__", "unknown"),
        "ducc0_version": ducc0.__version__,
        "pyspharm_cpp_version": spharm_cpp.__version__,
        "numpy": np.__version__,
        "warmups": args.warmups,
        "repeats": args.repeats,
    }
    print("BENCHMARK_METADATA=" + json.dumps(metadata), flush=True)

    columns = (
        "case", "gridtype", "nlon", "nlat", "ntrunc", "threads",
        "operation", "backend", "median_ms", "p90_ms",
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=columns)
        writer.writeheader()

        for case, gridtype, nlon, nlat, _ in CASES:
            if case not in args.cases:
                continue
            # Use the same recoverable band for scalar and spin-1 transforms.
            lmax = min(nlat - (2 if gridtype == "regular" else 1), (nlon - 1) // 2)
            rng = np.random.default_rng(20261024)
            grid, east, north = (
                rng.standard_normal((nlat, nlon)).astype(np.float32)
                for _ in range(3)
            )
            print(f"BENCHMARK_CASE={case} {nlon}x{nlat} {gridtype} lmax={lmax}", flush=True)

            candidates = (
                ("pyspharm", pyspharm.Spharmt(nlon, nlat, gridtype=gridtype),
                 grid, east, north),
                ("pyspharm-cpp",
                 spharm_cpp.Spharmt(nlon, nlat, gridtype=gridtype, dtype=np.float32),
                 grid, east, north),
                ("pyspharm-cpp-float64",
                 spharm_cpp.Spharmt(nlon, nlat, gridtype=gridtype, dtype=np.float64),
                 grid.astype(np.float64), east.astype(np.float64),
                 north.astype(np.float64)),
            )
            measurements = [
                (backend, pyspharm_operations(sht, scalar, u, v, lmax))
                for backend, sht, scalar, u, v in candidates
            ]
            measurements.append(
                ("ducc0", ducc_operations(
                    np, ducc0, grid, east, north, gridtype, lmax, threads
                ))
            )

            for backend, operations in measurements:
                for operation in OPERATIONS:
                    median, p90, samples = measure(
                        operations[operation], args.warmups, args.repeats
                    )
                    if args.print_samples:
                        print(
                            "BENCHMARK_SAMPLES="
                            + json.dumps({
                                "case": case,
                                "operation": operation,
                                "backend": backend,
                                "samples_ms": samples,
                            }),
                            flush=True,
                        )
                    row = {
                        "case": case, "gridtype": gridtype, "nlon": nlon,
                        "nlat": nlat, "ntrunc": lmax, "threads": threads,
                        "operation": operation, "backend": backend,
                        "median_ms": f"{median:.6f}", "p90_ms": f"{p90:.6f}",
                    }
                    writer.writerow(row)
                    stream.flush()
                    print("BENCHMARK_ROW=" + json.dumps(row), flush=True)


if __name__ == "__main__":
    main()
