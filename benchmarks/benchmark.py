# SPDX-License-Identifier: GPL-2.0-or-later
"""Compare float32 and float64 SHT latencies on five grid cases."""

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
from pathlib import Path

CASES = (
    ("cc-2p5", "regular", 144, 73),
    ("cc-2p5-1t", "regular", 144, 73),
    ("cc-0p25", "regular", 1440, 721),
    ("gl-128", "gaussian", 256, 128),
    ("gl-640", "gaussian", 1280, 640),
)
OPS = ("scalar_analysis", "scalar_synthesis", "wind_analysis")


def available_threads() -> int:
    cpus = os.cpu_count() or 1
    if hasattr(os, "sched_getaffinity"):
        cpus = min(cpus, len(os.sched_getaffinity(0)))
    try:
        value, period = Path("/sys/fs/cgroup/cpu.max").read_text().split()
        if value != "max":
            cpus = min(cpus, max(1, math.floor(int(value) / int(period))))
    except (OSError, ValueError):
        pass
    return max(1, min(8, cpus))


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
    return statistics.median(samples), samples[math.ceil(0.9 * repeats) - 1]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--threads", type=int, default=0)
    parser.add_argument(
        "--cases", nargs="+",
        choices=tuple(entry[0] for entry in CASES),
        default=["cc-2p5", "cc-0p25", "gl-128", "gl-640"],
        help="grid cases to benchmark",
    )
    parser.add_argument("--warmups", type=int, default=2)
    parser.add_argument("--repeats", type=int, default=7)
    args = parser.parse_args()
    if args.threads < 0 or args.warmups < 1 or args.repeats < 3:
        parser.error("threads >= 0; warmups >= 1; repeats >= 3")
    threads = args.threads or available_threads()
    if "cc-2p5-1t" in args.cases and threads != 1:
        parser.error("cc-2p5-1t requires --threads 1")
    os.environ["DUCC0_NUM_THREADS"] = str(threads)
    os.environ["OMP_NUM_THREADS"] = str(threads)
    os.environ["OPENBLAS_NUM_THREADS"] = "1"
    import numpy as np
    import ducc0
    import spharm
    import spharm_cpp

    metadata = {
        "threads": threads, "cpus": os.cpu_count(),
        "quota_cpus": available_threads(), "platform": platform.platform(),
        "python": platform.python_version(),
        "cpu_info": spharm_cpp.cpu_info(),
        "pyspharm_version": getattr(spharm, "__version__", "unknown"),
        "ducc0_version": ducc0.__version__,
        "pyspharm_cpp_version": spharm_cpp.__version__,
        "numpy": np.__version__, "warmups": args.warmups,
        "repeats": args.repeats,
    }
    print("BENCHMARK_METADATA=" + json.dumps(metadata, default=str), flush=True)
    columns = ("case", "gridtype", "nlon", "nlat", "ntrunc", "threads",
               "operation", "backend", "median_ms", "p90_ms")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", encoding="utf-8", newline="") as fp:
        writer = csv.DictWriter(fp, fieldnames=columns)
        writer.writeheader()
        for case, gridtype, nlon, nlat in CASES:
            if case not in args.cases:
                continue
            # Common mathematically supported band, avoiding vector Nyquist.
            lmax = min(nlat - (2 if gridtype == "regular" else 1),
                       (nlon - 1) // 2)
            rng = np.random.default_rng(20261024)
            grid = rng.standard_normal((nlat, nlon)).astype(np.float32)
            east = rng.standard_normal((nlat, nlon)).astype(np.float32)
            north = rng.standard_normal((nlat, nlon)).astype(np.float32)
            print(f"BENCHMARK_CASE={case} {nlon}x{nlat} "
                  f"{gridtype} lmax={lmax}", flush=True)
            # Match sampled input values when comparing candidate precisions.
            grid64 = grid.astype(np.float64)
            east64 = east.astype(np.float64)
            north64 = north.astype(np.float64)
            implementations = (
                ("PySpharm", spharm.Spharmt(nlon, nlat, gridtype=gridtype),
                 grid, east, north),
                ("pyspharm-cpp", spharm_cpp.Spharmt(
                    nlon, nlat, gridtype=gridtype, dtype=np.float32),
                 grid, east, north),
                ("pyspharm-cpp-float64", spharm_cpp.Spharmt(
                    nlon, nlat, gridtype=gridtype, dtype=np.float64),
                 grid64, east64, north64),
            )
            for backend, sht, scalar, zonal, meridional in implementations:
                coeff = sht.grdtospec(scalar, lmax)
                operations = {
                    "scalar_analysis": lambda: sht.grdtospec(scalar, lmax),
                    "scalar_synthesis": lambda: sht.spectogrd(coeff),
                    "wind_analysis": lambda: sht.getvrtdivspec(
                        zonal, meridional, lmax),
                }
                for operation in OPS:
                    median, p90 = measure(operations[operation],
                                          args.warmups, args.repeats)
                    row = dict(case=case, gridtype=gridtype, nlon=nlon,
                               nlat=nlat, ntrunc=lmax, threads=threads,
                               operation=operation, backend=backend,
                               median_ms=f"{median:.6f}", p90_ms=f"{p90:.6f}")
                    writer.writerow(row)
                    fp.flush()
                    print("BENCHMARK_ROW=" + json.dumps(row), flush=True)
            geometry = "CC" if gridtype == "regular" else "GL"
            alm = ducc0.sht.analysis_2d(
                map=grid[None, :, :], spin=0, lmax=lmax,
                geometry=geometry, nthreads=threads)
            winds = np.stack((east, north))
            operations = {
                "scalar_analysis": lambda: ducc0.sht.analysis_2d(
                    map=grid[None, :, :], spin=0, lmax=lmax,
                    geometry=geometry, nthreads=threads),
                "scalar_synthesis": lambda: ducc0.sht.synthesis_2d(
                    alm=alm, spin=0, lmax=lmax, geometry=geometry,
                    ntheta=nlat, nphi=nlon, nthreads=threads),
                "wind_analysis": lambda: ducc0.sht.analysis_2d(
                    map=winds, spin=1, lmax=lmax,
                    geometry=geometry, nthreads=threads),
            }
            for operation in OPS:
                median, p90 = measure(operations[operation],
                                      args.warmups, args.repeats)
                row = dict(case=case, gridtype=gridtype, nlon=nlon,
                           nlat=nlat, ntrunc=lmax, threads=threads,
                           operation=operation, backend="ducc0",
                           median_ms=f"{median:.6f}", p90_ms=f"{p90:.6f}")
                writer.writerow(row)
                fp.flush()
                print("BENCHMARK_ROW=" + json.dumps(row), flush=True)


if __name__ == "__main__":
    main()
