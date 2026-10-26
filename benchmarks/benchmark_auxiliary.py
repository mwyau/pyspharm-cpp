# SPDX-FileCopyrightText: 2026 Albert M. W. Yau
# SPDX-License-Identifier: GPL-2.0-or-later
"""Measure scalar-call and batched Legendre and interpolation kernels."""

from __future__ import annotations

import argparse
import csv
import gc
import hashlib
import json
import math
import os
import platform
import statistics
import subprocess
import time
from pathlib import Path

CASES = (
    ("small-low", 8, 16),
    ("medium-medium", 64, 64),
    ("large-medium", 128, 256),
    ("small-high", 16, 1024),
)
OPERATIONS = (
    "legendre_scalar",
    "legendre_batched",
    "specintrp_scalar",
    "specintrp_batched",
)
MAX_WORKING_SET = 256 * 1024**2


def source_sha() -> str:
    configured = os.environ.get("SPHARM_CPP_BENCHMARK_SOURCE_SHA")
    if configured:
        return configured
    repository = Path(__file__).resolve().parents[1]
    try:
        revision = subprocess.check_output(
            ["git", "rev-parse", "HEAD"],
            cwd=repository,
            text=True,
            stderr=subprocess.DEVNULL,
        ).strip()
        diff = subprocess.check_output(
            ["git", "diff", "--binary", "HEAD"],
            cwd=repository,
            stderr=subprocess.DEVNULL,
        )
    except (OSError, subprocess.CalledProcessError):
        return "unknown"
    if diff:
        return f"{revision}+dirty:{hashlib.sha256(diff).hexdigest()[:12]}"
    return revision


def cpu_model() -> str:
    cpuinfo = Path("/proc/cpuinfo")
    if cpuinfo.is_file():
        for line in cpuinfo.read_text(encoding="utf-8").splitlines():
            if line.lower().startswith("model name"):
                return line.split(":", 1)[1].strip()
    return platform.processor() or "unknown"


def cpu_quota_cores() -> float | None:
    cpu_max = Path("/sys/fs/cgroup/cpu.max")
    if not cpu_max.is_file():
        return None
    quota, period = cpu_max.read_text(encoding="utf-8").split()
    if quota == "max":
        return None
    return int(quota) / int(period)


def cpu_affinity() -> list[int] | None:
    if not hasattr(os, "sched_getaffinity"):
        return None
    return sorted(os.sched_getaffinity(0))


def compiler_version() -> str:
    compiler = os.environ.get("CXX", "c++")
    try:
        result = subprocess.run(
            [compiler, "--version"], capture_output=True, text=True, check=True
        )
    except (OSError, subprocess.CalledProcessError):
        return "unavailable"
    return result.stdout.splitlines()[0] if result.stdout else "unknown"


def measure(callback, warmups: int, repeats: int) -> tuple[float, float]:
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


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--threads", type=int, required=True, choices=(1, 2, 4))
    parser.add_argument("--warmups", type=int, default=2)
    parser.add_argument("--repeats", type=int, default=7)
    parser.add_argument("--case", choices=tuple(case[0] for case in CASES))
    parser.add_argument("--operation", choices=OPERATIONS)
    args = parser.parse_args()
    if args.warmups < 1 or args.repeats < 3:
        parser.error("warmups >= 1; repeats >= 3")

    os.environ["DUCC0_NUM_THREADS"] = str(args.threads)
    import numpy as np
    import spharm_cpp

    metadata = {
        "source_sha": source_sha(),
        "compiler": os.environ.get("CXX", "c++"),
        "compiler_version": compiler_version(),
        "cpu_model": cpu_model(),
        "cpu_affinity": cpu_affinity(),
        "cpu_quota_cores": cpu_quota_cores(),
        "threads": args.threads,
        "case": args.case or "all",
        "operation": args.operation or "all",
        "platform": platform.platform(),
        "python": platform.python_version(),
        "cpu_info": spharm_cpp.cpu_info(),
        "numpy": np.__version__,
        "warmups": args.warmups,
        "repeats": args.repeats,
        "max_working_set_bytes": MAX_WORKING_SET,
    }
    print("BENCHMARK_METADATA=" + json.dumps(metadata), flush=True)

    columns = (
        "case",
        "points",
        "lmax",
        "dtype",
        "threads",
        "operation",
        "median_ms",
        "p90_ms",
        "working_set_bytes",
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=columns)
        writer.writeheader()
        for case, points, lmax in CASES:
            if args.case is not None and case != args.case:
                continue
            latitudes = np.linspace(-90.0, 90.0, points, dtype=np.float64)
            longitudes = np.linspace(-720.0, 720.0, points, dtype=np.float64)
            ncoeff = (lmax + 1) * (lmax + 2) // 2
            for dtype, cdtype in (
                (np.float32, np.complex64),
                (np.float64, np.complex128),
            ):
                legendre_bytes = points * ncoeff * np.dtype(dtype).itemsize
                spectrum_bytes = ncoeff * np.dtype(cdtype).itemsize
                result_bytes = points * np.dtype(dtype).itemsize
                working_set_bytes = 2 * legendre_bytes + spectrum_bytes + result_bytes
                if working_set_bytes > MAX_WORKING_SET:
                    raise RuntimeError(
                        f"{case} {dtype.__name__} exceeds the working-set limit"
                    )
                legendre = spharm_cpp.legendre(latitudes, lmax, dtype=dtype)
                spectrum = (
                    np.random.default_rng(2468).standard_normal(ncoeff)
                    + 1j * np.random.default_rng(9753).standard_normal(ncoeff)
                ).astype(cdtype)
                spectrum[0] = spectrum[0].real
                callbacks = {
                    "legendre_scalar": lambda: [
                        spharm_cpp.legendre(float(latitude), lmax, dtype=dtype)
                        for latitude in latitudes
                    ],
                    "legendre_batched": lambda: spharm_cpp.legendre(
                        latitudes, lmax, dtype=dtype
                    ),
                    "specintrp_scalar": lambda: [
                        spharm_cpp.specintrp(float(longitude), spectrum, row)
                        for longitude, row in zip(longitudes, legendre)
                    ],
                    "specintrp_batched": lambda: spharm_cpp.specintrp(
                        longitudes, spectrum, legendre
                    ),
                }
                for operation in OPERATIONS:
                    if args.operation is not None and operation != args.operation:
                        continue
                    median, p90 = measure(
                        callbacks[operation], args.warmups, args.repeats
                    )
                    row = {
                        "case": case,
                        "points": points,
                        "lmax": lmax,
                        "dtype": np.dtype(dtype).name,
                        "threads": args.threads,
                        "operation": operation,
                        "median_ms": f"{median:.6f}",
                        "p90_ms": f"{p90:.6f}",
                        "working_set_bytes": working_set_bytes,
                    }
                    writer.writerow(row)
                    stream.flush()
                    print("BENCHMARK_ROW=" + json.dumps(row), flush=True)


if __name__ == "__main__":
    main()
