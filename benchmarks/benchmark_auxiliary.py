# SPDX-License-Identifier: GPL-2.0-or-later
"""Measure batched Legendre and spectral interpolation kernels."""

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
    ("small-low", 8, 16),
    ("medium-medium", 64, 64),
    ("large-medium", 128, 256),
    ("small-high", 16, 1024),
)
OPERATIONS = ("legendre_many", "specintrp_many")
MAX_WORKING_SET = 256 * 1024**2


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
    args = parser.parse_args()
    if args.warmups < 1 or args.repeats < 3:
        parser.error("warmups >= 1; repeats >= 3")

    os.environ["DUCC0_NUM_THREADS"] = str(args.threads)
    import numpy as np
    import spharm_cpp

    metadata = {
        "threads": args.threads,
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
                legendre = spharm_cpp.legendre_many(latitudes, lmax, dtype=dtype)
                spectrum = (
                    np.random.default_rng(2468).standard_normal(ncoeff)
                    + 1j * np.random.default_rng(9753).standard_normal(ncoeff)
                ).astype(cdtype)
                spectrum[0] = spectrum[0].real
                callbacks = {
                    "legendre_many": lambda: spharm_cpp.legendre_many(
                        latitudes, lmax, dtype=dtype
                    ),
                    "specintrp_many": lambda: spharm_cpp.specintrp_many(
                        longitudes, spectrum, legendre
                    ),
                }
                for operation in OPERATIONS:
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
