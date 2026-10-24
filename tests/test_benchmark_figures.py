# SPDX-License-Identifier: GPL-2.0-or-later
"""Benchmark plots require five cases and four measured series per operation."""

from __future__ import annotations

import csv
import subprocess
import sys
from pathlib import Path


def test_benchmark_figures_render_from_measurements(tmp_path: Path) -> None:
    directory = Path(__file__).resolve().parents[1] / "benchmarks"
    cases = ("cc-2p5", "cc-2p5-1t", "cc-0p25", "gl-128", "gl-640")
    operations = ("scalar_analysis", "scalar_synthesis", "wind_analysis")
    backends = (
        "PySpharm",
        "pyspharm-cpp",
        "ducc0",
        "pyspharm-cpp-float64",
    )
    csv_path = tmp_path / "measurements.csv"
    with csv_path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(
            stream, fieldnames=("case", "operation", "backend", "threads", "median_ms")
        )
        writer.writeheader()
        for case in cases:
            for operation in operations:
                for index, backend in enumerate(backends, start=1):
                    writer.writerow(
                        {
                            "case": case,
                            "operation": operation,
                            "backend": backend,
                            "threads": 1 if case.endswith("-1t") else 4,
                            "median_ms": f"{index * 0.125:.6f}",
                        }
                    )

    subprocess.run(
        [
            sys.executable,
            str(directory / "plot.py"),
            "--csv",
            str(csv_path),
            "--output-dir",
            str(tmp_path),
        ],
        check=True,
        capture_output=True,
        text=True,
    )
    for case in cases:
        figure = (tmp_path / f"timings-{case}.svg").read_text(encoding="utf-8")
        assert figure.startswith('<svg xmlns="http://www.w3.org/2000/svg"')
        assert figure.count("<rect ") == 13  # Background + 3 operations × 4 bars
        for label in ("PySpharm", "pyspharm-cpp", "DUCC", "f64"):
            assert label in figure
