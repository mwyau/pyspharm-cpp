# SPDX-License-Identifier: GPL-2.0-or-later
"""The checked-in benchmark figures must match their reproducible source."""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path


def test_benchmark_figures_match_recorded_medians(tmp_path: Path) -> None:
    directory = Path(__file__).resolve().parents[1] / "benchmarks"
    subprocess.run(
        [
            sys.executable,
            str(directory / "plot.py"),
            "--output-dir",
            str(tmp_path),
        ],
        check=True,
        capture_output=True,
        text=True,
    )
    for case in ("cc-2p5", "cc-0p25", "gl-128", "gl-640"):
        filename = f"timings-{case}.svg"
        assert (tmp_path / filename).read_text(encoding="utf-8") == (
            directory / filename
        ).read_text(encoding="utf-8")
