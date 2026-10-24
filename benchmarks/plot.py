# SPDX-License-Identifier: GPL-2.0-or-later
"""Draw five SHT latency charts from the paired benchmark CSV files."""

from __future__ import annotations

import argparse
import csv
import math
from html import escape
from pathlib import Path

from benchmark import CASES, OPERATIONS

TITLES = {
    "scalar_analysis": "Scalar analysis",
    "scalar_synthesis": "Scalar synthesis",
    "wind_analysis": "Wind analysis",
    "wind_synthesis": "Wind synthesis",
}
BACKENDS = (
    ("PySpharm", "pyspharm · f32", "#a25561"),
    ("pyspharm-cpp", "pyspharm-cpp · f32", "#286a9b"),
    ("ducc0", "ducc0 · f32", "#3d8462"),
    ("pyspharm-cpp-float64", "pyspharm-cpp · f64", "#7268a8"),
)
HEIGHT = 136 + 202 * len(OPERATIONS)


def load_measurements(paths: list[Path]):
    results = {}
    for path in paths:
        with path.open(newline="", encoding="utf-8") as stream:
            for row in csv.DictReader(stream):
                key = (row["case"], row["operation"], row["backend"])
                expected_threads = 1 if row["case"] == "cc-2p5-1t" else 4
                if int(row["threads"]) != expected_threads:
                    raise ValueError(
                        f"Incorrect thread count for {key}: "
                        f"expected {expected_threads}, got {row['threads']}"
                    )
                if key in results:
                    raise ValueError(f"Duplicate measurement {key}")
                median = float(row["median_ms"])
                if not math.isfinite(median) or median <= 0:
                    raise ValueError(f"Invalid median for {key}")
                results[key] = median

    expected = {
        (case, operation, backend)
        for case, *_ in CASES
        for operation in OPERATIONS
        for backend, *_ in BACKENDS
    }
    if results.keys() != expected:
        raise ValueError(
            f"Missing measurements: {sorted(expected - results.keys())}; "
            f"unexpected measurements: {sorted(results.keys() - expected)}"
        )
    return results


def render(case, label, measurements):
    left, width = 227, 486
    threads = 1 if case.endswith("-1t") else 4
    lines = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="920" height="{HEIGHT}" '
        f'viewBox="0 0 920 {HEIGHT}" role="img" aria-labelledby="t d">',
        '<title id="t">Measured spherical harmonic transform timings</title>',
        '<desc id="d">Scalar and wind analysis and synthesis latency; four backends '
        'on separate zero-based linear millisecond axes.</desc>',
        f'<rect width="920" height="{HEIGHT}" fill="white"/>',
        '<text x="27" y="40" font-size="21" font-family="sans-serif" '
        'font-weight="600" fill="#1f3546">Spherical harmonic transform latency</text>',
        f'<text x="27" y="66" font-size="14" font-family="sans-serif" '
        f'fill="#526478">{escape(label)} · pyspharm serial · ducc0/pyspharm-cpp '
        f'{threads} thread{"s" if threads != 1 else ""}</text>',
    ]

    for i, operation in enumerate(OPERATIONS):
        y0 = 98 + 202 * i
        values = [measurements[case, operation, name] for name, *_ in BACKENDS]
        maximum = max(values) * 1.13
        maximum = float(f"{maximum:.2g}") if maximum >= 1 else math.ceil(maximum * 10) / 10
        if maximum < max(values):
            maximum = max(values) * 1.08

        lines.append(
            f'<text x="27" y="{y0}" font-size="16" font-family="sans-serif" '
            f'font-weight="600" fill="#1f3546">{escape(TITLES[operation])}</text>'
        )
        for j, (name, caption, color) in enumerate(BACKENDS):
            y = y0 + 23 + 29 * j
            value = measurements[case, operation, name]
            bar_width = value / maximum * width
            lines.extend((
                f'<text x="27" y="{y + 17}" font-family="sans-serif" font-size="14" '
                f'fill="#233646">{escape(caption)}</text>',
                f'<rect x="{left}" y="{y}" width="{bar_width:.2f}" '
                f'height="23" fill="{color}" rx="2"/>',
                f'<text x="{left + bar_width + 9:.1f}" y="{y + 17}" '
                f'font-family="sans-serif" font-size="13" fill="#233646">'
                f'{value:.3f} ms</text>',
            ))

        axis_y = y0 + 156
        lines.append(
            f'<line x1="{left}" x2="{left + width}" y1="{axis_y}" '
            f'y2="{axis_y}" stroke="#9caeba"/>'
        )
        for fraction in (0, 0.5, 1):
            x = left + width * fraction
            lines.append(
                f'<text x="{x:.1f}" y="{axis_y + 18}" text-anchor="middle" '
                f'font-size="11" font-family="sans-serif" fill="#536779">'
                f'{maximum * fraction:.3g}</text>'
            )

    lines.extend((
        f'<text x="27" y="{HEIGHT - 17}" font-size="12" '
        'font-family="sans-serif" fill="#526478">'
        'Linear milliseconds; separate axis per operation. '
        'Read printed numbers to compare across operations.</text>',
        "</svg>",
    ))
    return "\n".join(lines) + "\n"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--csv", type=Path, nargs="+", required=True)
    parser.add_argument("--output-dir", type=Path, default=Path(__file__).parent)
    args = parser.parse_args()
    measurements = load_measurements(args.csv)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    for case, _, _, _, label in CASES:
        output = args.output_dir / f"timings-{case}.svg"
        output.write_text(render(case, label, measurements), encoding="utf-8")
        print(output)


if __name__ == "__main__":
    main()
