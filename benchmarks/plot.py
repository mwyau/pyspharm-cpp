# SPDX-License-Identifier: GPL-2.0-or-later
"""Render five comparable raw-latency SVG charts from benchmark.py CSVs.

Each operation has its own labeled linear millisecond axis, so the much slower
reference never forces candidate times into invisible slivers. Bars within an
operation are directly comparable; bar lengths across operations are not.
"""

from __future__ import annotations

import argparse
import csv
import math
from html import escape
from pathlib import Path

CASES = (
    ("cc-2p5", "Clenshaw–Curtis · 2.5° · 144 × 73"),
    ("cc-2p5-1t", "Clenshaw–Curtis · 2.5° · 144 × 73"),
    ("cc-0p25", "Clenshaw–Curtis · 0.25° · 1440 × 721"),
    ("gl-128", "Gaussian · 128 latitude rings · 256 × 128"),
    ("gl-640", "Gaussian · 640 latitude rings · 1280 × 640"),
)
OPERATIONS = (
    ("scalar_analysis", "Scalar analysis"),
    ("scalar_synthesis", "Scalar synthesis"),
    ("wind_analysis", "Wind analysis"),
)
BACKENDS = (
    ("PySpharm", "PySpharm · f32", "#a25561"),
    ("pyspharm-cpp", "pyspharm-cpp · f32", "#286a9b"),
    ("ducc0", "DUCC · f32", "#3d8462"),
    ("pyspharm-cpp-float64", "pyspharm-cpp · f64", "#7268a8"),
)


def load_measurements(paths: list[Path]):
    results = {}
    for path in paths:
        with path.open(newline="", encoding="utf-8") as fp:
            for row in csv.DictReader(fp):
                key = (row["case"], row["operation"], row["backend"])
                if key in results:
                    raise ValueError(f"Duplicate measurement {key}")
                v = float(row["median_ms"])
                if not math.isfinite(v) or v <= 0:
                    raise ValueError(f"Invalid median for {key}")
                results[key] = v
    for case, _ in CASES:
        for op, _ in OPERATIONS:
            for backend, _, _ in BACKENDS:
                if (case, op, backend) not in results:
                    raise ValueError(f"Missing measurement: {case}, {op}, {backend}")
    return results


def render(case, label, results):
    left, width = 227, 486
    doc = [
        '<svg xmlns="http://www.w3.org/2000/svg" width="920" height="740" '
        'viewBox="0 0 920 740" role="img" aria-labelledby="t d">',
        '<title id="t">Measured spherical harmonic transform timings</title>',
        '<desc id="d">Three operations comparing PySpharm and DUCC float32 '
        'against pyspharm-cpp float32 and float64 on separate zero-based '
        'linear millisecond axes.</desc>',
        '<rect width="920" height="740" fill="white"/>',
        '<text x="27" y="40" font-size="21" font-family="sans-serif" '
        'font-weight="600" fill="#1f3546">Spherical harmonic transform latency</text>',
        f'<text x="27" y="66" font-size="14" font-family="sans-serif" '
        f'fill="#526478">{escape(label)} · PySpharm serial · DUCC {1 if case.endswith("-1t") else 4} thread{"s" if not case.endswith("-1t") else ""}</text>',
    ]
    for i, (operation, title) in enumerate(OPERATIONS):
        y0 = 98 + 202 * i
        values = [results[case, operation, name] for name, _, _ in BACKENDS]
        maximum = max(values) * 1.13
        maximum = float(f"{maximum:.2g}") if maximum >= 1 else math.ceil(maximum * 10) / 10
        if maximum < max(values):
            maximum = max(values) * 1.08
        doc.append(
            f'<text x="27" y="{y0}" font-size="16" font-family="sans-serif" '
            f'font-weight="600" fill="#1f3546">{escape(title)}</text>'
        )
        for j, (name, caption, color) in enumerate(BACKENDS):
            y = y0 + 23 + 29 * j
            val = results[case, operation, name]
            w = val / maximum * width
            doc.extend([
                f'<text x="27" y="{y+17}" font-family="sans-serif" font-size="14" '
                f'fill="#233646">{escape(caption)}</text>',
                f'<rect x="{left}" y="{y}" width="{w:.2f}" height="23" '
                f'fill="{color}" rx="2"/>',
                f'<text x="{left+w+9:.1f}" y="{y+17}" font-family="sans-serif" '
                f'font-size="13" fill="#233646">{val:.3f} ms</text>',
            ])
        axis_y = y0 + 156
        doc.append(f'<line x1="{left}" x2="{left+width}" y1="{axis_y}" '
                   f'y2="{axis_y}" stroke="#9caeba"/>')
        for frac in (0, 0.5, 1):
            x = left + width * frac
            doc.append(f'<text x="{x:.1f}" y="{axis_y+18}" '
                       f'text-anchor="middle" font-size="11" '
                       f'font-family="sans-serif" fill="#536779">'
                       f'{maximum*frac:.3g}</text>')
    doc.append('<text x="27" y="723" font-size="12" '
               'font-family="sans-serif" fill="#526478">'
               'Linear milliseconds; separate axis per operation. '
               'Read printed numbers to compare across operations.</text>')
    doc.append("</svg>")
    return "\n".join(doc) + "\n"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--csv", type=Path, nargs="+", required=True,
        help="Both benchmark CSVs from the same runner",
    )
    parser.add_argument("--output-dir", type=Path, default=Path(__file__).parent)
    args = parser.parse_args()
    results = load_measurements(args.csv)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    for case, label in CASES:
        target = args.output_dir / f"timings-{case}.svg"
        target.write_text(render(case, label, results), encoding="utf-8")
        print(target)


if __name__ == "__main__":
    main()
