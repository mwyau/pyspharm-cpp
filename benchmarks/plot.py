# SPDX-License-Identifier: GPL-2.0-or-later
"""Render four comparable raw-latency SVG charts from benchmark.py CSV.

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
    ("PySpharm", "PySpharm", "#a25561"),
    ("pyspharm-cpp", "pyspharm-cpp", "#286a9b"),
    ("ducc0", "DUCC", "#3d8462"),
)


# Exact median observations from hosted GitHub Actions run 40101800670;
# chart regeneration does not depend on a checked-in raw CSV.
SNAPSHOT = {
    ("cc-2p5", "scalar_analysis", "PySpharm"): 0.069194,
    ("cc-2p5", "scalar_synthesis", "PySpharm"): 0.029174,
    ("cc-2p5", "wind_analysis", "PySpharm"): 0.177977,
    ("cc-2p5", "scalar_analysis", "pyspharm-cpp"): 0.208873,
    ("cc-2p5", "scalar_synthesis", "pyspharm-cpp"): 0.067051,
    ("cc-2p5", "wind_analysis", "pyspharm-cpp"): 0.385208,
    ("cc-2p5", "scalar_analysis", "ducc0"): 0.222424,
    ("cc-2p5", "scalar_synthesis", "ducc0"): 0.063245,
    ("cc-2p5", "wind_analysis", "ducc0"): 0.381092,
    ("cc-0p25", "scalar_analysis", "PySpharm"): 95.026696,
    ("cc-0p25", "scalar_synthesis", "PySpharm"): 76.560638,
    ("cc-0p25", "wind_analysis", "PySpharm"): 244.674017,
    ("cc-0p25", "scalar_analysis", "pyspharm-cpp"): 12.637532,
    ("cc-0p25", "scalar_synthesis", "pyspharm-cpp"): 4.047290,
    ("cc-0p25", "wind_analysis", "pyspharm-cpp"): 27.711609,
    ("cc-0p25", "scalar_analysis", "ducc0"): 14.259676,
    ("cc-0p25", "scalar_synthesis", "ducc0"): 3.875081,
    ("cc-0p25", "wind_analysis", "ducc0"): 29.121713,
    ("gl-128", "scalar_analysis", "PySpharm"): 0.633740,
    ("gl-128", "scalar_synthesis", "PySpharm"): 0.108904,
    ("gl-128", "wind_analysis", "PySpharm"): 1.160861,
    ("gl-128", "scalar_analysis", "pyspharm-cpp"): 0.131477,
    ("gl-128", "scalar_synthesis", "pyspharm-cpp"): 0.127902,
    ("gl-128", "wind_analysis", "pyspharm-cpp"): 0.264737,
    ("gl-128", "scalar_analysis", "ducc0"): 0.127501,
    ("gl-128", "scalar_synthesis", "ducc0"): 0.117677,
    ("gl-128", "wind_analysis", "ducc0"): 0.227181,
    ("gl-640", "scalar_analysis", "PySpharm"): 184.648172,
    ("gl-640", "scalar_synthesis", "PySpharm"): 47.849009,
    ("gl-640", "wind_analysis", "PySpharm"): 612.969222,
    ("gl-640", "scalar_analysis", "pyspharm-cpp"): 3.449693,
    ("gl-640", "scalar_synthesis", "pyspharm-cpp"): 3.064896,
    ("gl-640", "wind_analysis", "pyspharm-cpp"): 7.590132,
    ("gl-640", "scalar_analysis", "ducc0"): 3.199028,
    ("gl-640", "scalar_synthesis", "ducc0"): 2.963785,
    ("gl-640", "wind_analysis", "ducc0"): 6.404595,
}

def load_measurements(path: Path):
    results = {}
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
        '<svg xmlns="http://www.w3.org/2000/svg" width="920" height="702" '
        'viewBox="0 0 920 702" role="img" aria-labelledby="t d">',
        '<title id="t">Measured spherical harmonic transform timings</title>',
        '<desc id="d">Three operations, each comparing PySpharm, pyspharm-cpp '
        'and DUCC on a zero-based linear millisecond scale. Each operation '
        'uses a different labeled axis extent.</desc>',
        '<rect width="920" height="702" fill="white"/>',
        '<text x="27" y="40" font-size="21" font-family="sans-serif" '
        'font-weight="600" fill="#1f3546">Spherical harmonic transform latency</text>',
        f'<text x="27" y="66" font-size="14" font-family="sans-serif" '
        f'fill="#526478">{escape(label)} · float32 · PySpharm serial · DUCC 4 threads</text>',
    ]
    for i, (operation, title) in enumerate(OPERATIONS):
        y0 = 98 + 190 * i
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
            y = y0 + 23 + 34 * j
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
        axis_y = y0 + 138
        doc.append(f'<line x1="{left}" x2="{left+width}" y1="{axis_y}" '
                   f'y2="{axis_y}" stroke="#9caeba"/>')
        for frac in (0, 0.5, 1):
            x = left + width * frac
            doc.append(f'<text x="{x:.1f}" y="{axis_y+18}" '
                       f'text-anchor="middle" font-size="11" '
                       f'font-family="sans-serif" fill="#536779">'
                       f'{maximum*frac:.3g}</text>')
    doc.append('<text x="27" y="685" font-size="12" '
               'font-family="sans-serif" fill="#526478">'
               'Linear milliseconds; separate axis per operation. '
               'Read printed numbers to compare across operations.</text>')
    doc.append("</svg>")
    return "\n".join(doc) + "\n"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--csv", type=Path)
    parser.add_argument("--output-dir", type=Path, default=Path(__file__).parent)
    args = parser.parse_args()
    results = load_measurements(args.csv) if args.csv else SNAPSHOT
    args.output_dir.mkdir(parents=True, exist_ok=True)
    for case, label in CASES:
        target = args.output_dir / f"timings-{case}.svg"
        target.write_text(render(case, label, results), encoding="utf-8")
        print(target)


if __name__ == "__main__":
    main()
