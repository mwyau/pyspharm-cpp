# SPDX-License-Identifier: GPL-2.0-or-later
"""Plot three measured SHT speedups from the benchmark CSV, using only stdlib.

With no --csv argument, reproduce the checked-in historical figure from the
summary medians below. These are observations on one host, not performance
guarantees. Supplying --csv replaces the snapshots with newly measured values.
"""

from __future__ import annotations

import argparse
import csv
from html import escape
from pathlib import Path

# Observed float32, 512x256 Gaussian, 1 thread; CPython 3.15.0rc2,
# AVX2 host, PySpharm 1.0.9 and the multiarch pyspharm-cpp 0.1.0 wheel.
SNAPSHOT = {
    "scalar_analysis": (1.113594, 28.545004),
    "wind_analysis": (2.863599, 20.102418),
    "multi_field_analysis": (9.432003, 225.173798),
}
OPERATIONS = [
    ("scalar_analysis", "Scalar analysis"),
    ("wind_analysis", "Wind analysis"),
    ("multi_field_analysis", "8-field analysis"),
]


def load_csv(path: Path) -> dict[str, tuple[float, float]]:
    """Read candidate and PySpharm medians from benchmark.py's CSV."""
    values: dict[str, dict[str, float]] = {}
    with path.open(newline="", encoding="utf-8") as stream:
        for row in csv.DictReader(stream):
            if (
                row["nlon"] != "512"
                or row["nlat"] != "256"
                or row["dtype"] != "float32"
                or row["threads"] != "1"
            ):
                continue
            if row["operation"] not in dict(OPERATIONS):
                continue
            if row["backend"] not in ("pyspharm-cpp", "PySpharm"):
                continue
            values.setdefault(row["operation"], {})[row["backend"]] = float(
                row["median_ms"]
            )
    pairs = {}
    for operation, _ in OPERATIONS:
        try:
            a = values[operation]["pyspharm-cpp"]
            b = values[operation]["PySpharm"]
        except KeyError as exc:
            raise ValueError(
                f"missing PySpharm/candidate measurements for {operation}; "
                "install both packages into the benchmark environment"
            ) from exc
        if a <= 0 or b <= 0:
            raise ValueError("median latencies must be positive")
        pairs[operation] = (a, b)
    return pairs


def render_svg(pairs: dict[str, tuple[float, float]]) -> str:
    """Render a deterministic, accessible SVG without plotting dependencies."""
    max_speedup = max(reference / candidate for candidate, reference in pairs.values())
    axis_max = max(10, ((int(max_speedup) // 10) + 1) * 10)
    start, span = 242, 440
    lines = [
        '<svg xmlns="http://www.w3.org/2000/svg" width="860" height="345" '
        'viewBox="0 0 860 345" role="img" '
        'aria-labelledby="chart-title chart-desc">',
        '<title id="chart-title">Observed PySpharm-to-pyspharm-cpp speedup</title>',
        '<desc id="chart-desc">Three operations at 512 by 256, float32, '
        'one thread. Larger ratios mean lower pyspharm-cpp latency.</desc>',
        '<rect x="0" y="0" width="860" height="345" fill="white"/>',
        '<text x="26" y="41" font-family="sans-serif" font-size="20" '
        'font-weight="600" fill="#1e293b">Observed SHT speedup</text>',
        '<text x="26" y="66" font-family="sans-serif" font-size="13" '
        'fill="#475569">PySpharm 1.0.9 latency / pyspharm-cpp latency · '
        '512×256 Gaussian · float32 · one thread</text>',
    ]
    for tick in range(0, axis_max + 1, 10):
        x = start + tick / axis_max * span
        lines.append(
            f'<line x1="{x:.1f}" y1="92" x2="{x:.1f}" y2="289" '
            'stroke="#dce2e8" stroke-width="1"/>'
        )
        lines.append(
            f'<text x="{x:.1f}" y="308" text-anchor="middle" '
            'font-family="sans-serif" font-size="12" fill="#475569">'
            f'{tick}×</text>'
        )
    for i, (key, label) in enumerate(OPERATIONS):
        value, baseline = pairs[key]
        ratio = baseline / value
        y = 110 + i * 70
        width = ratio / axis_max * span
        lines.extend(
            [
                f'<text x="26" y="{y+23}" font-family="sans-serif" '
                f'font-size="14" fill="#1e293b">{escape(label)}</text>',
                f'<rect x="{start}" y="{y}" width="{width:.1f}" height="33" '
                'fill="#35698c" rx="2"/>',
                f'<text x="{start+width+9:.1f}" y="{y+23}" '
                'font-family="sans-serif" font-weight="600" '
                f'font-size="14" fill="#1e293b">{ratio:.1f}×</text>',
            ]
        )
    lines.extend(
        [
            '<text x="26" y="335" font-family="sans-serif" font-size="11" '
            'fill="#64748b">Historical single-host medians; workload-dependent. '
            'Re-run benchmark.py for your system.</text>',
            "</svg>",
        ]
    )
    return "\n".join(lines) + "\n"


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--csv", type=Path, help="benchmark.py CSV containing PySpharm and candidate"
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path(__file__).with_name("speedup.svg"),
    )
    args = parser.parse_args()
    pairs = load_csv(args.csv) if args.csv else SNAPSHOT
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(render_svg(pairs), encoding="utf-8")
    print(f"Wrote {args.output}")


if __name__ == "__main__":
    main()
