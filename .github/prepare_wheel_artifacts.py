# SPDX-License-Identifier: GPL-2.0-or-later
"""Name downloadable test wheels by their compilation profile."""

import argparse
from pathlib import Path
from shutil import copy2

parser = argparse.ArgumentParser()
parser.add_argument("source", type=Path)
parser.add_argument("output", type=Path)
parser.add_argument("profile")
parser.add_argument("--expected", type=int)
args = parser.parse_args()

wheels = sorted(args.source.rglob("*.whl"))
assert wheels, f"No wheels in {args.source}"
if args.expected is not None:
    assert len(wheels) == args.expected, f"Expected {args.expected} wheels"
args.output.mkdir(parents=True, exist_ok=True)
for wheel in wheels:
    distribution, python, abi, platform = wheel.stem.rsplit("-", 3)
    name = f"{distribution}-0_{args.profile}-{python}-{abi}-{platform}.whl"
    copy2(wheel, args.output / name)
print(f"Prepared {len(wheels)} {args.profile} wheel artifacts")
