# SPDX-FileCopyrightText: 2026 Albert M. W. Yau
# SPDX-License-Identifier: GPL-2.0-or-later
"""Check the ABI and platform tags of a single cibuildwheel distribution."""

import os
from pathlib import Path

wheels = list(Path("wheelhouse").glob("*.whl"))
assert len(wheels) == 1, f"Expected one wheel, found {len(wheels)}"

build = os.environ["CIBW_BUILD"].split("-", 1)[0]
expected = {
    "cp311": ("cp311", "cp311"),
    "cp312": ("cp312", "abi3"),
    "cp314t": ("cp314", "cp314t"),
    "cp315t": ("cp315", "cp315t"),
}[build]
family = os.environ["WHEEL_FAMILY"]
arch = os.environ["EXPECTED_WHEEL_ARCH"]
deployment = os.environ.get("MACOS_DEPLOYMENT_TARGET", "")
wheel = wheels[0]
_, python, abi, platform = wheel.stem.rsplit("-", 3)
assert (python, abi) == expected, wheel.name
tags = platform.split(".")

if family in {"manylinux", "musllinux"}:
    assert all(
        tag.startswith(f"{family}_") and tag.endswith(f"_{arch}") for tag in tags
    ), wheel.name
elif family == "macos":
    prefix = f"macosx_{deployment.replace('.', '_')}_"
    assert any(tag.startswith(prefix) and tag.endswith(f"_{arch}") for tag in tags), (
        wheel.name
    )
    assert all(tag.endswith(f"_{arch}") for tag in tags), wheel.name
else:
    assert family == "windows"
    assert tags == [f"win_{arch}"], wheel.name

print(f"Validated {wheel.name} ({family})")
