# SPDX-License-Identifier: GPL-2.0-or-later
"""Check the wheel ABI and platform matrix produced by cibuildwheel."""

import os
from pathlib import Path

wheels = sorted(Path("wheelhouse").rglob("*.whl"))
target = os.environ["TARGET"]
arch = os.environ["EXPECTED_WHEEL_ARCH"]
deployment = os.environ.get("MACOS_DEPLOYMENT_TARGET", "")
expected_pairs = {
    ("cp311", "cp311"),
    ("cp312", "abi3"),
    ("cp314", "cp314t"),
    ("cp315", "cp315t"),
}
if target.startswith("linux-"):
    families = {"manylinux", "musllinux"}
elif target.startswith("macos-"):
    families = {"macos"}
else:
    families = {"windows"}
expected = {(python, abi, family) for python, abi in expected_pairs for family in families}

actual = set()
for wheel in wheels:
    _, python, abi, platform = wheel.stem.rsplit("-", 3)
    tags = platform.split(".")
    if target.startswith("linux-"):
        for tag in tags:
            assert tag.endswith(f"_{arch}"), wheel.name
        family_set = {
            "manylinux" if tag.startswith("manylinux") else "musllinux"
            for tag in tags
        }
        assert all(tag.startswith(("manylinux", "musllinux")) for tag in tags)
    elif target.startswith("macos-"):
        prefix = f"macosx_{deployment.replace('.', '_')}_"
        assert any(
            tag.startswith(prefix) and tag.endswith(f"_{arch}") for tag in tags
        ), wheel.name
        family_set = {"macos"}
    else:
        assert any(tag == f"win_{arch}" for tag in tags), wheel.name
        family_set = {"windows"}
    assert len(family_set) == 1, wheel.name
    actual.add((python, abi, next(iter(family_set))))

assert len(wheels) == len(expected), (
    f"Expected {len(expected)} wheels, got {len(wheels)}"
)
assert actual == expected, (
    f"Wheel matrix mismatch: missing={sorted(expected - actual)}, "
    f"unexpected={sorted(actual - expected)}"
)
print(f"Validated {len(wheels)} wheels for {target}")
