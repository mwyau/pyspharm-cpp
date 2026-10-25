# SPDX-FileCopyrightText: 2026 Albert M. W. Yau
# SPDX-License-Identifier: GPL-2.0-or-later
"""Reference loading, status collection, and numerical error recording."""

from __future__ import annotations

import importlib.metadata
import json
import os
import platform
import subprocess
import sys
import sysconfig
from pathlib import Path
from typing import Any

import numpy as np
import pytest

from .manifest import APPROVED_DIFFERENCES, PUBLIC_CONTRACTS, UNAPPROVED_DIFFERENCES

REFERENCE = None
REFERENCE_IMPORT_ERROR: str | None = None
ENVIRONMENT: dict[str, Any] = {}
COMPARISONS: list[dict[str, Any]] = []
INVESTIGATIONS: dict[str, dict[str, Any]] = {}
CONTRACT_STATUSES: dict[str, list[dict[str, Any]]] = {}
DIFFERENCE_STATUSES: dict[str, list[dict[str, Any]]] = {}
UNAPPROVED_DIFFERENCE_STATUSES: dict[str, list[dict[str, Any]]] = {}


def _report_directory() -> Path:
    path = Path(
        os.environ.get("PYSPHARM_VALIDATION_DIR", "/tmp/pyspharm-validation-artifacts")
    )
    path.mkdir(parents=True, exist_ok=True)
    return path


def _package_version(name: str) -> str | None:
    try:
        return importlib.metadata.version(name)
    except importlib.metadata.PackageNotFoundError:
        return None


def _git_worktree_clean(root: Path) -> bool | None:
    try:
        status = subprocess.check_output(
            ["git", "-C", str(root), "status", "--porcelain"], text=True
        )
    except (OSError, subprocess.CalledProcessError):
        return None
    return not status.strip()


def pytest_configure(config: pytest.Config) -> None:
    config.addinivalue_line(
        "markers", "contract(name): public API item covered by this test"
    )
    config.addinivalue_line(
        "markers", "approved_difference(name): tested approved API difference"
    )
    config.addinivalue_line(
        "markers", "unapproved_difference(name): observed difference requiring approval"
    )


def pytest_sessionstart(session: pytest.Session) -> None:
    global REFERENCE, REFERENCE_IMPORT_ERROR
    root = Path(__file__).resolve().parents[2]
    try:
        ENVIRONMENT.update(
            {
                "python_version": sys.version.split()[0],
                "python_implementation": platform.python_implementation(),
                "platform": platform.platform(),
                "architecture": platform.machine(),
                "numpy_version": np.__version__,
                "pyspharm_conda_distribution": _package_version("pyspharm"),
                "pyspharm_cpp_distribution": _package_version("pyspharm-cpp"),
                "git_commit": subprocess.check_output(
                    ["git", "-C", str(root), "rev-parse", "HEAD"], text=True
                ).strip(),
                "git_worktree_clean": _git_worktree_clean(root),
            }
        )
        if sys.version_info[:2] != (3, 14):
            raise RuntimeError(f"validation requires CPython 3.14, got {sys.version}")
        if platform.python_implementation() != "CPython":
            raise RuntimeError(
                f"validation requires CPython, got {platform.python_implementation()}"
            )
        if sys.platform != "linux" or platform.machine().lower() not in {
            "x86_64",
            "amd64",
        }:
            raise RuntimeError(
                "reference validation requires Linux x86-64; got "
                f"{sys.platform}/{platform.machine()}"
            )

        import spharm

        import spharm_cpp
        from spharm_cpp import _spharm_cpp

        REFERENCE = spharm
        if getattr(spharm, "__version__", None) != "1.0.9":
            raise RuntimeError(
                f"unexpected spharm module version {getattr(spharm, '__version__', None)!r}"
            )
        reference_path = Path(spharm.__file__).resolve()
        native_cpp_path = Path(_spharm_cpp.__file__).resolve()
        site_packages = Path(sysconfig.get_path("purelib")).resolve()
        if not reference_path.is_file():
            raise RuntimeError(f"spharm has no source file: {spharm.__file__!r}")
        if site_packages not in reference_path.parents:
            raise RuntimeError(
                "reference import shadowing detected: spharm must come from "
                f"the active environment's site-packages; got {reference_path}"
            )
        if not native_cpp_path.name.endswith((".so", ".pyd")):
            raise RuntimeError(f"C++ backend is not native: {native_cpp_path}")
        if importlib.metadata.version("pyspharm") != "1.0.9":
            raise RuntimeError("validation requires conda-forge pyspharm==1.0.9")
        if importlib.metadata.version("pyspharm-cpp") != spharm_cpp.__version__:
            raise RuntimeError(
                "spharm_cpp.__version__ does not match distribution metadata"
            )

        source_revision = subprocess.check_output(
            ["git", "-C", str(root), "rev-parse", "HEAD"], text=True
        ).strip()
        ducc_revision = subprocess.check_output(
            ["git", "-C", str(root / "third_party" / "ducc"), "rev-parse", "HEAD"],
            text=True,
        ).strip()
        if len(ducc_revision) != 40:
            raise RuntimeError(
                f"could not determine pinned DUCC revision: {ducc_revision!r}"
            )

        ENVIRONMENT.update(
            {
                "reference_module_version": spharm.__version__,
                "reference_module_path": str(reference_path),
                "pyspharm_upstream_commit": "a464ebfed334188fc0ef275fad4aafd2d9fa716e",
                "cpp_module_version": spharm_cpp.__version__,
                "cpp_module_path": str(Path(spharm_cpp.__file__).resolve()),
                "cpp_native_path": str(native_cpp_path),
                "ducc_revision": ducc_revision,
                "cpu_info": spharm_cpp.cpu_info(),
                "module_imports_independent": True,
            }
        )

        print("PySpharm parity environment:")
        print(
            f"  Python: {sys.version.split()[0]} ({platform.python_implementation()})"
        )
        print(f"  NumPy: {np.__version__}")
        print(f"  pyspharm: {importlib.metadata.version('pyspharm')}")
        print(f"  spharm: {spharm.__version__} at {reference_path}")
        print(
            f"  pyspharm-cpp: {spharm_cpp.__version__} at {Path(spharm_cpp.__file__).resolve()}"
        )
        print(f"  native extension: {native_cpp_path}")
        print(f"  source commit: {source_revision}")
        print(f"  DUCC revision: {ducc_revision}")
        print(f"  CPU profile: {spharm_cpp.cpu_info()}", flush=True)
    except Exception as error:  # Fail the run visibly; never skip parity tests.
        REFERENCE_IMPORT_ERROR = f"{type(error).__name__}: {error}"
        pytest.exit(
            f"BLOCKED: parity environment check failed: {REFERENCE_IMPORT_ERROR}",
            returncode=pytest.ExitCode.TESTS_FAILED,
        )


def pytest_runtest_makereport(item: pytest.Item, call: pytest.CallInfo[Any]):
    outcome = yield
    report = outcome.get_result()
    if report.when != "call":
        return
    status = "PASS" if report.passed else ("BLOCKED" if report.skipped else "FAIL")
    entry = {"test": item.nodeid, "status": status}
    if report.failed:
        entry["detail"] = report.longreprtext[:4000]
    for marker in item.iter_markers("contract"):
        for name in marker.args:
            CONTRACT_STATUSES.setdefault(str(name), []).append(entry)
    for marker in item.iter_markers("approved_difference"):
        for name in marker.args:
            difference_status = "ACCEPTED_DIFFERENCE" if report.passed else status
            DIFFERENCE_STATUSES.setdefault(str(name), []).append(
                {**entry, "status": difference_status}
            )
    for marker in item.iter_markers("unapproved_difference"):
        for name in marker.args:
            UNAPPROVED_DIFFERENCE_STATUSES.setdefault(str(name), []).append(entry)


pytest_runtest_makereport = pytest.hookimpl(hookwrapper=True)(pytest_runtest_makereport)


def _rollup(entries: list[dict[str, Any]]) -> dict[str, Any]:
    order = {"PASS": 0, "ACCEPTED_DIFFERENCE": 1, "BLOCKED": 2, "FAIL": 3}
    status = max((entry["status"] for entry in entries), key=order.__getitem__)
    return {"status": status, "tests": entries}


def pytest_sessionfinish(session: pytest.Session, exitstatus: int) -> None:
    contracts = {
        name: _rollup(CONTRACT_STATUSES[name])
        if CONTRACT_STATUSES.get(name)
        else {"status": "BLOCKED", "tests": [], "detail": "no parity test collected"}
        for name in PUBLIC_CONTRACTS
    }
    differences = {
        name: _rollup(DIFFERENCE_STATUSES[name])
        if DIFFERENCE_STATUSES.get(name)
        else {
            "status": "BLOCKED",
            "tests": [],
            "detail": "approved difference has no behavior test",
        }
        for name in APPROVED_DIFFERENCES
    }
    statuses = [item["status"] for item in contracts.values()]
    statuses.extend(item["status"] for item in differences.values())
    unresolved = {}
    for name in sorted(
        set(UNAPPROVED_DIFFERENCES) | set(UNAPPROVED_DIFFERENCE_STATUSES)
    ):
        detail = UNAPPROVED_DIFFERENCES.get(
            name, "unapproved difference is absent from the validation manifest"
        )
        tests = UNAPPROVED_DIFFERENCE_STATUSES.get(name, [])
        unresolved[name] = {
            "status": "FAIL" if tests else "BLOCKED",
            "detail": detail,
            "tests": tests,
        }
    statuses.extend(item["status"] for item in unresolved.values())
    unknown_contracts = {
        name: _rollup(entries)
        for name, entries in CONTRACT_STATUSES.items()
        if name not in PUBLIC_CONTRACTS
    }
    unknown_approved_differences = {
        name: _rollup(entries)
        for name, entries in DIFFERENCE_STATUSES.items()
        if name not in APPROVED_DIFFERENCES
    }
    statuses.extend(["FAIL"] * len(unknown_contracts))
    statuses.extend(["FAIL"] * len(unknown_approved_differences))
    counts = {
        status: statuses.count(status)
        for status in ("PASS", "ACCEPTED_DIFFERENCE", "FAIL", "BLOCKED")
    }
    comparison_failures = [item for item in COMPARISONS if item.get("status") == "FAIL"]
    if counts["FAIL"] or counts["BLOCKED"] or comparison_failures:
        session.exitstatus = pytest.ExitCode.TESTS_FAILED
    comparison_counts = {
        status: sum(item.get("status") == status for item in COMPARISONS)
        for status in ("PASS", "ACCEPTED_DIFFERENCE", "FAIL", "BLOCKED")
    }
    payload = {
        "schema_version": 1,
        "git_commit": ENVIRONMENT.get("git_commit"),
        "environment": ENVIRONMENT,
        "reference_import_error": REFERENCE_IMPORT_ERROR,
        "contract_results": contracts,
        "approved_difference_results": differences,
        "approved_differences": APPROVED_DIFFERENCES,
        "unapproved_difference_results": unresolved,
        "unapproved_differences": UNAPPROVED_DIFFERENCES,
        "unknown_contract_results": unknown_contracts,
        "unknown_approved_difference_results": unknown_approved_differences,
        "comparisons": COMPARISONS,
        "investigations": INVESTIGATIONS,
        "counts": counts,
        "comparison_counts": comparison_counts,
        "pytest_exit_status": int(session.exitstatus),
    }
    try:
        (_report_directory() / "validation-results.json").write_text(
            json.dumps(payload, indent=2, sort_keys=True, allow_nan=False) + "\n",
            encoding="utf-8",
        )
    except (OSError, TypeError, ValueError) as error:
        session.exitstatus = pytest.ExitCode.TESTS_FAILED
        print(
            f"Could not write parity validation results: {type(error).__name__}: {error}"
        )
    print(
        "Parity contract status: "
        f"{counts['PASS']} PASS, "
        f"{counts['ACCEPTED_DIFFERENCE']} accepted differences, "
        f"{counts['FAIL']} FAIL, {counts['BLOCKED']} BLOCKED"
    )


def record_investigation(name: str, result: dict[str, Any]) -> None:
    INVESTIGATIONS[name] = result


def assert_close(
    operation: str,
    precision: str,
    actual: Any,
    expected: Any,
    *,
    rtol: float,
    atol: float,
) -> None:
    if not record_comparison(
        operation, precision, actual, expected, rtol=rtol, atol=atol
    ):
        comparison = COMPARISONS[-1]
        if comparison["maximum_absolute_error"] is None:
            raise AssertionError(
                f"{operation} shape mismatch: {comparison['actual_shape']} != "
                f"{comparison['expected_shape']}"
            )
        relative_error = comparison["maximum_relative_error"]
        relative_text = (
            "not meaningful below the absolute-error threshold"
            if relative_error is None
            else f"{relative_error:.9g}"
        )
        raise AssertionError(
            f"{operation} {precision} max abs="
            f"{comparison['maximum_absolute_error']:.9g}, max rel={relative_text}, "
            f"rtol={rtol:g}, atol={atol:g}"
        )


def record_comparison(
    operation: str,
    precision: str,
    actual: Any,
    expected: Any,
    *,
    rtol: float,
    atol: float,
    approved_difference: str | None = None,
) -> bool:
    actual_array = np.asarray(actual)
    expected_array = np.asarray(expected)
    if actual_array.shape != expected_array.shape:
        COMPARISONS.append(
            {
                "operation": operation,
                "precision": precision,
                "status": "FAIL",
                "actual_shape": list(actual_array.shape),
                "expected_shape": list(expected_array.shape),
                "maximum_absolute_error": None,
                "maximum_relative_error": None,
                "approved_difference": approved_difference,
                "rtol": rtol,
                "atol": atol,
            }
        )
        return False
    difference = np.abs(actual_array - expected_array)
    absolute_error = float(np.max(difference, initial=0.0))
    meaningful = np.abs(expected_array) > atol
    relative_error = None
    if np.any(meaningful):
        relative_error = float(
            np.max(difference[meaningful] / np.abs(expected_array[meaningful]))
        )
    passed = bool(np.allclose(actual_array, expected_array, rtol=rtol, atol=atol))
    COMPARISONS.append(
        {
            "operation": operation,
            "precision": precision,
            "status": (
                "PASS"
                if passed
                else "ACCEPTED_DIFFERENCE"
                if approved_difference in APPROVED_DIFFERENCES
                else "FAIL"
            ),
            "maximum_absolute_error": absolute_error,
            "maximum_relative_error": relative_error,
            "relative_error_threshold": atol,
            "approved_difference": approved_difference,
            "rtol": rtol,
            "atol": atol,
            "shape": list(actual_array.shape),
        }
    )
    return passed
