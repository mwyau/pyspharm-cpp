# SPDX-FileCopyrightText: 2026 Albert M. W. Yau
# SPDX-License-Identifier: GPL-2.0-or-later
"""Check SPDX copyright and license headers in first-party text files."""

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
HEADERS = ("SPDX-FileCopyrightText:", "SPDX-License-Identifier: GPL-2.0-or-later")
SUFFIXES = {
    ".c",
    ".cff",
    ".cpp",
    ".h",
    ".html",
    ".hpp",
    ".md",
    ".py",
    ".pyi",
    ".toml",
    ".txt",
    ".yaml",
    ".yml",
}
SPECIAL_FILES = {
    ".clang-format",
    ".clang-tidy",
    ".gitignore",
    ".gitmodules",
    "CMakeLists.txt",
}
SOURCE_DIRECTORIES = (
    ".github",
    "benchmarks",
    "docs",
    "src",
    "tests",
    "typing",
)
EXCLUDED_PARTS = {
    ".git",
    ".venv",
    "__pycache__",
    "build",
    "dist",
    "third_party",
}


def candidate_files() -> set[Path]:
    files = {
        path
        for path in ROOT.iterdir()
        if path.is_file() and (path.suffix in SUFFIXES or path.name in SPECIAL_FILES)
    }
    for directory in SOURCE_DIRECTORIES:
        base = ROOT / directory
        if base.is_dir():
            files.update(
                path
                for path in base.rglob("*")
                if path.is_file()
                and path.name != "py.typed"
                and path.suffix in SUFFIXES
                and not EXCLUDED_PARTS.intersection(path.relative_to(ROOT).parts)
            )
    return files


def missing_headers(path: Path) -> list[str]:
    try:
        lines = path.read_text(encoding="utf-8").splitlines()[:8]
    except (OSError, UnicodeDecodeError):
        return list(HEADERS)
    return [header for header in HEADERS if not any(header in line for line in lines)]


def main() -> int:
    missing = {
        path.relative_to(ROOT): headers
        for path in candidate_files()
        if (headers := missing_headers(path))
    }
    if missing:
        print("Missing SPDX headers in the first eight lines:", file=sys.stderr)
        for path, headers in sorted(missing.items()):
            print(f"  {path}: {', '.join(headers)}", file=sys.stderr)
        return 1
    print(f"SPDX headers present in {len(candidate_files())} first-party files.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
