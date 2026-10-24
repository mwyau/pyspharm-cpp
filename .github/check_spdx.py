# SPDX-License-Identifier: GPL-2.0-or-later
"""Check SPDX identifiers in first-party source and documentation files."""

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
IDENTIFIER = "SPDX-License-Identifier: GPL-2.0-or-later"
SUFFIXES = {
    ".c",
    ".cff",
    ".cpp",
    ".h",
    ".hpp",
    ".md",
    ".py",
    ".pyi",
    ".toml",
    ".yaml",
    ".yml",
}
SPECIAL_FILES = {
    ".clang-format",
    ".clang-tidy",
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


def has_identifier(path: Path) -> bool:
    try:
        lines = path.read_text(encoding="utf-8").splitlines()[:8]
    except (OSError, UnicodeDecodeError):
        return False
    return any(IDENTIFIER in line for line in lines)


def main() -> int:
    missing = sorted(
        path.relative_to(ROOT) for path in candidate_files() if not has_identifier(path)
    )
    if missing:
        print("Missing SPDX identifier in the first eight lines:", file=sys.stderr)
        for path in missing:
            print(f"  {path}", file=sys.stderr)
        return 1
    print(f"SPDX identifiers present in {len(candidate_files())} first-party files.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
