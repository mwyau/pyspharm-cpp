<!-- SPDX-FileCopyrightText: 2026 Albert M. W. Yau -->
<!-- SPDX-License-Identifier: GPL-2.0-or-later -->

# Contributing

Contributions and numerical bug reports are welcome. For numerical changes, include a small reproducible case and, where possible, an analytic result or an independent reference comparison. Numerical background and software references are listed in [the documentation](docs/references.md). Describe the grid, truncation, precision, and expected behavior.

## Development

Initialize the ducc0 submodule, install the locked environment, and run the tests:

```sh
git submodule update --init --recursive
uv sync --locked --python 3.14
uv run --locked python -m pytest -q
```

Before submitting a pull request, run the configured checks:

```sh
uv run --locked prek run --all-files
```

CI additionally tests supported Python versions and Linux build profiles, checks documentation, and measures native C++ coverage. After a push to `main`, or on manual dispatch, CI also builds and tests Linux, macOS, and Windows wheels. The release workflow independently builds distribution wheels from the version tag.

## Code and numerical conventions

First-party C++ code is in `src/cpp/`. ducc0 is an unmodified Git submodule. The [numerical conventions](docs/conventions.md) describe coordinate order, normalization, and truncation. Preserve both precision modes and add regression tests for changed numerical behavior.

[AGENTS.md](AGENTS.md) records repository-specific implementation constraints. Keep documentation focused on behavior that users can observe.

## Releases

After the main-branch CI and PySpharm parity checks pass for a commit, push a version tag. The [publishing workflow](.github/workflows/python-publish.yml) builds the distributions automatically; `vX.Y.Z.devN` tags publish to TestPyPI and `vX.Y.Z` tags publish to PyPI. Publication waits for approval in the corresponding GitHub environment.

Configure required reviewers for both `pypi` and `testpypi` environments, enable release immutability in the repository settings, and protect `v*` tags against updates or deletion. Register `python-publish.yml` as a [Trusted Publisher](https://docs.pypi.org/trusted-publishers/adding-a-publisher/) on each package index with its matching environment name. After publishing, the workflow creates a GitHub Release containing the wheels and source distribution. Do not add publishing credentials to the repository.

Changes to the public API should update the generated typing stub, [API documentation](docs/api.md), and [PySpharm comparison](docs/pyspharm-compatibility.md) when appropriate.
