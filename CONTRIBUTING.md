<!-- SPDX-License-Identifier: GPL-2.0-or-later -->

# Contributing

Contributions and numerical bug reports are welcome. For numerical changes, include a small reproducible case and, where possible, an analytic result or an independent reference comparison. Numerical background and software references are listed in [the documentation](docs/references.md). Describe the grid, truncation, precision, and expected behavior.

## Development

Initialize DUCC, install the locked environment, and run the tests:

```sh
git submodule update --init --recursive
uv sync --locked --python 3.14
uv run --locked python -m pytest -q
```

Before submitting a pull request, run the configured checks:

```sh
uv run --locked prek run --all-files
```

CI additionally builds the extension, tests supported Python versions and build profiles, checks documentation, and measures native C++ coverage.

## Code and numerical conventions

First-party C++ code is in `src/cpp/`. DUCC is an unmodified Git submodule. The [numerical conventions](docs/conventions.md) describe coordinate order, normalization, and truncation. Preserve both precision modes and add regression tests for changed numerical behavior.

[AGENTS.md](AGENTS.md) records repository-specific implementation constraints. Keep documentation focused on behavior that users can observe.

## Releases

Package releases are manually triggered through GitHub Actions and require approval. See [the publishing workflow](.github/workflows/python-publish.yml) for checks and [PyPI's trusted-publisher instructions](https://docs.pypi.org/trusted-publishers/adding-a-publisher/) for initial setup. Do not add credentials to the repository.

Changes to the public API should update the generated typing stub, [API documentation](docs/api.md), and [PySpharm comparison](docs/pyspharm-compatibility.md) when appropriate.
