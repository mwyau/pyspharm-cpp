# Contributing

For bug reports, include a reproducible example, grid type and dimensions,
spectral truncation, dtype, and expected and observed results. For numerical
differences, identify the analytic result or reference implementation used for
comparison.

## Development

Initialize the ducc0 submodule and install the locked environment:

```sh
git submodule update --init --recursive
uv sync --locked --python 3.14
```

Run the tests and configured checks before opening a pull request:

```sh
uv run --locked python -m pytest -q
uv run --locked prek run --all-files
```

## Changes

First-party C++ code is in `src/cpp/`; the ducc0 submodule is unmodified.
Follow the [numerical conventions](docs/conventions.md) for coordinates,
normalization, and truncation. Add regression tests for changed behavior,
including relevant grids and precisions. Use analytic results or independent
reference comparisons where possible.

For public API changes, update the generated typing stub and the relevant
[API documentation](docs/api.md) and [PySpharm comparison](docs/pyspharm-compatibility.md).
