<!-- SPDX-License-Identifier: GPL-2.0-or-later -->

# Contributing

## Development environment

Initialize the pinned DUCC submodule and install the locked development
environment with `uv`:

```sh
git submodule update --init --recursive
uv sync --locked --python 3.14 --all-groups
```

Run the full test suite and checks before submitting changes:

```sh
uv run --locked python -m pytest -q
uv run --locked ruff check src tests typing docs .github/check_spdx.py
uv run --locked ruff format --check src tests typing docs .github/check_spdx.py
uv run --locked ty check typing/consumer.py
uv run --locked clang-format --dry-run --Werror src/cpp/*.cpp src/cpp/*.hpp
uv lock --check
```

The `quality` dependency group contains Ruff, Ty, clang-format, mdformat,
`prek` and zizmor. Sphinx and its theme are in `docs`; gcovr is in `coverage`.
The CI workflows run the SPDX header check, documentation build, C++ analysis,
coverage and packaging checks. SPDX checks are kept in CI rather than as a
pre-commit hook.

The Read the Docs configuration targets
<https://pyspharm-cpp.readthedocs.io/>. Register `mwyau/pyspharm-cpp` with
Read the Docs and verify a hosted build before linking the site as published
documentation. The `.readthedocs.yaml` file fetches the pinned
DUCC submodule and installs the locked build and documentation groups.

The native coverage workflow tests an instrumented wheel, generates LCOV
and XML reports, and uploads LCOV to Codecov. Codecov rejected the
October 24 upload with HTTP 400 because it received no valid token; the native
tests and report generation passed. Register the repository with Codecov and
configure a valid `CODECOV_TOKEN` Actions secret. The workflow keeps
`fail_ci_if_error: true` and uploads the coverage files as Actions artifacts.
Do not commit the token.

## Source and numerical changes

Keep first-party C++ under `src/cpp/`. Do not edit or reformat `third_party/ducc`;
it is an unmodified submodule at the commit recorded by Git. Follow the
constraints in [`AGENTS.md`](AGENTS.md), especially the regular-grid analysis
limit, longitude sampling limit, coefficient layout and use of DUCC's
transforms and SIMD abstractions.

Numerical changes should add an independent reference or analytic test where
possible. Preserve both float32/complex64 and float64/complex128 behavior,
strided and read-only inputs, and the free-threaded concurrency tests. Keep
benchmarks descriptive; measured speedups are not correctness criteria.

When changing a nanobind signature or docstring, regenerate the stub and check
the result:

```sh
uv run --locked python -m nanobind.stubgen \
  -m spharm_cpp._spharm_cpp \
  -O /tmp/pyspharm-stub \
  -M /tmp/pyspharm-stub/py.typed \
  -p .github/stubgen-patterns.txt
diff -u src/spharm_cpp/_spharm_cpp.pyi /tmp/pyspharm-stub/_spharm_cpp.pyi
```

## Pull requests and releases

Keep workflow permissions minimal, pin each GitHub Action to an exact
published patch-version tag (for example, `actions/checkout@v7.0.1`), and use
maintained runner labels in the release matrix. Do not use moving major-version
tags such as `@v7` or branch references such as `@main` or `@master`. Release workflows must preserve recursive DUCC submodule checkout,
disabled checkout credentials, Python dependency locking and
architecture-specific wheel checks.

Do not publish packages, create release tags, or enable unattended publishing
as part of ordinary development. See [`CHANGELOG.md`](CHANGELOG.md) for release
notes and [`CITATION.cff`](CITATION.cff) for citation metadata.

### Trusted Publishing setup

The release workflow is manual and accepts only version tags. It routes
`vX.Y.Z.devN` to TestPyPI and `vX.Y.Z` to PyPI, then checks the selected index
against the tag and project version before publication.

Before the first package upload, configure both package indexes:

1. On GitHub, open **Settings → Environments** and create `testpypi` and
   `pypi`. Add trusted maintainers or a release team as **Required reviewers**
   and enable **Prevent self-review**. Restrict deployment refs to version
   tags (`v*`). Keep these approvals enabled so a tag cannot publish without
   a human review.
1. On TestPyPI, open the `pyspharm-cpp` project’s **Manage → Publishing** page
   and add a GitHub Actions trusted publisher with owner `mwyau`, repository
   `pyspharm-cpp`, workflow `python-publish.yml`, and environment `testpypi`. If the
   project does not exist yet, create its pending publisher from account
   settings.
1. On PyPI, configure the same publisher values with environment `pypi`.
1. In **Settings → Rules → Rulesets**, protect version tags (`v*`) so only
   maintainers can create or update them.

PyPI and TestPyPI use separate publisher registrations and environments. The
workflow uses GitHub OIDC and stores no package-index tokens. For details, see
[PyPI’s Trusted Publisher setup](https://docs.pypi.org/trusted-publishers/adding-a-publisher/)
and [GitHub’s environment protection settings](https://docs.github.com/en/actions/how-tos/deploy/configure-and-manage-deployments/manage-environments).
