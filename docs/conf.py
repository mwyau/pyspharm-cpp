# SPDX-FileCopyrightText: 2026 Albert M. W. Yau
# SPDX-License-Identifier: GPL-2.0-or-later
import sys
import tomllib
from datetime import date
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

project = "pyspharm-cpp"
author = "Albert M. W. Yau"
copyright = f"{date.today().year}, {author}"
pyproject = Path(__file__).resolve().parents[1] / "pyproject.toml"
release = tomllib.loads(pyproject.read_text())["project"]["version"]

extensions = [
    "myst_parser",
    "sphinx.ext.autodoc",
    "sphinx.ext.intersphinx",
    "sphinx.ext.mathjax",
    "sphinx.ext.napoleon",
]

templates_path = ["_templates"]
html_sidebars = {"**": ["site-collapse.html", "site-navigation.html"]}
source_suffix = {".md": "markdown"}
root_doc = "index"
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]
nitpicky = True
nitpick_ignore_regex = [(r"py:class", r"numpy(\..*)?")]
nitpick_ignore = [
    ("py:class", "ndarray"),
    ("py:class", "optional"),
    ("py:class", "numpy dtype"),
]
todo_include_todos = False

html_theme = "pydata_sphinx_theme"
html_title = "pyspharm-cpp"
html_static_path = []
html_theme_options = {
    "show_toc_level": 2,
    "navigation_depth": 3,
}

intersphinx_mapping = {
    "numpy": ("https://numpy.org/doc/stable/", None),
    "python": ("https://docs.python.org/3/", None),
}
