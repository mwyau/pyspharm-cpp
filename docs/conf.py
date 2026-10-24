# SPDX-License-Identifier: GPL-2.0-or-later
import sys
from datetime import date
from pathlib import Path
import tomllib

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

project = "pyspharm-cpp"
author = "Albert Yau"
copyright = f"{date.today().year}, {author}"
release = tomllib.loads((Path(__file__).resolve().parents[1] / "pyproject.toml").read_text())["project"]["version"]

extensions = [
    "myst_parser",
    "sphinx.ext.autodoc",
    "sphinx.ext.intersphinx",
    "sphinx.ext.mathjax",
    "sphinx.ext.napoleon",
]

source_suffix = {".md": "markdown"}
root_doc = "index"
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]
nitpicky = True
nitpick_ignore_regex = [(r"py:class", r"numpy(\..*)?")]
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
