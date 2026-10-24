# SPDX-License-Identifier: GPL-2.0-or-later
"""Check that the generated docs include navigation and inline mathematics."""

from pathlib import Path

from bs4 import BeautifulSoup

pages = (
    "index",
    "quickstart",
    "api",
    "conventions",
    "pyspharm-compatibility",
    "benchmarks",
    "references",
)
all_links = {f"{page}.html" for page in pages[1:]}
root = Path("docs/_build/html")
for page in pages:
    soup = BeautifulSoup((root / f"{page}.html").read_text(), "html.parser")
    nav = soup.select_one("nav.bd-docs-nav")
    assert nav is not None, f"{page}: missing section navigation"
    links = {
        item.get("href", "").split("#")[0].rsplit("/", 1)[-1]
        for item in nav.select("a[href]")
    }
    required = all_links - {f"{page}.html"}
    assert required <= links, f"{page}: missing links {required - links}"

for page, formula in (
    ("conventions", r"\bar P_{00}"),
    ("conventions", r"\theta=\pi/2-\phi"),
    ("pyspharm-compatibility", r"e^{im\lambda_j}"),
):
    soup = BeautifulSoup((root / f"{page}.html").read_text(), "html.parser")
    content = " ".join(node.get_text() for node in soup.select("span.math"))
    assert formula in content, f"{page}: missing rendered math {formula}"
