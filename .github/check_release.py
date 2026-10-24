# SPDX-FileCopyrightText: 2026 Albert M. W. Yau
# SPDX-License-Identifier: GPL-2.0-or-later
"""Verify a version-tag release and its prerequisite CI runs."""

import argparse
import json
import os
import re
import subprocess
import tomllib
from pathlib import Path


def gh_api(path: str) -> dict:
    url = f"repos/{os.environ['GITHUB_REPOSITORY']}/{path}"
    result = subprocess.check_output(["gh", "api", url], text=True)
    return json.loads(result)


def release_commit(tag: str) -> str:
    obj = gh_api(f"git/ref/tags/{tag}")["object"]
    seen = set()
    while obj["type"] == "tag":
        sha = obj["sha"]
        if sha in seen:
            raise SystemExit("Release tag contains a reference cycle")
        seen.add(sha)
        obj = gh_api(f"git/tags/{sha}")["object"]
    if obj["type"] != "commit":
        raise SystemExit("Release tag must resolve to a commit")
    return obj["sha"]


def check_latest_run(workflow: str, event: str, sha: str) -> None:
    data = gh_api(
        f"actions/workflows/{workflow}/runs?head_sha={sha}&event={event}&per_page=100"
    )
    runs = [
        run
        for run in data["workflow_runs"]
        if run["head_sha"] == sha
        and run["event"] == event
        and run["head_branch"] == "main"
    ]
    if not runs:
        raise SystemExit(f"No {workflow} {event} run for {sha}")
    latest = max(runs, key=lambda run: run["run_started_at"] or run["created_at"])
    if latest["status"] != "completed" or latest["conclusion"] != "success":
        raise SystemExit(
            f"{workflow} has no successful latest {event} run: {latest['html_url']}"
        )
    print(f"{workflow}: {latest['html_url']}")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--evidence", action="store_true")
    args = parser.parse_args()

    tag = os.environ["GITHUB_REF_NAME"]
    if os.environ["GITHUB_REF_TYPE"] != "tag":
        raise SystemExit("Release workflow must run on a version tag")
    if not re.fullmatch(
        r"v(?:0|[1-9]\d*)\.(?:0|[1-9]\d*)\.(?:0|[1-9]\d*)"
        r"(?:\.dev(?:0|[1-9]\d*))?",
        tag,
    ):
        raise SystemExit("Expected a tag of the form vX.Y.Z or vX.Y.Z.devN")
    version = tomllib.loads(Path("pyproject.toml").read_text())["project"]["version"]
    if tag != f"v{version}":
        raise SystemExit(f"Tag {tag} does not match package version {version}")
    sha = os.environ["GITHUB_SHA"]
    if release_commit(tag) != sha:
        raise SystemExit(f"Tag {tag} no longer resolves to {sha}")

    if args.evidence:
        index = "testpypi" if ".dev" in version else "pypi"
        print(f"Validated {tag} for {index}")
        check_latest_run("ci.yml", "push", sha)
        check_latest_run("validation.yml", "push", sha)


if __name__ == "__main__":
    main()
