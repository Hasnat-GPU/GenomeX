"""The version is written in three files; they must agree.

`genomex.json` records `genomex_version` in every run, and `provenance.json` is
meant to make a result reproducible. A version string that disagrees with the
tag it was released under makes that record actively misleading — the run claims
a version whose code is something else. Three hand-edited copies drift; this
pins them together.
"""

import re
from pathlib import Path

import genomex

ROOT = Path(__file__).resolve().parents[1]
SEMVER = re.compile(r"^\d+\.\d+\.\d+$")


def _pyproject_version() -> str:
    text = (ROOT / "pyproject.toml").read_text(encoding="utf-8")
    m = re.search(r'^version\s*=\s*"([^"]+)"', text, re.MULTILINE)
    assert m, "pyproject.toml has no version"
    return m.group(1)


def _citation_version() -> str:
    text = (ROOT / "CITATION.cff").read_text(encoding="utf-8")
    m = re.search(r"^version:\s*(\S+)", text, re.MULTILINE)
    assert m, "CITATION.cff has no version"
    return m.group(1).strip('"')


def test_version_is_semver():
    assert SEMVER.match(genomex.__version__), genomex.__version__


def test_all_three_version_strings_agree():
    assert genomex.__version__ == _pyproject_version() == _citation_version(), (
        f"drift: package {genomex.__version__}, pyproject {_pyproject_version()}, "
        f"citation {_citation_version()}"
    )


def test_changelog_documents_the_current_version():
    """A release whose changelog stops at the previous version tells a reader
    the code has not changed since then."""
    text = (ROOT / "CHANGELOG.md").read_text(encoding="utf-8")
    assert f"## [{genomex.__version__}]" in text, (
        f"CHANGELOG.md has no section for {genomex.__version__}"
    )


def test_changelog_links_resolve_to_a_definition():
    """`[0.2.0]` in a heading with no link definition renders as literal brackets."""
    text = (ROOT / "CHANGELOG.md").read_text(encoding="utf-8")
    headings = set(re.findall(r"^## \[([^\]]+)\]", text, re.MULTILINE))
    defined = set(re.findall(r"^\[([^\]]+)\]:\s*http", text, re.MULTILINE))
    assert headings <= defined, f"undefined changelog links: {sorted(headings - defined)}"
