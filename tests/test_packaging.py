"""The package tells type checkers that it is typed.

Without a `py.typed` marker a package is **invisible to mypy** however
thoroughly it is annotated: PEP 561 says a checker must ignore inline types
from an installed package that does not ship one. This package is fully
annotated and checked under `mypy` in its own CI, and consumers still saw
`Module "wiggles_em" has no attribute ...` or nothing at all — protean carries
an `ignore_missing_imports` override for exactly this reason.

**What this test proves and does not prove.** It proves the marker is present in
the source tree, so deleting it fails a test. It does *not* prove the marker
reaches the built wheel — that depends on the build backend, and asserting it
means building, which costs seconds against a suite that runs in one. Hatchling
includes every file under the configured package directory, and the wheel was
built and inspected once by hand when this landed:

    $ uv build --wheel && python -c "import zipfile,glob; \\
        print([n for n in zipfile.ZipFile(sorted(glob.glob('dist/*.whl'))[-1]) \\
               .namelist() if 'py.typed' in n])"
    ['wiggles_em/py.typed']

If the build configuration in `pyproject.toml` ever stops including package
data, this test will still pass and the marker will still be missing from the
wheel. That is the gap; it is recorded rather than papered over.
"""

from __future__ import annotations

from pathlib import Path

import wiggles_em


def test_the_py_typed_marker_is_present():
    package = Path(wiggles_em.__file__).parent
    marker = package / "py.typed"
    assert marker.exists(), (
        f"{marker} is missing. Without it every consumer's type checker ignores "
        f"this package's annotations entirely, however complete they are."
    )


def test_the_marker_is_empty():
    """PEP 561's marker is a flag, not a manifest. A `partial\\n` in it means
    something quite different — that the package's types are incomplete and a
    checker should keep looking in typeshed — and is easy to write by accident
    when a file is created with an editor that adds a newline."""
    package = Path(wiggles_em.__file__).parent
    assert (package / "py.typed").read_text() == ""
