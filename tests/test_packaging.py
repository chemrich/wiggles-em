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


# ── the declared tool surface ───────────────────────────────────────────────


def test_every_declared_tool_is_exported_and_callable():
    """`TOOLS` is what `SPEC.md` reconciles against, so a name in it that does
    not resolve would make the spec agree with nothing."""
    for name in wiggles_em.TOOLS:
        assert name in wiggles_em.__all__, f"{name} is in TOOLS but not exported"
        assert callable(getattr(wiggles_em, name, None)), f"{name} is not callable"


#: Exports that are deliberately not tools. Value types, conversions, the port
#: protocol, and plumbing that undoes another tool's side effect.
#:
#: The **exclusion** list is the maintained one, not the inclusion list. That
#: direction matters: forgetting to classify a new export fails the test below,
#: with a message naming it. Guarding only `_view` names instead let a
#: non-view tool — a second loader, an exporter — reach `__all__` unclassified,
#: and the first symptom was `check_spec.py` reporting "SPEC.md marks `X` as
#: BUILT, but wiggles-em does not offer it… the tick is wrong", sending the
#: maintainer to edit the spec when the defect was a missing TOOLS entry.
NOT_TOOLS = {
    "TOOLS",
    # value types
    "Atom",
    "Ensemble",
    "MapHeader",
    "MapStats",
    "DensityJob",
    "Method",
    "Populations",
    "Provenance",
    "StatsSource",
    "WeightSource",
    # the port protocol and its implementations
    "BridgePort",
    "FakePort",
    "PortError",
    "PymolPort",
    "SendRequestPort",
    # conversions and readers, below the level of a tool
    "contains_absence_claim",
    "read_density_job",
    "read_deconvolved_weights",
    "fetch_atoms",
    "grid_differences",
    "loaded_ensemble",
    "provenance_of",
    "read_map_header",
    "to_absolute",
    "to_sigma",
    "declare",
    # plumbing: undoes a side effect of three tools, argues for nothing itself
    "clear_stash",
    "restore_bfactors",
}


def test_every_export_is_classified_as_tool_or_not_tool():
    """The load-bearing half. A declared list rots the moment someone adds an
    entry point and forgets it, and the symptom would be `check_spec.py`
    reconciling a smaller surface than the package offers — the exact failure
    re-pointing it away from MCPymol was meant to end.

    Exact equality, so a new export must be put in one bucket or the other.
    """
    exported = set(wiggles_em.__all__)
    unclassified = exported - NOT_TOOLS - set(wiggles_em.TOOLS)
    assert not unclassified, (
        f"{sorted(unclassified)} exported but classified neither way. Add each "
        f"to wiggles_em.TOOLS if the research argued for it as a tool, or to "
        f"NOT_TOOLS here if it is a value type, a conversion or plumbing."
    )

    stale = set(wiggles_em.TOOLS) - exported
    assert not stale, f"{sorted(stale)} in TOOLS but no longer exported"


def test_the_declared_surface_holds_no_duplicates():
    assert len(wiggles_em.TOOLS) == len(set(wiggles_em.TOOLS))
