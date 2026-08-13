"""What a lowered selection actually matches, against a running PyMOL.

    pytest -m live tests/test_selection_live.py

WARNING: this drives the PyMOL session you have open. It only adds and deletes
objects under the ``_wsel_`` prefix and touches nothing already loaded, but
``pytest -m live`` as a whole clears the session.

Why this file exists
--------------------
``test_backends.py`` asserted that ``render_selection(Sel.residues([("A",
"-3")]))`` contained ``resi "-3"`` — that the string had been quoted. It said
nothing about what PyMOL does with that string, and the answer turned out to be
"reads it as the range 1-3 anyway". The bug the quoting was added to fix was
still there, under a passing test, in this package and in MCPymol both, because
the test and the bug shared a mental model.

Everything above ``backends/`` is pure and testable with ``FakePort``, which is
the design working as intended — but it means nothing in the suite can see
PyMOL's *grammar*. This file is the seam where that has to be checked against
the real parser rather than reasoned about.

So it asserts **atom counts from a live session**: an object whose residue
numbering makes the range reading and the literal reading return different
answers, selected through the real lowering path, then counted.

The structure has chain A residues -3, 1 and 2 (two atoms each) and a
blank-chain residue 5. Selecting residue -3 must match 2. Under the range
reading it matches all 6.
"""

from __future__ import annotations

import tempfile
from pathlib import Path

import pytest

from wiggles_em.backends.pymol import render_selection
from wiggles_em.port import BridgePort, call
from wiggles_em.scene import Sel

pytestmark = pytest.mark.live

OBJ = "_wsel_probe"

# (chain, resseq, atom name) — two atoms per residue.
ROWS = [
    ("A", -3, "N"),
    ("A", -3, "CA"),
    ("A", 1, "N"),
    ("A", 1, "CA"),
    ("A", 2, "N"),
    ("A", 2, "CA"),
    ("", 5, "N"),
    ("", 5, "CA"),
]


def _pdb() -> str:
    lines = []
    for serial, (chain, resseq, name) in enumerate(ROWS, start=1):
        name_field = f" {name:<3}"
        lines.append(
            f"{'ATOM':<6}{serial:>5} {name_field:<4}{' ':1}{'ALA':>3} "
            f"{chain:1}{resseq:>4}{' ':1}   "
            f"{float(serial):>8.3f}{0.0:>8.3f}{0.0:>8.3f}"
            f"{1.0:>6.2f}{0.0:>6.2f}          {'C':>2}"
        )
    return "\n".join(lines) + "\nEND\n"


def _scoped(sel: Sel) -> str:
    """Lower ``sel`` and scope it to the probe object."""
    return f"({OBJ}) and {render_selection(sel)}"


@pytest.fixture
def probe():
    """The probe object, loaded into the live session and removed afterwards."""
    port = BridgePort()
    with tempfile.TemporaryDirectory() as tmp:
        path = Path(tmp) / "wsel.pdb"
        path.write_text(_pdb())
        call(port, "delete", OBJ)
        call(port, "load", str(path), OBJ)
        try:
            yield port
        finally:
            call(port, "delete", OBJ)


def test_the_probe_is_numbered_the_way_the_assertions_assume(probe):
    """If PyMOL ever stops reading these residue numbers as written, every
    other assertion in this file would pass or fail for the wrong reason."""
    assert call(probe, "count_atoms", f"({OBJ})") == 8
    assert call(probe, "count_atoms", f'({OBJ}) and chain "A"') == 6


def test_a_negative_residue_matches_only_itself(probe):
    """The whole point: two atoms, not the six that residues -3, 1 and 2 hold
    between them."""
    assert call(probe, "count_atoms", _scoped(Sel.residues([("A", "-3")]))) == 2


def test_the_prop_form_is_escaped_too(probe):
    """`Sel.prop` lowers through the same `quote`, and a fix that reached
    `Sel.residues` only is exactly how the last two review rounds went."""
    scoped = f'({OBJ}) and chain "A" and {render_selection(Sel.prop("resi", "-3"))}'

    assert call(probe, "count_atoms", scoped) == 2


def test_an_unescaped_negative_residue_really_does_over_match(probe):
    """The bug this guards against, asserted directly, so the tests above
    cannot quietly start passing for some unrelated reason.

    The first string is what this package emitted before the fix."""
    assert call(probe, "count_atoms", f'({OBJ}) and chain "A" and resi "-3"') == 6
    assert call(probe, "count_atoms", f'({OBJ}) and chain "A" and resi -3') == 6


def test_a_negative_residue_survives_the_plus_list(probe):
    """The compact `+` list is what keeps 1123 scored residues out of a
    1123-term disjunction. A quoted, escaped value has to compose into it, or
    the negative case needs its own path again."""
    sel = Sel.residues([("A", "-3"), ("A", "1")])

    assert "+" in render_selection(sel), render_selection(sel)
    assert call(probe, "count_atoms", _scoped(sel)) == 4


def test_ordinary_residues_are_unaffected(probe):
    for resi in ("1", "2"):
        assert call(probe, "count_atoms", _scoped(Sel.residues([("A", resi)]))) == 2


def test_a_blank_chain_stays_scoped_to_the_object(probe):
    """The other half of `quote`: `chain  and resi 5` would take `and` as the
    chain name and stop being scoped to the object."""
    assert call(probe, "count_atoms", _scoped(Sel.residues([("", "5")]))) == 2


def test_a_multi_chain_selection_matches_every_chain_it_names(probe):
    """The per-chain grouping builds one term per chain; a chain dropped or
    over-matched there is invisible to a single-chain test."""
    sel = Sel.residues([("A", "-3"), ("", "5")])

    assert call(probe, "count_atoms", _scoped(sel)) == 4


# ── atom identity under removal (REVIEW #2) ─────────────────────────────────


ROWS_WITH_HYDROGENS = [
    # (chain, resseq, name, element) -- hydrogens sort *between* the heavy
    # atoms in index order, which is what makes their removal renumber.
    ("A", 1, "N", "N"),
    ("A", 1, "H", "H"),
    ("A", 1, "CA", "C"),
    ("A", 1, "HA", "H"),
    ("A", 1, "C", "C"),
    ("A", 2, "N", "N"),
    ("A", 2, "CA", "C"),
    ("A", 2, "C", "C"),
]

HYDRO_OBJ = "_wsel_hydro"


def _hydro_pdb() -> str:
    lines = []
    for serial, (chain, resseq, name, element) in enumerate(ROWS_WITH_HYDROGENS, start=1):
        name_field = f" {name:<3}"
        lines.append(
            f"{'ATOM':<6}{serial:>5} {name_field:<4}{' ':1}{'ALA':>3} "
            f"{chain:1}{resseq:>4}{' ':1}   "
            f"{float(serial):>8.3f}{0.0:>8.3f}{0.0:>8.3f}"
            f"{1.0:>6.2f}{float(serial):>6.2f}          {element:>2}"
        )
    return "\n".join(lines) + "\nEND\n"


@pytest.fixture
def hydro():
    """A structure with hydrogens, loaded live and removed afterwards."""
    port = BridgePort()
    with tempfile.TemporaryDirectory() as tmp:
        path = Path(tmp) / "hydro.pdb"
        path.write_text(_hydro_pdb())
        call(port, "delete", HYDRO_OBJ)
        call(port, "load", str(path), HYDRO_OBJ)
        try:
            yield port
        finally:
            call(port, "delete", HYDRO_OBJ)


def _read(port, expr="(rank, index, b)"):
    """Rows keyed by b-factor, which is unique per atom here and never moves."""
    rows = call(port, "iterate_to_list", HYDRO_OBJ, expr)
    flat = [r[0] if isinstance(r, list) and len(r) == 1 and isinstance(r[0], list) else r for r in rows]
    return {row[2]: (row[0], row[1]) for row in flat}


def test_removal_renumbers_index_but_not_rank(hydro):
    """REVIEW #2, established by observation rather than from documentation.

    `Atom.key` keys on (model, rank). The reason has to be checked against a
    real session, because both fields are unique within an object and the
    uniqueness alone is what made `index` look adequate.
    """
    before = _read(hydro)
    call(hydro, "remove", f"{HYDRO_OBJ} and hydro")
    after = _read(hydro)

    heavy = sum(1 for *_rest, element in ROWS_WITH_HYDROGENS if element != "H")
    survivors = set(after) & set(before)
    assert len(survivors) == heavy, (
        f"expected {heavy} heavy atoms to survive, got {sorted(survivors)}"
    )

    moved_index = {b for b in survivors if before[b][1] != after[b][1]}
    moved_rank = {b for b in survivors if before[b][0] != after[b][0]}

    assert moved_index, (
        "no atom's index moved, so this structure does not exercise the "
        "finding at all — the test would pass for the wrong reason"
    )
    assert not moved_rank, (
        f"rank moved for atoms {sorted(moved_rank)}, so it is not the stable "
        f"field this package keys identity on. before={before} after={after}"
    )


def test_rank_is_unique_across_the_object(hydro):
    """Stability is worth nothing without uniqueness: a repeated key hands one
    atom another's value, which is the collision the key exists to remove."""
    ranks = [rank for rank, _index in _read(hydro).values()]

    assert len(set(ranks)) == len(ranks), ranks


def test_the_package_asks_for_the_field_it_keys_on(hydro):
    """ATOM_EXPR is what the session is actually queried with, so a mismatch
    between it and `Atom.key` would surface here rather than in a mock."""
    from wiggles_em.atoms import ATOM_EXPR, fetch_atoms

    assert "rank" in ATOM_EXPR
    atoms = fetch_atoms(hydro, HYDRO_OBJ)

    assert len(atoms) == len(ROWS_WITH_HYDROGENS)
    assert len({a.key for a in atoms}) == len(atoms), "keys collided in a live read"
