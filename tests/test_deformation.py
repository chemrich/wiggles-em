"""Tests for deformation_view — the well-supported tool, and its one boundary.

Motion survives blind testing; populations do not. So the assertions here are
about getting the motion right, plus one that matters more than it looks: the
report must never let an arrow be read as a population. An arrow is a statement
about where density moved, not about how many particles moved that way, and
nothing in the picture distinguishes those.
"""

from __future__ import annotations

import pytest
from conftest import draw, make_atoms

from wiggles_em.deformation import (
    deformation_view,
    per_atom_displacement,
    per_residue_displacement,
    read_uncertainty_table,
)
from wiggles_em.port import PortError
from wiggles_em.scene import Arrows, ColorByScalar

# chain, resi, resn, name, alt, q, b
ROWS = [
    ("A", "1", "ALA", "N", "", 1.0, 20.0),
    ("A", "1", "ALA", "CA", "", 1.0, 20.0),
    ("A", "2", "GLY", "CA", "", 1.0, 30.0),
    ("A", "3", "SER", "CA", "", 1.0, 40.0),
]


def _coords(shifts):
    """State-1 coords at the origin, state-2 shifted along x by ``shifts``."""
    start = [(0.0, float(i), 0.0) for i in range(len(ROWS))]
    end = [(dx, float(i), 0.0) for i, dx in enumerate(shifts)]
    return start, end


def _view(shifts, *, n_states=2, **kwargs):
    """Run deformation_view over a shift pattern and render what it returned."""
    start, end = _coords(shifts)
    kwargs.setdefault("start_state", 1)
    kwargs.setdefault("end_state", min(2, n_states))
    return draw(deformation_view, ROWS, start, end, "obj", n_states, **kwargs)


# ── the arithmetic ───────────────────────────────────────────────────────────


def test_displacement_is_the_distance_moved():
    got = per_atom_displacement([(0.0, 0.0, 0.0)], [(3.0, 4.0, 0.0)])
    assert got == pytest.approx([5.0])


def test_mismatched_atom_counts_are_refused_not_approximated():
    """The same refusal morph_states makes, for the same reason: unpaired
    atoms make the field meaningless rather than merely imprecise."""
    with pytest.raises(PortError, match="cannot be paired"):
        per_atom_displacement([(0.0, 0.0, 0.0)], [(1.0, 0.0, 0.0), (2.0, 0.0, 0.0)])


def test_the_refusal_points_at_the_tool_that_needs_no_pairing():
    with pytest.raises(PortError, match="ensemble_spread_view"):
        per_atom_displacement([(0.0, 0.0, 0.0)], [])


def test_per_residue_is_the_mean_over_its_atoms():
    atoms = list(ROWS)

    got = per_residue_displacement(make_atoms(ROWS), [2.0, 4.0, 1.0, 9.0])
    assert got[("A", "1")] == pytest.approx(3.0)
    assert got[("A", "2")] == pytest.approx(1.0)
    assert atoms  # the fixture really did hand back rows


def test_diverged_lengths_are_refused():

    atoms = make_atoms(ROWS)
    with pytest.raises(PortError, match="orders have diverged"):
        per_residue_displacement(atoms, [1.0])


# ── refusals ─────────────────────────────────────────────────────────────────


def test_a_single_state_object_is_refused():
    d = _view([0.0] * 4, n_states=1)
    out = d.report
    assert "REFUSED" in out
    assert "nothing here to difference" in out


def test_an_out_of_range_state_is_refused():
    with pytest.raises(PortError, match=r"outside the object's 1\.\.2 states"):
        _view([1.0] * 4, end_state=7)


def test_a_state_against_itself_is_refused():
    with pytest.raises(PortError, match="zero everywhere"):
        _view([1.0] * 4, start_state=1, end_state=1)


# ── the view ─────────────────────────────────────────────────────────────────


def test_displacement_reaches_the_bfactor_column():
    d = _view([0.0, 0.0, 3.0, 6.0], arrows=False)

    # On the scene: the values are the displacement, keyed by residue.
    (op,) = d.scene.of(ColorByScalar)
    by_residue = dict(zip(op.field.keys, op.field.values, strict=True))
    assert by_residue[("A", "2")] == pytest.approx(3.0)
    assert by_residue[("A", "3")] == pytest.approx(6.0)

    # And in the lowering: one dictionary push, not one alter per residue.
    pushed = next(c for c in d.port.commands if "stored.wiggles_scalar" in c)
    assert "3.0" in pushed and "6.0" in pushed


def test_the_report_names_the_biggest_movers():
    d = _view([0.0, 0.0, 3.0, 6.0], arrows=False)
    out = d.report
    assert "A/3 6.00Å" in out
    assert "Moved most" in out


def test_bfactors_are_stashed_before_being_overwritten():
    d = _view([1.0] * 4, arrows=False)
    # The caveat is the backend's, not the view's: protean re-sends a display
    # copy and has nothing to restore, so it must not inherit the advice.
    assert "restore_bfactors" in d.notes
    assert "B-factor" not in d.report


def test_arrows_are_drawn_as_cgo():
    d = _view([0.0, 0.0, 3.0, 6.0])

    # The scene carries geometry, not opcodes — Mol* has no CGO.
    (op,) = d.scene.of(Arrows)
    assert op.name == "obj_arrows"
    assert op.segments and all(a.start != a.end for a in op.segments)

    calls = d.port.calls("load_cgo")
    assert len(calls) == 1, d.port.call_log
    cgo, name = calls[0][0]
    assert name == "obj_arrows"
    assert cgo[0] == 9.0, "a cylinder opcode starts the shaft"
    assert 27.0 in cgo, "a cone opcode makes the head"


def test_only_ca_atoms_get_arrows():
    """One arrow per residue, at its CA — an arrow per atom is a thicket."""
    out = _view([5.0, 5.0, 5.0, 5.0]).report
    assert "3 drawn" in out, out


def test_motions_below_the_noise_floor_get_no_arrow():
    out = _view([0.0, 0.01, 0.01, 5.0]).report

    assert "1 drawn" in out
    assert "moved less than" in out
    assert "assert a motion nobody measured" in out


def test_the_arrow_ceiling_is_reported_never_silent():
    out = _view([0.0, 5.0, 6.0, 7.0], max_arrows=1).report

    assert "2 further residue(s)" in out
    assert "truncation of the largest movers" in out


def test_exaggeration_is_labelled():
    d = _view([0.0, 0.0, 3.0, 6.0], arrow_scale=3.0)
    out = d.report
    assert "EXAGGERATED 3x" in out
    assert "no longer" in out


def test_unexaggerated_output_says_nothing_about_exaggeration():
    d = _view([0.0, 0.0, 3.0, 6.0])
    out = d.report
    assert "EXAGGERATED" not in out


# ── half-set uncertainty ─────────────────────────────────────────────────────


def test_missing_uncertainty_is_stated_loudly():
    """A uniformly orange field is not evidence of uniform reliability."""
    d = _view([0.0, 0.0, 3.0, 6.0])
    out = d.report

    assert "NO UNCERTAINTY SUPPLIED" in out
    assert "do not have identical confidence" in out


def test_an_uncertainty_table_is_read_and_reported(tmp_path):
    path = tmp_path / "unc.txt"
    path.write_text("# chain resi value\nA 2 0.1\nA 3 0.9\n")

    out = _view([0.0, 0.0, 3.0, 6.0], uncertainty=str(path)).report
    assert "coloured by half-set uncertainty" in out
    assert "2 residues" in out
    assert "NO UNCERTAINTY SUPPLIED" not in out


def test_uncertainty_changes_the_arrow_colours(tmp_path):
    """Not just "a table was accepted" — the geometry must be identical and the
    colours must not be, or the table is being read and thrown away."""
    path = tmp_path / "unc.txt"
    path.write_text("A 2 0.0\nA 3 1.0\n")

    plain = _view([0.0, 0.0, 3.0, 6.0])
    coloured = _view([0.0, 0.0, 3.0, 6.0], uncertainty=str(path))

    # On the scene the claim is direct: same geometry, different colours.
    pa = plain.scene.of(Arrows)[0].segments
    ca = coloured.scene.of(Arrows)[0].segments
    assert [(a.start, a.end) for a in pa] == [(a.start, a.end) for a in ca]
    assert [a.colour for a in pa] != [a.colour for a in ca]

    a = plain.port.calls("load_cgo")[0][0][0]
    b = coloured.port.calls("load_cgo")[0][0][0]
    assert len(a) == len(b), "same arrows, so the same CGO length"
    assert a != b, "the uncertainty table should have changed the colours"
    # The shaft endpoints are geometry and must be untouched; CYLINDER lays out
    # opcode + 6 coordinates before the radius and colours begin.
    assert a[1:7] == b[1:7], "geometry must not depend on the uncertainty table"


def test_separators_are_flexible(tmp_path):
    path = tmp_path / "unc.csv"
    path.write_text("A,2,0.4\nA\t3\t0.6\n")
    assert read_uncertainty_table(path) == {("A", "2"): 0.4, ("A", "3"): 0.6}


def test_an_unparseable_table_says_what_it_wanted(tmp_path):
    path = tmp_path / "unc.txt"
    path.write_text("nonsense\n")
    with pytest.raises(ValueError, match="chain resi value"):
        read_uncertainty_table(path)


# ── the claim boundary ───────────────────────────────────────────────────────


def test_the_report_forbids_reading_a_population_off_an_arrow():
    d = _view([0.0, 0.0, 3.0, 6.0])
    out = d.report

    assert "MOTION IS THE WELL-SUPPORTED CLAIM" in out
    assert "Do NOT read them as: this fraction of particles moved" in out


# ── the states the coordinates actually came from ───────────────────────────


def test_the_reported_states_must_be_the_ones_the_coords_came_from():
    """Before the seam the view fetched the coordinates itself, so the state
    numbers it printed were the states it had read. Now the host supplies both
    and nothing links them: an `end_state` that defaults to the last state
    describes a 1->20 transition over a 1->2 displacement, and the user cites a
    magnitude for a motion that was never computed."""
    start, end = _coords([0.0, 0.0, 3.0, 6.0])
    with pytest.raises(TypeError):
        # The host read states 1 and 2 of a 20-state object; it must say so
        # rather than letting the view guess.
        deformation_view(make_atoms(ROWS), start, end, "obj", 20)


def test_states_are_reported_as_given():
    start, end = _coords([0.0, 0.0, 3.0, 6.0])
    report, _ = deformation_view(
        make_atoms(ROWS), start, end, "obj", 20, start_state=1, end_state=2
    )
    assert "state 1 -> 2" in report
    assert "1 -> 20" not in report
