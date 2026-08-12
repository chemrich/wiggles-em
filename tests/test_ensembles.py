"""Tests for multi-state ensemble views."""

from __future__ import annotations

import pytest
from conftest import draw, render

from wiggles_em.atoms import Atom
from wiggles_em.bfactors import has_stash
from wiggles_em.ensembles import (
    SPREAD_LEGEND,
    ensemble_spread_view,
    morph_states,
    per_atom_spread,
    per_residue_spread,
)
from wiggles_em.port import FakePort, PortError
from wiggles_em.scene import ColorByScalar, Morph, SizeByScalar


def spread(coords, rows, **kw):
    """Run ensemble_spread_view over coordinates and rows, and render it."""
    return draw(ensemble_spread_view, rows, coords, "obj", **kw)


ROWS_2 = [
    ("A", "1", "MET", "CA", "", 1.0, 20.0),
    ("A", "2", "ALA", "CA", "", 1.0, 20.0),
]


# -- the maths -------------------------------------------------------------


def test_no_motion_gives_zero_spread():
    spread = per_atom_spread([[(0.0, 0.0, 0.0)], [(0.0, 0.0, 0.0)]])
    assert spread == [0.0]


def test_spread_is_rms_about_the_centroid():
    """Two states at ±1 along x: centroid 0, RMS deviation 1."""
    spread = per_atom_spread([[(-1.0, 0.0, 0.0)], [(1.0, 0.0, 0.0)]])
    assert spread[0] == pytest.approx(1.0)


def test_spread_needs_two_states():
    with pytest.raises(PortError, match="at least two states"):
        per_atom_spread([[(0.0, 0.0, 0.0)]])


def test_differing_atom_counts_make_spread_undefined():
    with pytest.raises(PortError, match="differing atom counts"):
        per_atom_spread([[(0.0, 0.0, 0.0)], [(0.0, 0.0, 0.0), (1.0, 0.0, 0.0)]])


def test_residue_spread_averages_its_atoms():
    atoms = [
        Atom("A", "1", "MET", "N", "", 1.0, 0.0),
        Atom("A", "1", "MET", "CA", "", 1.0, 0.0),
        Atom("A", "2", "ALA", "CA", "", 1.0, 0.0),
    ]
    out = per_residue_spread(atoms, [1.0, 3.0, 5.0])
    assert out[("A", "1")] == pytest.approx(2.0)
    assert out[("A", "2")] == pytest.approx(5.0)


def test_misaligned_atoms_and_spread_is_an_error():
    """iterate order and coordinate order diverging is silent corruption."""
    with pytest.raises(PortError, match="orders have diverged"):
        per_residue_spread([Atom("A", "1", "MET", "CA", "", 1.0, 0.0)], [1.0, 2.0])


# -- ensemble_spread_view --------------------------------------------------


ONE_ATOM = [("A", "1", "MET", "CA", "", 1.0, 20.0)]


def test_single_state_object_declines_and_redirects():
    d = spread([[(0.0, 0.0, 0.0)]], ONE_ATOM)
    assert "needs at least 2" in d.report
    assert "altloc_view" in d.report  # points at the right tool for the other case
    assert not d.scene.draws, d.scene
    assert not d.port.queried("spectrum")


def test_spread_view_colours_and_putties():
    d = spread(
        [[(0.0, 0.0, 0.0), (0.0, 0.0, 0.0)], [(2.0, 0.0, 0.0), (0.0, 0.0, 0.0)]],
        ROWS_2,
    )
    assert "2 states, 2 atoms, 2 residues" in d.report
    assert d.scene.has(ColorByScalar) and d.scene.has(SizeByScalar), d.scene
    assert d.port.called("cartoon", "putty", "(obj)"), d.port.call_log
    assert SPREAD_LEGEND in d.report


def test_spread_domain_is_scaled_to_this_object():
    """Unlike occupancy, spread has no absolute scale — so the domain is the
    observed maximum, and the report has to say so or the scaling is implicit."""
    d = spread([[(0.0, 0.0, 0.0)], [(2.0, 0.0, 0.0)]], ONE_ATOM)
    (op,) = d.scene.of(ColorByScalar)
    assert op.domain == (0.0, 1.0), op  # RMS of +/-1 about the centroid
    assert "scaled to this object" in d.report


def test_colour_and_thickness_share_one_field():
    """Two encodings of one quantity must not disagree about its scale."""
    d = spread([[(0.0, 0.0, 0.0)], [(2.0, 0.0, 0.0)]], ONE_ATOM)
    (colour,) = d.scene.of(ColorByScalar)
    (size,) = d.scene.of(SizeByScalar)
    assert colour.field == size.field
    assert colour.domain == size.domain


def test_spread_legend_denies_being_an_error_bar():
    out = spread([[(0.0, 0.0, 0.0)], [(1.0, 0.0, 0.0)]], ONE_ATOM).report
    assert "NOT a calibrated uncertainty" in out
    assert "not mean the true position is less well determined" in out


def test_spread_view_preserves_bfactors():
    d = spread([[(0.0, 0.0, 0.0)], [(1.0, 0.0, 0.0)]], ONE_ATOM)
    assert has_stash("obj")
    assert "restore_bfactors" in d.full_report


def test_putty_can_be_disabled():
    d = spread([[(0.0, 0.0, 0.0)], [(1.0, 0.0, 0.0)]], ONE_ATOM, as_putty=False)
    assert not d.scene.has(SizeByScalar)
    assert not d.port.queried("cartoon")


# -- morph_states — the refusal is the point -------------------------------


def test_morph_refuses_differing_topology_and_explains():
    d = render(morph_states([1, 2], "obj"))
    assert "REFUSED" in d.report
    assert "state 1: 1" in d.report and "state 2: 2" in d.report
    assert "ensemble_spread_view" in d.report  # offers the tool that does work
    assert not d.scene.has(Morph), d.scene
    assert not d.port.queried("morph"), d.port.call_log


def test_morph_refuses_single_state():
    d = render(morph_states([1], "obj"))
    assert "REFUSED" in d.report
    assert not d.port.queried("morph")


def test_morph_proceeds_on_shared_topology():
    d = render(morph_states([1, 1], "obj", steps=10))
    assert "Topology check passed" in d.report
    assert d.port.called("morph", "obj_morph", "obj", refinement=0, steps=10), d.port.call_log


def test_morph_states_says_interpolated_frames_are_invented():
    """SPEC invariant I1 in spirit: generated content must be labelled."""
    report, _ = morph_states([1, 1], "obj")
    assert "generated, not observed" in report


def test_validate_only_creates_nothing():
    d = render(morph_states([1, 1], "obj", validate_only=True))
    assert "no morph created" in d.report
    assert not d.scene.has(Morph)
    assert not d.port.queried("morph")


def test_morph_honours_custom_name():
    d = render(morph_states([1, 1], "obj", name="swing"))
    assert d.port.calls("morph")[0][0][:2] == ("swing", "obj"), d.port.call_log


def test_morph_degrades_when_pymol_is_open_source():
    """cmd.morph is Incentive-only. Most users are on open-source PyMOL, so
    for most users the interpolation half of this tool does not exist. The
    validation half does, and it is the part carrying a judgement."""

    class RefusingPort(FakePort):
        def query(self, action, *args, **kwargs):
            if action == "morph":
                raise PortError(
                    "morph failed: PyMOL execution error: Incentive-Only-Error: "
                    '"morph" is not available in Open-Source PyMOL'
                )
            return super().query(action, *args, **kwargs)

    d = render(morph_states([1, 1], "obj"), port=RefusingPort({}))

    # The judgement is the view's and needs no licence...
    assert "Topology check passed" in d.report
    # ...while the licence caveat is the backend's, because protean
    # interpolates natively and must not inherit PyMOL's advice.
    assert "Incentive-only" not in d.report
    assert "Incentive-only" in d.notes
    assert "set all_states, on" in d.notes  # a way to see the motion regardless
    assert "ensemble_spread_view" in d.notes  # and the tool that needs no licence


def test_other_morph_errors_still_raise():
    """Only the licence case degrades — a real failure must not be swallowed."""

    class BrokenPort(FakePort):
        def query(self, action, *args, **kwargs):
            if action == "morph":
                raise PortError("morph failed: Error: Invalid selection name")
            return super().query(action, *args, **kwargs)

    with pytest.raises(PortError, match="Invalid selection"):
        render(morph_states([1, 1], "obj"), port=BrokenPort({}))
