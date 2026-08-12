"""Tests for the sense-1 occupancy views.

These run with no viewer at all. The views are pure — atoms in, a report and a
:class:`~wiggles_em.scene.Scene` out — so most of these assert on the value the
view returned. Where a test asserts on ``port``, it is checking the PyMOL
*lowering*, which is a different claim and worth keeping separate.
"""

from __future__ import annotations

import json

import pytest
from conftest import draw, make_atoms

from wiggles_em.atoms import fetch_atoms
from wiggles_em.bfactors import has_stash, restore_bfactors
from wiggles_em.occupancy import SENSE_1_LEGEND, altloc_view, occupancy_view
from wiggles_em.port import FakePort, PortError
from wiggles_em.scene import ColorByScalar, ColorFlat, Label, Legend, Sense, Show

FULL = [
    ("A", "1", "MET", "CA", "", 1.0, 20.0),
    ("A", "2", "ALA", "CA", "", 1.0, 22.0),
]

PARTIAL = [
    ("A", "1", "MET", "CA", "", 1.0, 20.0),
    ("A", "2", "SER", "CA", "A", 0.6, 25.0),
    ("A", "2", "SER", "CA", "B", 0.4, 31.0),
]


# -- the invariant that matters most --------------------------------------


@pytest.mark.parametrize("view", [occupancy_view, altloc_view], ids=["occupancy", "altloc"])
@pytest.mark.parametrize("rows", [FULL, PARTIAL], ids=["full", "partial"])
def test_every_report_names_the_sense(view, rows):
    """SPEC: occupancy_view and composition_view stay separate forever, and a
    reader must never have to guess which sense a render shows."""
    d = draw(view, rows, "obj")
    assert SENSE_1_LEGEND in d.report
    assert "SENSE 1" in d.report


@pytest.mark.parametrize("view", [occupancy_view, altloc_view], ids=["occupancy", "altloc"])
@pytest.mark.parametrize("rows", [FULL, PARTIAL], ids=["full", "partial"])
def test_the_sense_is_declared_on_the_scene_not_only_in_prose(view, rows):
    """The legend is a value, so the claim can be checked rather than grepped.

    A report is a string and "SENSE 1" appears in it by convention. A
    ``Legend`` carrying ``Sense.ATOM_OCCUPANCY`` is the view *declaring* which
    quantity it drew, which a substring match cannot fake.
    """
    d = draw(view, rows, "obj")
    legends = [op for op in d.ops if isinstance(op, Legend)]
    assert legends, d.scene
    assert all(legend.sense is Sense.ATOM_OCCUPANCY for legend in legends)


@pytest.mark.parametrize("rows", [FULL, PARTIAL], ids=["full", "partial"])
def test_never_claims_compositional_occupancy(rows):
    """No sense-2 language may leak into a sense-1 view."""
    out = draw(occupancy_view, rows, "obj").report.lower()
    # The legend explicitly disclaims these; check they only appear as denial.
    assert "fraction of imaged particles" in out  # the disclaimer
    assert "particles containing" not in out.replace("particles containing this subunit", "")


# -- occupancy_view --------------------------------------------------------


def test_fully_occupied_model_is_not_coloured():
    """A spectrum across a constant is a rainbow that means nothing."""
    d = draw(occupancy_view, FULL, "obj")
    assert "All 2 atoms are fully occupied" in d.report
    assert not d.scene.draws, d.scene
    assert not d.scene.has(ColorByScalar), d.scene
    # And nothing reached PyMOL either — a refusal that still drew would be
    # the failure this asserts against.
    assert not d.port.queried("spectrum"), d.port.call_log
    assert not d.port.queried("alter"), d.port.call_log


def test_partial_occupancy_is_coloured_on_a_fixed_scale():
    d = draw(occupancy_view, PARTIAL, "obj")
    assert "2 of 3 atoms partially occupied" in d.report

    (op,) = d.scene.of(ColorByScalar)
    # Occupancy is absolute, so autoscaling would make a 0.9-occupied model
    # look as partial as a 0.1-occupied one. The domain is the assertion.
    assert op.domain == (0.0, 1.0), op
    assert op.field.values == (1.0, 0.6, 0.4), op.field
    assert d.port.calls("spectrum")[0][1] == {"minimum": 0.0, "maximum": 1.0}, d.port.call_log


def test_the_scalar_field_carries_a_key_per_value():
    """A field with a length mismatch colours by an offset — plausibly, wrongly."""
    (op,) = draw(occupancy_view, PARTIAL, "obj").scene.of(ColorByScalar)
    assert len(op.field.keys) == len(op.field.values) == 3
    assert op.field.keys[1] == ("m", "2")


def test_bfactors_are_preserved_by_default():
    """Assert the values actually come back, not merely that a command ran."""
    d = draw(occupancy_view, PARTIAL, "obj")

    assert has_stash("obj")
    assert "restore_bfactors" in d.full_report

    restore = FakePort()
    restore_bfactors(restore, "obj")
    payload = next(c for c in restore.commands if "stored.wiggles_b" in c)
    saved = json.loads(payload.split("= ", 1)[1])
    # PARTIAL's original B-factors, not the occupancies written over them.
    assert sorted(saved.values()) == [20.0, 25.0, 31.0]


def test_bfactor_destruction_is_warned_when_opted_into():
    d = draw(occupancy_view, PARTIAL, "obj", preserve_bfactors=False)
    assert not has_stash("obj")
    assert "WARNING" in d.notes and "not preserved" in d.notes


def test_the_bfactor_caveat_comes_from_the_backend_not_the_view():
    """It is a fact about PyMOL, not about occupancy.

    protean builds a display copy and re-sends it, so it has nothing to
    restore and must not inherit a warning telling a user to restore it. The
    view's own report therefore says nothing about B-factors at all.
    """
    d = draw(occupancy_view, PARTIAL, "obj")
    assert "B-factor" not in d.report
    assert "B-factor" in d.notes


def test_qfit_floor_is_flagged():
    rows = [("A", "1", "SER", "CA", "A", 0.05, 20.0), ("A", "1", "SER", "CA", "B", 0.95, 20.0)]
    out = draw(occupancy_view, rows, "obj").report
    assert "below q=0.1" in out
    assert "qFit removes" in out


def test_float_noise_counts_as_fully_occupied():
    """Refined occupancies are float32; exact 1.0 comparison mislabels them."""
    rows = [("A", "1", "MET", "CA", "", 0.9999, 20.0)]
    assert "fully occupied" in draw(occupancy_view, rows, "obj").report


# -- altloc_view -----------------------------------------------------------


def test_no_altlocs_says_so_without_overclaiming():
    d = draw(altloc_view, FULL, "obj")
    assert "No alternate conformations" in d.report
    # Must not imply the molecule is rigid — only that this file doesn't model it.
    assert "says nothing about" in d.report
    assert not d.scene.draws, d.scene
    assert not d.port.queried("color"), d.port.call_log


def test_altloc_groups_get_distinct_colours():
    d = draw(altloc_view, PARTIAL, "obj")
    assert "2 altloc group(s)" in d.report

    shown = [op.sel for op in d.scene.of(Show)]
    alts = {s.value for sel in shown for s in sel.walk() if s.key == "alt"}
    assert alts == {"A", "B"}, shown

    used = {op.colour for op in d.scene.of(ColorFlat) if op.colour != "grey70"}
    assert len(used) == 2, f"expected 2 distinct colours, got {used}"


def test_altloc_view_reports_occupancies_not_just_colours():
    """Drawing every alternate at equal weight discards the measurement."""
    out = draw(altloc_view, PARTIAL, "obj").report
    assert "0.6" in out and "0.4" in out
    assert "would discard the measurement" in out


def test_altloc_labels_carry_the_occupancy_as_a_field():
    """The label's value is data, not an embedded PyMOL expression.

    ``"%.2f" % q`` is a Python expression PyMOL evaluates per atom. Mol\\* has
    to build the string itself, so a raw expression would work on exactly one
    backend — the field name travels instead.
    """
    labels = draw(altloc_view, PARTIAL, "obj").scene.of(Label)
    assert labels
    assert all(op.fields == ("q",) for op in labels), labels


def test_labels_can_be_suppressed():
    d = draw(altloc_view, PARTIAL, "obj", label=False)
    assert not d.scene.has(Label)
    assert not d.port.queried("label")


def test_blank_altloc_is_not_a_group():
    """A blank altloc is the absence of grouping, not a group named ''."""
    rows = [
        ("A", "1", "MET", "CA", "", 1.0, 20.0),
        ("A", "2", "SER", "CA", "A", 0.5, 20.0),
        ("A", "2", "SER", "CA", "B", 0.5, 20.0),
    ]
    assert "2 altloc group(s)" in draw(altloc_view, rows, "obj").report


# -- failure behaviour -----------------------------------------------------
#
# These moved down a layer with the seam. Reading atoms is no longer part of a
# view, so an empty or malformed read is fetch_atoms' contract to keep — and
# testing it here would have tested nothing about the view.


def test_empty_selection_is_an_error_not_an_empty_view():
    """MCPymol issue #15: an empty result reported as success is the bug."""
    with pytest.raises(PortError, match="matched no atoms"):
        fetch_atoms(FakePort({"iterate_to_list": []}), "obj")


def test_malformed_atom_row_is_an_error():
    with pytest.raises(PortError, match="malformed atom row"):
        fetch_atoms(FakePort({"iterate_to_list": [("A", "1")]}), "obj")


def test_unanticipated_query_raises_rather_than_returning_none():
    with pytest.raises(KeyError):
        fetch_atoms(FakePort({}), "obj")


def test_a_view_cannot_touch_a_viewer_at_all():
    """The seam, asserted directly: a view given atoms needs no port.

    If a view still reached for one this would raise rather than quietly
    working against whatever was passed.
    """
    report, scene = occupancy_view(make_atoms(PARTIAL), "obj")
    assert report and len(scene)


# -- key identity ----------------------------------------------------------


def test_atoms_differing_only_by_insertion_code_keep_separate_values():
    """The collision bfactors._key documents, asserted end to end.

    Antibody models are full of insertion codes, and the pre-seam code could
    not get this wrong: `alter obj, b=q` let PyMOL read each atom's own q. A
    keyed dictionary can, and a collision does not fail — it colours one
    residue with its neighbour's occupancy on a fixed 0-1 ramp.
    """
    rows = [
        ("A", "100", "SER", "CA", "", 0.30, 20.0, "m", 1),
        ("A", "100A", "SER", "CA", "", 0.90, 21.0, "m", 2),
    ]
    d = draw(occupancy_view, rows, "obj")

    (op,) = d.scene.of(ColorByScalar)
    assert len(set(op.field.keys)) == 2, f"keys collided: {op.field.keys}"

    pushed = next(c for c in d.port.commands if "stored.wiggles_scalar" in c)
    payload = json.loads(pushed.split("= ", 1)[1])
    assert sorted(payload.values()) == [0.30, 0.90], payload


def test_one_selection_spanning_two_models_keeps_them_apart():
    """`obj` can be a named selection covering two loaded structures, where
    chain/resi/name repeat by construction."""
    rows = [
        ("A", "1", "MET", "CA", "", 0.40, 20.0, "first", 1),
        ("A", "1", "MET", "CA", "", 0.80, 20.0, "second", 1),
    ]
    d = draw(occupancy_view, rows, "both")

    (op,) = d.scene.of(ColorByScalar)
    assert len(set(op.field.keys)) == 2, f"keys collided: {op.field.keys}"
