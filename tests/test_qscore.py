"""Tests for Q-score parsing and the view built on it."""

from __future__ import annotations

import gzip

import pytest
from conftest import draw, make_atoms

from wiggles_em.bfactors import has_stash
from wiggles_em.port import PortError
from wiggles_em.qscore import NO_DATA_COLOUR, QSCORE_LEGEND, parse_validation_xml, qscore_view
from wiggles_em.scene import ColorByScalar, ColorFlat

VALIDATION = """<?xml version="1.0" encoding="UTF-8"?>
<wwPDB-validation-information>
  <Entry pdbid="7abc"/>
  <ModelledSubgroup chain="A" resnum="1" resname="MET" Q_score="0.92"/>
  <ModelledSubgroup chain="A" resnum="2" resname="SER" Q_score="0.41"/>
  <ModelledSubgroup chain="B" resnum="1" resname="GLY" Q_score="0.77"/>
  <ModelledSubgroup chain="A" resnum="3" resname="HOH"/>
</wwPDB-validation-information>
"""

ROWS = [
    ("A", "1", "MET", "CA", "", 1.0, 20.0),
    ("A", "2", "SER", "CA", "", 1.0, 30.0),
    ("B", "1", "GLY", "CA", "", 1.0, 25.0),
]


@pytest.fixture
def validation(tmp_path):
    p = tmp_path / "7abc_validation.xml"
    p.write_text(VALIDATION)
    return p


# -- parsing ---------------------------------------------------------------


def test_parses_per_residue_scores(validation):
    scores = parse_validation_xml(validation)
    assert scores[("A", "1")] == pytest.approx(0.92)
    assert scores[("A", "2")] == pytest.approx(0.41)
    assert scores[("B", "1")] == pytest.approx(0.77)


def test_residue_without_a_score_is_omitted_not_defaulted(validation):
    """Absence must stay distinguishable from a low score."""
    scores = parse_validation_xml(validation)
    assert ("A", "3") not in scores
    assert len(scores) == 3


def test_gzipped_validation_file(tmp_path):
    p = tmp_path / "v.xml.gz"
    p.write_bytes(gzip.compress(VALIDATION.encode()))
    assert len(parse_validation_xml(p)) == 3


def test_alternate_attribute_spellings(tmp_path):
    p = tmp_path / "v.xml"
    p.write_text('<r><ModelledSubgroup chain="A" resnum="1" q_score="0.5"/></r>')
    assert parse_validation_xml(p)[("A", "1")] == pytest.approx(0.5)


def test_xray_validation_file_gives_a_useful_error(tmp_path):
    """No Q-score is the normal case for non-3DEM entries; say why."""
    p = tmp_path / "v.xml"
    p.write_text('<r><ModelledSubgroup chain="A" resnum="1" resname="MET"/></r>')
    with pytest.raises(ValueError, match="3DEM"):
        parse_validation_xml(p)


def test_malformed_xml_raises(tmp_path):
    p = tmp_path / "v.xml"
    p.write_text("<not closed")
    with pytest.raises(ValueError, match="not valid XML"):
        parse_validation_xml(p)


# -- the view --------------------------------------------------------------


def test_colours_scored_residues_and_reports_stats(validation):
    d = draw(qscore_view, ROWS, "obj", validation)
    assert "3 of 3 residues scored" in d.report
    assert "0.41" in d.report and "0.92" in d.report
    (op,) = d.scene.of(ColorByScalar)
    # Fixed 0-1: Q-score is a similarity, so autoscaling would make a
    # uniformly poor map look well resolved.
    assert op.domain == (0.0, 1.0), op
    assert d.port.calls("spectrum")[0][1] == {"minimum": 0.0, "maximum": 1.0}
    assert QSCORE_LEGEND in d.report


def test_unscored_residues_are_greyed_not_zeroed(validation):
    """The correctness point: absent data is not a bad score."""
    rows = [*ROWS, ("A", "99", "HOH", "O", "", 1.0, 40.0)]
    d = draw(qscore_view, rows, "obj", validation)

    assert "1 unscored" in d.report
    assert f"{NO_DATA_COLOUR}" in d.report
    assert "deliberately off the scale rather than scored 0" in d.report

    # The unscored residue is coloured flat, never put on the scale...
    (flat,) = d.scene.of(ColorFlat)
    assert flat.colour == NO_DATA_COLOUR
    assert ("A", "99") in {r for s in flat.sel.walk() if s.kind == "residues" for r in s.value}

    # ...and it is absent from the scalar field entirely, which is the stronger
    # claim: padding it with a zero would put absent data on the scale however
    # the selection was written afterwards.
    (scalar,) = d.scene.of(ColorByScalar)
    assert ("A", "99") not in scalar.field.keys
    assert "resi 99" not in d.port.calls("spectrum")[0][0][2]


def test_no_matching_residues_is_a_clear_error(validation):
    rows = [("Z", "500", "MET", "CA", "", 1.0, 20.0)]
    with pytest.raises(PortError, match="author numbering"):
        qscore_view(make_atoms(rows), "obj", validation)


def test_bfactors_preserved_by_default(validation):
    d = draw(qscore_view, ROWS, "obj", validation)
    assert has_stash("obj")
    assert "restore_bfactors" in d.full_report


def test_legend_states_it_is_a_map_model_metric(validation):
    out = draw(qscore_view, ROWS, "obj", validation).report
    assert "map–MODEL agreement metric" in out
    assert "voxel size is wrong" in out  # links back to map_info


def test_worst_residues_are_surfaced(validation):
    out = draw(qscore_view, ROWS, "obj", validation).report
    assert "Least resolvable" in out
    assert "A/2 0.41" in out  # the lowest scorer leads
