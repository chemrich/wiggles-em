"""Tests for composition_view — sense 2, and the line it must never cross.

The compendium's central design decision is that `occupancy_view` (sense 1,
per-atom `q`) and `composition_view` (sense 2, fraction of particles) stay
separately named forever and neither infers the other. The most important test
in this file is therefore a negative one: this tool must never read an atom
property at all. If it ever did, the render would look identical and the number
would mean something else.
"""

from __future__ import annotations

import pytest
from conftest import render

from wiggles_em.composition import composition_view, parse_composition_table
from wiggles_em.port import ITERATE_TO_LIST, PortError
from wiggles_em.scene import ColorFlat, Label, Opacity, Sense


def _counts(table, **overrides):
    """Atom counts per table selection, as the host would have counted them.

    Counting moved out of the view with the Scene seam: it is a viewer read,
    and the view's job is to refuse when a count is zero rather than to go and
    get it. The refusal is what these tests are about, and it is unchanged.
    """
    return {sel: overrides.get(sel, 10) for sel in table}


def _view(table, **kwargs):
    counts = kwargs.pop("counts", None)
    return render(composition_view(counts or _counts(table), "obj", table, **kwargs))


# ── the table ────────────────────────────────────────────────────────────────


def test_a_dict_is_taken_as_given():
    assert parse_composition_table({"chain A": 0.4}) == {"chain A": 0.4}


def test_an_inline_string_is_parsed():
    got = parse_composition_table("chain A=0.4, chain B=1.0")
    assert got == {"chain A": 0.4, "chain B": 1.0}


def test_percentages_are_accepted():
    assert parse_composition_table("chain A=40%") == {"chain A": 0.4}


def test_a_file_of_tab_separated_rows_is_read(tmp_path):
    path = tmp_path / "comp.tsv"
    path.write_text("# part\tfraction\nchain A\t1.0\nchain B and resi 1-40\t0.38\n")

    got = parse_composition_table(path)
    assert got == {"chain A": 1.0, "chain B and resi 1-40": 0.38}


def test_selections_may_contain_spaces(tmp_path):
    """Which is why the delimiter is '=' or a tab, never whitespace."""
    got = parse_composition_table("chain B and resi 1-40=0.5")
    assert got == {"chain B and resi 1-40": 0.5}


def test_an_empty_table_is_refused():
    with pytest.raises(ValueError, match="no entries found"):
        parse_composition_table("# just a comment")


@pytest.mark.parametrize("value", [1.4, -0.2, 40.0])
def test_a_fraction_outside_zero_to_one_is_refused(value):
    with pytest.raises(ValueError, match=r"must lie in \[0, 1\]"):
        parse_composition_table({"chain A": value})


def test_the_refusal_guesses_the_likely_mistake():
    with pytest.raises(ValueError, match="raw particle count or a percentage"):
        parse_composition_table({"chain A": 40.0})


# ── the line this tool must not cross ────────────────────────────────────────


def test_no_atom_property_is_ever_read():
    """The whole design decision, as an assertion. Deriving sense 2 from sense 1
    would present a number that means one thing as though it meant another, and
    the picture would be identical."""
    d = _view({"chain A": 0.4})

    assert not d.port.queried(ITERATE_TO_LIST), d.port.call_log
    assert not any("q" in str(args) for args, _ in d.port.calls("alter")), d.port.call_log
    # And on the scene: nothing here carries a scalar field at all, so there is
    # no per-atom quantity for sense 1 to leak in through.
    assert not d.fake.fields, d.scene


def test_the_legend_names_sense_2_and_disclaims_q():
    d = _view({"chain A": 0.4})
    out = d.report

    assert "SENSE 2" in out
    assert "NOT the per-atom crystallographic occupancy q" in out
    assert "NOT derived from q" in out
    assert "occupancy_view" in out, "the reader must be pointed at sense 1"
    # Declared on the scene as well, so the two senses can never be told apart
    # by prose alone.
    assert all(
        legend.sense is Sense.PARTICLE_COMPOSITION for legend in d.scene.legends
    ), d.scene


def test_the_estimate_is_not_presented_as_measured():
    out = _view({"chain A": 0.4}).report
    assert "came from the table you supplied" in out
    assert "carries the uncertainty of whichever method made it" in out


# ── the empty-selection trap ─────────────────────────────────────────────────


def test_a_selection_matching_nothing_is_refused():
    """PyMOL would accept it, colour nothing, and leave a render that looks like
    a fully-present structure."""
    with pytest.raises(PortError, match="matched no atoms"):
        composition_view({"chain Z": 0}, "obj", {"chain Z": 0.5})


def test_the_refusal_explains_what_the_render_would_have_looked_like():
    with pytest.raises(PortError, match="looks like a fully-present structure"):
        composition_view({"chain Z": 0}, "obj", {"chain Z": 0.5})


def test_nothing_is_drawn_when_a_selection_is_empty():
    table = {"chain A": 1.0, "chain Z": 0.5}
    with pytest.raises(PortError):
        composition_view(_counts(table, **{"chain Z": 0}), "obj", table)


# ── the view ─────────────────────────────────────────────────────────────────


def test_each_part_is_coloured():
    d = _view({"chain A": 1.0, "chain B": 0.3})

    coloured = [args[1] for args, _ in d.port.calls("color")]
    assert "((obj) and (chain A))" in coloured
    assert "((obj) and (chain B))" in coloured
    # The table's text is the caller's, in the viewer's own dialect, so it
    # travels as raw and a backend that cannot parse it refuses.
    assert all(op.sel.dialects == {"pymol"} for op in d.scene.of(ColorFlat))


def test_transparency_is_the_complement_of_presence():
    """A part in 40% of particles is drawn 60% transparent."""
    d = _view({"chain A": 0.4})

    (op,) = d.scene.of(Opacity)
    # Stated as opacity on the scene — transparency is its inverse, and every
    # viewer picks a different one of the two.
    assert op.value == pytest.approx(0.4)
    settings = {(args[0], args[2]): args[1] for args, _ in d.port.calls("set")}
    assert settings[("transparency", "((obj) and (chain A))")] == pytest.approx(0.6)


def test_transparency_can_be_switched_off():
    d = _view({"chain A": 0.4}, transparency=False)
    assert not d.scene.has(Opacity)
    assert not d.port.calls("set"), d.port.call_log


def test_labels_carry_the_percentage():
    d = _view({"chain A": 0.38})
    assert any("38%" in op.text for op in d.scene.of(Label)), d.scene
    assert any("38%" in str(args) for args, _ in d.port.calls("label")), d.port.call_log


def test_labels_can_be_switched_off():
    d = _view({"chain A": 0.4}, label=False)
    assert not d.scene.has(Label)
    assert not d.port.calls("label"), d.port.call_log


def test_the_report_orders_parts_by_presence():
    out = _view({"chain A": 1.0, "chain B": 0.2}).report
    assert out.index("chain B") < out.index("chain A"), "least present first"
    assert "Least present: chain B at 20%" in out


def test_the_colour_direction_is_stated():
    out = _view({"chain A": 0.4}).report
    assert "red (rarely present) → blue (always present)" in out
    assert "opposite to the usual" in out


def test_atom_counts_reach_the_report():
    table = {"chain A": 0.4}
    out = _view(table, counts=_counts(table, **{"chain A": 137})).report
    assert "137 atoms" in out
