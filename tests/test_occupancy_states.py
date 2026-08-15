"""Tests for state occupancy — the third sense.

Two claims carry this view, and both are about what it refuses to draw:

1. Occupancy is never encoded as visibility. Fading a rare state in proportion
   to its population would hide precisely the state these methods lose.
2. No colour scale is drawn for weights that are not quantitative, because a
   ramp asserts the differences between the numbers mean something.

Everything else is bookkeeping, and the sense declaration is the guard that
keeps this from being read as either of the other two occupancies.
"""

from __future__ import annotations

import pytest
from conftest import render
from test_mapinfo import write_map

from wiggles_em.heterogeneity import forget_ensemble, load_ensemble
from wiggles_em.occupancy_states import state_occupancy_view
from wiggles_em.populations import Populations, WeightSource
from wiggles_em.port import FakePort, PortError
from wiggles_em.scene import ColorFlat, Isosurface, Legend, Opacity, Sense


@pytest.fixture(autouse=True)
def _clean():
    forget_ensemble()
    yield
    forget_ensemble()


@pytest.fixture
def ensemble(tmp_path):
    def _make(name="ens", *, method_marker="recovar", rms=(0.5, 0.6, 0.8)):
        root = tmp_path / name
        root.mkdir()
        if method_marker:
            (root / method_marker).write_text("x")
        for i, r in enumerate(rms, start=1):
            write_map(root, f"vol_{i}.mrc", rms=r)
        names = [f"{name}_f{i:02d}" for i in range(1, len(rms) + 1)]
        load_ensemble(FakePort({"get_names": names}), root, name)
        return FakePort({"get_names": names}), name

    return _make


def deconvolved(weights, **kw):
    return Populations.declare(weights, WeightSource.DECONVOLVED, **kw)


# ── occupancy is never encoded as visibility ─────────────────────────────────


def test_no_state_is_faded_by_its_population(ensemble):
    """The load-bearing refusal.

    A rare intermediate drawn at 5% opacity is invisible, and that state is the
    one these methods historically lose — the picture would reproduce the
    failure instead of revealing it.
    """
    port, name = ensemble()
    d = render(state_occupancy_view(name, deconvolved([0.9, 0.09, 0.01])), port=port)

    assert not d.scene.of(Opacity), "occupancy must not be drawn as transparency"
    assert not d.port.calls("set_transparency"), d.port.call_log
    assert "EVERY STATE IS DRAWN EQUALLY VISIBLE" in d.report


def test_the_rarest_state_is_still_drawn(ensemble):
    port, name = ensemble()
    d = render(state_occupancy_view(name, deconvolved([0.98, 0.015, 0.005])), port=port)

    assert len(d.scene.of(Isosurface)) == 3, "including the 0.5% one"
    assert "0.50%" in d.report


# ── a ramp only for quantitative weights ─────────────────────────────────────


def test_quantitative_weights_are_coloured_on_a_ramp(ensemble):
    port, name = ensemble()
    d = render(state_occupancy_view(name, deconvolved([0.6, 0.3, 0.1])), port=port)

    colours = [op.colour for op in d.scene.of(ColorFlat)]
    assert len(set(colours)) == 3, "each state's occupancy shows as a distinct colour"


def test_non_quantitative_weights_get_one_flat_colour_and_a_reason(ensemble):
    port, name = ensemble()
    counted = Populations.declare([0.6, 0.3, 0.1], WeightSource.LATENT_HISTOGRAM)
    d = render(state_occupancy_view(name, counted), port=port)

    colours = {op.colour for op in d.scene.of(ColorFlat)}
    assert len(colours) == 1, "a ramp would assert these differences are measurements"
    assert "NO COLOUR SCALE IS DRAWN" in d.report
    assert "not quantitative" in d.report
    # The numbers are still reported — the user asked what the method said.
    assert "60.00%" in d.report


def test_unknown_weights_are_treated_as_non_quantitative(ensemble):
    port, name = ensemble()
    d = render(state_occupancy_view(name, Populations.declare([0.6, 0.3, 0.1])), port=port)

    assert len({op.colour for op in d.scene.of(ColorFlat)}) == 1
    assert "NOT DECLARED" in d.report


# ── the third sense ──────────────────────────────────────────────────────────


def test_every_legend_declares_state_occupancy(ensemble):
    """Without this the render is readable as either of the other two senses."""
    port, name = ensemble()
    d = render(state_occupancy_view(name, deconvolved([0.5, 0.3, 0.2])), port=port)

    declared = [lg.sense for lg in d.scene.legends if lg.sense is not None]
    assert declared, "at least one legend must name the sense"
    assert set(declared) == {Sense.STATE_OCCUPANCY}


def test_the_report_distinguishes_it_from_the_other_two(ensemble):
    port, name = ensemble()
    d = render(state_occupancy_view(name, deconvolved([0.5, 0.5, 0.0])), port=port).report

    assert "sense 3" in d
    assert "not per-atom occupancy" in d.lower()
    assert "100% intact" in d, "the case that makes the distinction concrete"


# ── refusals ─────────────────────────────────────────────────────────────────


def test_a_weight_count_mismatch_is_refused(ensemble):
    _, name = ensemble()
    with pytest.raises(ValueError, match="paired by position"):
        state_occupancy_view(name, deconvolved([0.5, 0.5]))


def test_an_unidentified_method_refuses_to_render(ensemble):
    port, name = ensemble(method_marker=None)
    d = render(state_occupancy_view(name, deconvolved([0.5, 0.3, 0.2])), port=port)

    assert "REFUSED" in d.report
    assert "I2" in d.report
    assert not d.port.calls("isosurface"), d.port.call_log
    assert "stronger claim than a shape" in d.report


def test_an_unloaded_ensemble_is_refused():
    with pytest.raises(PortError, match="was not loaded through load_ensemble"):
        state_occupancy_view("nope", deconvolved([1.0]))


def test_bad_units_are_refused(ensemble):
    _, name = ensemble()
    with pytest.raises(ValueError, match="units must be"):
        state_occupancy_view(name, deconvolved([0.5, 0.3, 0.2]), units="angstrom")


# ── free energy in the report ────────────────────────────────────────────────


def test_free_energy_appears_only_when_a_temperature_was_declared(ensemble):
    port, name = ensemble()
    without = render(state_occupancy_view(name, deconvolved([0.6, 0.3, 0.1])), port=port).report
    assert "free energy" not in without.lower()

    forget_ensemble()
    port, name = ensemble("ens2")
    with_t = render(
        state_occupancy_view(name, deconvolved([0.6, 0.3, 0.1], temperature_k=298.15)),
        port=port,
    ).report
    assert "Relative free energy at 298.15 K" in with_t


def test_an_unbounded_error_is_labelled_not_printed_as_a_number(ensemble):
    """A state whose uncertainty reaches its own population has no upper bound."""
    port, name = ensemble()
    pops = deconvolved([0.9, 0.08, 0.02], uncertainty=[0.01, 0.01, 0.02], temperature_k=298.15)
    d = render(state_occupancy_view(name, pops), port=port).report

    assert "+unbounded" in d
    assert "arbitrarily unfavourable" in d


# ── the contour trap, inherited from latent_traverse_view ────────────────────


def test_one_absolute_level_becomes_a_different_sigma_per_state(ensemble):
    """Same trap as a traversal: PyMOL normalises each map on load, so holding
    sigma constant would contour each state on its own scale."""
    port, name = ensemble(rms=(0.5, 1.0, 2.0))
    d = render(
        state_occupancy_view(name, deconvolved([0.5, 0.3, 0.2]), level=1.0, units="absolute"),
        port=port,
    )

    assert [op.level for op in d.scene.of(Isosurface)] == pytest.approx([2.0, 1.0, 0.5])


def test_a_state_with_no_usable_rms_is_skipped_and_the_total_is_flagged(ensemble):
    port, name = ensemble(rms=(0.5, 0.0, 0.8))
    d = render(state_occupancy_view(name, deconvolved([0.5, 0.3, 0.2])), port=port).report

    assert "1 state(s) not contoured" in d
    assert "not sum to 100%" in d


def test_legends_are_scene_values_not_only_report_text(ensemble):
    """A viewer host lowering the Scene never sees the report string."""
    port, name = ensemble()
    d = render(state_occupancy_view(name, deconvolved([0.5, 0.3, 0.2])), port=port)

    texts = [op.text for op in d.scene.ops if isinstance(op, Legend)]
    assert any("STATE OCCUPANCY (sense 3)" in t for t in texts)
    assert any("EQUALLY VISIBLE" in t for t in texts)
