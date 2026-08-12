"""Tests for latent_traverse_view, and for invariants I2 and I3.

These two invariants are the reason tier 3 is worth building rather than
reading about. I2 says no latent rendering may go out unlabelled; I3 says no
output may turn a gap into an absence. Both are asserted here directly, and I3
is asserted over *every* tier-3 report rather than only this one, because the
invariant is about what this package says, not about one function.

The fixture stamps each frame's rms differently on purpose. That is what makes
the per-frame sigma conversion observable: with identical statistics the whole
mechanism would be invisible and the tests would pass whether it worked or not.
"""

from __future__ import annotations

import pytest
from conftest import render
from test_mapinfo import write_map

from wiggles_em.heterogeneity import forget_ensemble, load_ensemble, loaded_ensemble
from wiggles_em.latent import (
    ABSENCE_CLAIMS,
    contains_absence_claim,
    frame_levels,
    latent_traverse_view,
)
from wiggles_em.port import FakePort, PortError
from wiggles_em.scene import Frames, Isosurface, Unit


@pytest.fixture(autouse=True)
def _clean():
    forget_ensemble()
    yield
    forget_ensemble()


@pytest.fixture
def ensemble(tmp_path):
    """An ensemble whose frames have deliberately different density statistics."""

    def _make(name="ens", *, method_marker="z.pkl", rms=(0.5, 0.6, 0.8), dmean=0.0):
        root = tmp_path / name
        root.mkdir()
        if method_marker:
            (root / method_marker).write_text("x")
        for i, r in enumerate(rms, start=1):
            write_map(root, f"vol_{i}.mrc", rms=r, dmean=dmean)
        names = [f"{name}_f{i:02d}" for i in range(1, len(rms) + 1)]
        load_ensemble(FakePort({"get_names": names}), root, name)
        return FakePort({"get_names": names}), name

    return _make


# ── I2: no unlabelled latent plot ────────────────────────────────────────────


def test_an_unidentified_method_refuses_to_render(ensemble):
    """The invariant working, not the tool failing."""
    port, name = ensemble(method_marker=None)
    d = render(latent_traverse_view(name), port=port)
    out = d.report

    assert "REFUSED" in out
    assert "I2" in out
    assert not d.port.calls("isosurface"), d.port.call_log


def test_the_refusal_explains_why_the_caveats_are_not_interchangeable(ensemble):
    port, name = ensemble(method_marker=None)
    d = render(latent_traverse_view(name), port=port)
    out = d.report

    assert "cryoDRGN" in out and "3DVA" in out
    assert "looks equally convincing either way" in out


def test_the_refusal_names_the_way_forward(ensemble):
    port, name = ensemble(method_marker=None)
    d = render(latent_traverse_view(name), port=port)
    out = d.report

    assert "method='cryodrgn'" in out
    assert "dynamight" in out


def test_an_identified_method_renders_and_carries_its_caveat(ensemble):
    port, name = ensemble()
    d = render(latent_traverse_view(name), port=port)
    out = d.report

    assert "REFUSED" not in out
    assert "cryoDRGN" in out
    assert "NOT RELIABLY MEANINGFUL" in out, "the method's own caveat must be present"
    assert len(d.port.calls("isosurface")) == 3, d.port.call_log


def test_every_method_has_a_distinct_caveat(ensemble):
    """I2 is worthless if the caveats are boilerplate."""
    from wiggles_em.heterogeneity import Method

    caveats = {m: m.caveat for m in Method}
    assert len(set(caveats.values())) == len(caveats)
    for method, caveat in caveats.items():
        assert len(caveat) > 80, f"{method} has a stub caveat"


# ── I3: a gap is not an absence ──────────────────────────────────────────────


def test_the_gap_legend_is_present(ensemble):
    port, name = ensemble()
    d = render(latent_traverse_view(name), port=port)
    out = d.report

    assert "A GAP IS NOT AN ABSENCE" in out
    assert "UNKNOWN OCCUPANCY" in out


def test_no_population_is_drawn_and_the_report_says_why(ensemble):
    port, name = ensemble()
    d = render(latent_traverse_view(name), port=port)
    out = d.report

    assert "NO POPULATION IS SHOWN" in out
    assert not d.port.calls("ramp_new"), "a density ramp would be a population claim"


def test_the_absence_detector_catches_its_phrases():
    for phrase in ABSENCE_CLAIMS:
        assert contains_absence_claim(f"this region is {phrase} by the molecule") == phrase
    assert contains_absence_claim("these conformations are supported") is None


def test_the_detector_is_case_insensitive():
    assert contains_absence_claim("A FORBIDDEN STATE") == "forbidden state"


# ── the per-frame normalisation trap ─────────────────────────────────────────


def test_one_absolute_level_becomes_a_different_sigma_per_frame(ensemble):
    """The whole point. PyMOL normalises each map independently, so holding a
    sigma level constant would contour each frame on its own scale and rescale
    away the density change the traversal exists to show."""
    port, name = ensemble(rms=(0.5, 1.0, 2.0))
    d = render(latent_traverse_view(name, level=1.0, units="absolute"), port=port)

    sigmas = [args[2] for args, _ in d.port.calls("isosurface")]
    assert sigmas == pytest.approx([2.0, 1.0, 0.5]), d.port.call_log

    # On the scene, each frame carries its own level and names its unit, so a
    # backend that contours in absolute values is handed the right thing too.
    ops = d.scene.of(Isosurface)
    assert [op.level for op in ops] == pytest.approx([2.0, 1.0, 0.5])
    assert all(op.unit is Unit.SIGMA for op in ops)


def test_a_sigma_level_is_anchored_to_the_first_frame(ensemble):
    port, name = ensemble(rms=(0.5, 1.0))
    d = render(latent_traverse_view(name, level=2.0, units="sigma"), port=port)

    sigmas = [args[2] for args, _ in d.port.calls("isosurface")]
    assert sigmas[0] == pytest.approx(2.0)
    assert sigmas[1] == pytest.approx(1.0), "same absolute value, this frame's sigma"


def test_a_varying_sigma_spread_is_called_out_as_the_signal(ensemble):
    port, name = ensemble(rms=(0.5, 1.0, 2.0))
    d = render(latent_traverse_view(name), port=port)
    out = d.report

    assert "IS the density change" in out
    assert "flattened exactly this signal" in out


def test_near_identical_frames_are_described_as_such(ensemble):
    port, name = ensemble(rms=(0.5, 0.5, 0.5))
    d = render(latent_traverse_view(name), port=port)
    out = d.report

    assert "barely differ" in out


def test_frame_levels_reports_none_for_an_undefined_sigma(ensemble):
    _, name = ensemble(rms=(0.5, 0.0))
    levels = frame_levels(loaded_ensemble(name), 1.0)

    assert levels[0] is not None
    assert levels[1] is None


def test_a_frame_with_zero_rms_is_skipped_and_reported(ensemble):
    port, name = ensemble(rms=(0.5, 0.0, 0.8))
    d = render(latent_traverse_view(name), port=port)
    out = d.report

    assert len(d.port.calls("isosurface")) == 2, d.port.call_log
    assert "1 frame(s) were skipped" in out
    assert "rms=0" in out


def test_all_frames_unusable_is_refused(ensemble):
    port, name = ensemble(rms=(0.0, 0.0))
    with pytest.raises(PortError, match=r"no frame .* has a usable RMS"):
        render(latent_traverse_view(name), port=port)


# ── the traversal itself ─────────────────────────────────────────────────────


def test_frames_are_wired_to_the_movie_timeline(ensemble):
    port, name = ensemble()
    d = render(latent_traverse_view(name), port=port)

    (frames,) = d.scene.of(Frames)
    assert frames.build_timeline
    assert len(frames.names) == 3
    assert d.port.called("mset", "1 x3"), d.port.call_log
    assert len(d.port.calls("mdo")) == 3
    first_frame, command = d.port.calls("mdo")[0][0]
    assert first_frame == 1
    assert "enable ens_surf_01" in command


def test_the_movie_can_be_switched_off(ensemble):
    port, name = ensemble()
    d = render(latent_traverse_view(name, build_movie=False), port=port)

    # The frames still exist as a sequence; only the viewer playback is off.
    (frames,) = d.scene.of(Frames)
    assert not frames.build_timeline
    assert not d.port.calls("mset"), d.port.call_log
    assert len(d.port.calls("isosurface")) == 3


def test_one_colour_for_every_frame(ensemble):
    """A per-frame spectrum would encode frame index as if it were measured."""
    port, name = ensemble()
    d = render(latent_traverse_view(name, color="yellow"), port=port)

    colours = {args[0] for args, _ in d.port.calls("color")}
    assert colours == {"yellow"}


def test_surface_names_can_be_overridden(ensemble):
    port, name = ensemble()
    d = render(latent_traverse_view(name, name="traj"), port=port)
    assert d.port.calls("isosurface")[0][0][0] == "traj_01"


def test_an_unloaded_ensemble_is_refused_with_the_reason(ensemble):
    with pytest.raises(PortError, match="not loaded through load_ensemble"):
        latent_traverse_view("never_loaded")


def test_bad_units_are_rejected(ensemble):
    port, name = ensemble()
    with pytest.raises(ValueError, match="must be 'sigma' or 'absolute'"):
        render(latent_traverse_view(name, level=1.0, units="angstrom"), port=port)


def test_the_provenance_banner_is_carried(ensemble):
    port, name = ensemble()
    d = render(latent_traverse_view(name), port=port)
    out = d.report
    assert "Provenance: GENERATED" in out
