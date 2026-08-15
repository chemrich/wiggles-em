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
    FRAME_QUALITY_LEGEND,
    contains_absence_claim,
    frame_levels,
    latent_traverse_view,
)
from wiggles_em.port import FakePort, PortError
from wiggles_em.scene import Frames, Isosurface, Legend, Unit


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


# ── frame quality: the thing this viewer cannot measure ──────────────────────


def test_the_frame_quality_caveat_reaches_the_picture_and_the_report(ensemble):
    """Both surfaces, because they are read by different people.

    A caller sees the report; a viewer host lowering the Scene sees only the
    ops. A caveat present in one and not the other is absent for half the
    audience.
    """
    port, name = ensemble()
    d = render(latent_traverse_view(name), port=port)

    assert "NOT NECESSARILY OF EQUAL QUALITY" in d.report
    legends = [op.text for op in d.scene.ops if isinstance(op, Legend)]
    assert FRAME_QUALITY_LEGEND in legends


def test_the_frame_quality_caveat_does_not_vary_with_the_data(ensemble):
    """It states a limit of the viewer, not a finding about this ensemble.

    Resolution is not in a map header and the FSC that would measure it needs a
    reference no session has, so any per-ensemble wording here would be a
    measurement we did not make. Two ensembles with deliberately different
    density statistics must therefore get the identical sentence — if this ever
    becomes conditional on the frames, that is the defect, and this is what
    catches it.
    """
    port_tight, tight = ensemble("tight", rms=(0.5, 0.5, 0.5))
    port_spread, spread = ensemble("spread", rms=(0.2, 0.9, 2.4))

    tight_report = render(latent_traverse_view(tight), port=port_tight).report
    spread_report = render(latent_traverse_view(spread), port=port_spread).report

    # The sigma spread genuinely differs between them — otherwise this test
    # would pass on two ensembles the tool could not tell apart anyway.
    assert "differ by" in spread_report
    assert "barely differ" in tight_report

    assert FRAME_QUALITY_LEGEND in tight_report
    assert FRAME_QUALITY_LEGEND in spread_report


def test_the_frame_quality_caveat_makes_no_absence_claim():
    """It is new text under an existing invariant, so I3 must still hold."""
    assert contains_absence_claim(FRAME_QUALITY_LEGEND) is None


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

    # The view resolves "yellow" to RGB before it reaches the Scene, so PyMOL
    # is asked to define a colour for that triple and then to use it. The claim
    # is unchanged, and now pins the colour itself rather than a name PyMOL
    # would have had to interpret.
    defined = {tuple(args[1]) for args, _ in d.port.calls("set_color")}
    assert defined == {(1.0, 1.0, 0.0)}, defined

    colours = {args[0] for args, _ in d.port.calls("color")}
    assert len(colours) == 1, f"one colour for every frame, got {colours}"


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


# ── surface names must name the frame they came from ────────────────────────


def test_a_skipped_frame_does_not_renumber_the_ones_after_it(ensemble):
    """Frame 3 is unusable, so the surfaces are 1, 2, 4, 5 — never 1, 2, 3, 4.

    Numbering over the surviving frames makes `_03` hold frame 4's density.
    Nothing in the report says which frame was dropped, so the user reads each
    density against the wrong latent coordinate — the off-by-one that
    heterogeneity._natural_key exists to prevent, reintroduced downstream.
    """
    port, name = ensemble(rms=(0.5, 0.6, 0.0, 0.8, 1.0))
    d = render(latent_traverse_view(name), port=port)

    made = [args[0] for args, _ in d.port.calls("isosurface")]
    assert made == [f"{name}_surf_{i:02d}" for i in (1, 2, 4, 5)], made

    # And the source volume must match the number in the name.
    volumes = [args[1] for args, _ in d.port.calls("isosurface")]
    for surface, volume in zip(made, volumes, strict=True):
        assert surface.rsplit("_", 1)[1].lstrip("0") == volume.rsplit("f", 1)[1].lstrip("0")

    assert "frame 3" in d.report, "the report must name which frame was skipped"


def test_the_anchor_skips_frames_whose_rms_cannot_define_sigma(ensemble):
    """A sigma level is interpreted against the first *usable* frame.

    Anchoring on frame 0 regardless meant a leading rms=0 or rms=-1 header
    yielded an anchor of dmean — a number unrelated to the level asked for —
    which was then applied to every other frame.
    """
    port, name = ensemble(rms=(0.0, 0.5, 1.0))
    d = render(latent_traverse_view(name, level=2.0, units="sigma"), port=port)

    # Frame 1 is unusable, so the anchor is frame 2: 2 sigma there.
    sigmas = [args[2] for args, _ in d.port.calls("isosurface")]
    assert sigmas[0] == pytest.approx(2.0), d.port.call_log
    assert sigmas[1] == pytest.approx(1.0), "same absolute value, frame 3's sigma"


def test_no_frame_with_a_usable_rms_is_refused(ensemble):
    _, name = ensemble(rms=(-1.0, 0.0))
    with pytest.raises(PortError, match=r"no frame .* has a usable RMS"):
        latent_traverse_view(name, level=2.0, units="sigma")


class TestTheAnchorFrameIsAlwaysAttributed:
    """A sigma level is meaningless without the header it was read against.

    The anchor is the first frame whose header carries a usable rms, which is
    frame 1 only when frame 1's statistics exist. The report used to claim "the
    FIRST frame" unconditionally; the fix for that then attributed the anchor
    only on the defaulted path, leaving a caller-supplied sigma silently
    interpreted against a frame the caller never named — the same defect, one
    branch over.
    """

    def _gapped(self, tmp_path):
        """Frame 1 has rms=0, so the anchor is frame 2."""
        root = tmp_path / "ens"
        root.mkdir()
        (root / "z.pkl").write_text("x")
        write_map(root, "vol_1.mrc", rms=0.0, dmean=10.0)
        write_map(root, "vol_2.mrc", rms=0.5, dmean=10.0)
        write_map(root, "vol_3.mrc", rms=0.6, dmean=10.0)
        names = [f"ens_f{i:02d}" for i in range(1, 4)]
        load_ensemble(FakePort({"get_names": names, "iterate_to_list": []}), root, "ens")

    def test_a_defaulted_level_names_its_anchor(self, tmp_path):
        self._gapped(tmp_path)
        report, _ = latent_traverse_view("ens")

        assert "frame 2" in report, report
        assert "FIRST frame" not in report

    def test_an_explicit_sigma_level_names_its_anchor(self, tmp_path):
        """The branch the previous fix missed."""
        self._gapped(tmp_path)
        report, _ = latent_traverse_view("ens", level=1.5, units="sigma")

        assert "frame 2" in report, (
            "a caller-supplied sigma was interpreted against frame 2 and the "
            f"report never says so:\n{report}"
        )

    def test_an_absolute_level_claims_no_anchor(self, tmp_path):
        """An absolute level needs no anchor, so attributing one would be a
        claim about a conversion that never happened."""
        self._gapped(tmp_path)
        report, _ = latent_traverse_view("ens", level=10.75, units="absolute")

        assert "interpreted against frame" not in report, report
        assert "sigma against frame" not in report, report

    def test_the_anchor_is_not_explained_away_when_it_is_frame_one(self, tmp_path):
        """The "not frame 1" paragraph must not appear when it *is* frame 1."""
        root = tmp_path / "ok"
        root.mkdir()
        (root / "z.pkl").write_text("x")
        for i, rms in enumerate((0.5, 0.6, 0.8), start=1):
            write_map(root, f"vol_{i}.mrc", rms=rms)
        names = [f"ok_f{i:02d}" for i in range(1, 4)]
        load_ensemble(FakePort({"get_names": names, "iterate_to_list": []}), root, "ok")

        report, _ = latent_traverse_view("ok", level=1.5, units="sigma")

        assert "frame 1" in report, report
        assert "not frame 1" not in report, report


def test_a_frame_skipped_for_rms_minus_one_is_not_described_as_flat(tmp_path):
    """MRC writes rms=-1 for "statistics never computed" — mrcfile's default
    without update_header_stats(). That is a different problem from a flat map,
    with a different remedy, and the report hard-coded the flat one.

    The same defect was fixed in localres.py; this site was missed because the
    completeness grep was for the predicate `header.rms` and not for the string.
    """
    root = tmp_path / "ens"
    root.mkdir()
    (root / "z.pkl").write_text("x")
    write_map(root, "vol_1.mrc", rms=0.5)
    write_map(root, "vol_2.mrc", rms=-1.0)  # statistics never computed
    write_map(root, "vol_3.mrc", rms=0.6)
    names = [f"ens_f{i:02d}" for i in range(1, 4)]
    load_ensemble(FakePort({"get_names": names, "iterate_to_list": []}), root, "ens")

    report, _ = latent_traverse_view("ens")

    assert "frame 2" in report, report
    assert "rms=0" not in report, (
        "a frame whose statistics were never computed is reported as having "
        f"rms=0, which says the map is flat:\n{report}"
    )
    assert "never computed" in report, (
        f"the report should say what rms=-1 actually means:\n{report}"
    )
