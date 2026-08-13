"""Tests for local_resolution_view.

Two things carry the judgement in this tool and both are asserted here: the
refusal when the volumes do not share a grid, and the Ångström-to-sigma
conversion of the ramp breakpoints against the *resolution* map's header.

The resolution-map fixture stamps dmin/dmax/dmean/rms = 2/6/4/1, so a
breakpoint in Å converts to sigma by subtracting 4 — 2 Å is -2 sigma, 6 Å is
+2 sigma. The default five-colour palette then lands on exactly
[-2, -1, 0, 1, 2], which makes a wrong conversion impossible to miss.
"""

from __future__ import annotations

import pytest
from conftest import render
from test_mapinfo import write_map

from wiggles_em.backends.pymol import PymolBackend, normalisation_state
from wiggles_em.density import DEFAULT_SIGMA
from wiggles_em.localres import (
    DEFAULT_PALETTE,
    grid_differences,
    local_resolution_view,
)
from wiggles_em.mapinfo import read_map_header
from wiggles_em.maps import forget_map, load_map
from wiggles_em.port import FakePort, PortError
from wiggles_em.provenance import Provenance
from wiggles_em.scene import ColorSurfaceByMap, Isosurface

#: A resolution field runs 2–6 Å. These are resolutions, not densities, which
#: is the whole reason the conversion needs this map's own statistics.
RES_STATS = {"dmin": 2.0, "dmax": 6.0, "dmean": 4.0, "rms": 1.0}


@pytest.fixture(autouse=True)
def _clean_maps():
    forget_map()
    yield
    forget_map()


@pytest.fixture
def session(tmp_path):
    """Two matched volumes in the session: a density map and a resolution map.

    ``get`` is what ``normalize_ccp4_maps`` comes back as. The default "OK" is
    what FakePort returns for an unstubbed command, i.e. PyMOL declining to
    say — which is a path the tool has to handle, so it is the default here.
    """

    def _make(*, main_kw=None, res_kw=None, get="OK"):
        main_path = write_map(tmp_path, "main.mrc", **(main_kw or {}))
        res_path = write_map(tmp_path, "locres.mrc", **{**RES_STATS, **(res_kw or {})})

        loader = FakePort({"get_names": ["main", "locres"]})
        load_map(loader, main_path, "main", provenance=Provenance.MEASURED)
        load_map(loader, res_path, "locres")  # UNKNOWN, the realistic default

        port = FakePort({"get_names": ["main", "locres"], "get": get})
        return port, read_map_header(main_path), read_map_header(res_path)

    return _make


def localres(port, *args, **kwargs):
    """Run local_resolution_view the way a host would.

    The host reads the viewer's normalisation setting and hands the answer to
    the view for its report; the backend reads the same setting for itself when
    it converts. Both go through normalisation_state, so they cannot disagree
    about what PyMOL said — only about what to do with it, which is the point
    of the split.
    """
    # One read, handed to both — which is what a host must do, and what the
    # required `normalised` argument exists to force.
    normalised = kwargs.pop("normalised", None)
    if normalised is None and "normalised" not in kwargs:
        normalised = normalisation_state(port)
    return render(
        local_resolution_view(*args, normalised=normalised, **kwargs),
        port=port,
        normalised=normalised,
    )


# ── the grid check ──────────────────────────────────────────────────────────


def test_identical_grids_have_no_differences(tmp_path):
    a = read_map_header(write_map(tmp_path, "a.mrc"))
    b = read_map_header(write_map(tmp_path, "b.mrc", **RES_STATS))
    assert grid_differences(a, b) == []


def test_voxel_size_is_compared_with_a_tolerance(tmp_path):
    """EMD-30913 reports 0.7999967 where another map reports 0.8. Exact
    equality would refuse a pair that matches."""
    a = read_map_header(write_map(tmp_path, "a.mrc", voxel=0.8))
    b = read_map_header(write_map(tmp_path, "b.mrc", voxel=0.7999967))
    assert grid_differences(a, b) == []


@pytest.mark.parametrize(
    ("kw", "expected"),
    [
        ({"nx": 128, "mx": 256, "cella": (256.0, 256.0, 256.0)}, "extent differs"),
        ({"mx": 512, "cella": (512.0, 256.0, 256.0)}, "grid sampling differs"),
        ({"voxel": 1.2}, "voxel size along X differs"),
        ({"origin": (0.0, 0.0, 12.0)}, "origin along Z differs"),
        ({"nzstart": 5}, "start position differs"),
        ({"mapc": 3, "mapr": 1, "maps": 2}, "axis order differs"),
    ],
)
def test_each_way_the_grids_can_differ_is_named(tmp_path, kw, expected):
    a = read_map_header(write_map(tmp_path, "a.mrc"))
    b = read_map_header(write_map(tmp_path, "b.mrc", **kw))
    differences = grid_differences(a, b)
    assert any(expected in d for d in differences), differences


def test_undefined_voxel_size_is_a_difference_not_a_crash(tmp_path):
    a = read_map_header(write_map(tmp_path, "a.mrc"))
    b = read_map_header(write_map(tmp_path, "b.mrc", cella=(0.0, 256.0, 256.0)))
    assert any("undefined" in d for d in grid_differences(a, b))


def test_mismatched_grids_are_refused_and_nothing_is_drawn(session):
    port, _, _ = session(res_kw={"origin": (0.0, 0.0, 12.0)})
    d = localres(port, "main", "locres")
    out = d.report

    assert "REFUSED" in out
    assert "origin along Z differs" in out
    assert not d.port.calls("isosurface"), d.port.call_log
    assert not d.port.calls("ramp_new"), d.port.call_log


def test_the_refusal_explains_why_a_wrong_grid_is_dangerous(session):
    """A mismatch does not render visibly broken — that is the point."""
    port, _, _ = session(res_kw={"voxel": 1.2})
    d = localres(port, "main", "locres")
    out = d.report
    assert "smooth, plausible, and wrong" in out


def test_the_same_object_twice_is_refused(session):
    port, _, _ = session()
    d = localres(port, "main", "main")
    out = d.report
    assert "REFUSED" in out
    assert "Colouring a map by itself" in out
    assert not d.port.calls("isosurface"), d.port.call_log


# ── the conversion, which is the other half ─────────────────────────────────


def test_default_breakpoints_span_the_headers_range(session):
    port, _, _ = session()
    d = localres(port, "main", "locres")
    out = d.report
    assert "2–6 Å" in out or "2-6 Å" in out


def test_breakpoints_reach_pymol_in_the_resolution_maps_sigma(session):
    """The central claim: PyMOL normalised the volume, so ramp_new needs sigma
    against *this* map's mean and rms, not Ångström and not the main map's."""
    port, _, _ = session()
    d = localres(port, "main", "locres")

    # On the scene the breakpoints are still Angstrom, the unit they were
    # measured in — converting is the backend's job because it is PyMOL's
    # normalisation that makes it necessary.
    (op,) = d.scene.of(ColorSurfaceByMap)
    assert op.volume == "locres"
    assert op.breakpoints == pytest.approx([2.0, 3.0, 4.0, 5.0, 6.0])

    args, _ = d.port.calls("ramp_new")[0]
    assert args[1] == "locres"
    assert args[2] == pytest.approx([-2.0, -1.0, 0.0, 1.0, 2.0]), d.port.call_log
    assert args[3] == list(DEFAULT_PALETTE)


def test_both_units_appear_in_the_report(session):
    """Stating one unit without the other is how the two get confused."""
    port, _, _ = session()
    d = localres(port, "main", "locres")
    out = d.report

    assert "Resolution ramp" in out
    rows = [line for line in out.splitlines() if "->" in line and "Å" in line]
    assert len(rows) == len(DEFAULT_PALETTE), out
    assert "2.00 Å" in rows[0] and "-2" in rows[0], rows[0]
    assert "6.00 Å" in rows[-1] and "2" in rows[-1], rows[-1]


def test_a_pymol_that_will_not_answer_the_setting_query_is_handled(session):
    """Older plugins may not expose `get` at all. An unanswerable question is
    not a failure — it is the unknown case, which is already handled."""

    def _refuse(*args, **kwargs):
        raise PortError("unknown action 'get'")

    port, _, _ = session(get=_refuse)
    d = localres(port, "main", "locres")
    out = d.report

    assert "would not report normalize_ccp4_maps" in out
    assert d.port.calls("ramp_new")[0][0][2] == pytest.approx([-2.0, -1.0, 0.0, 1.0, 2.0])


def test_normalisation_off_sends_angstrom_unconverted(session):
    """With normalize_ccp4_maps off the stored values are still resolutions,
    so converting would be the error."""
    port, _, _ = session(get="off")
    d = localres(port, "main", "locres")

    args, _ = d.port.calls("ramp_new")[0]
    assert args[2] == pytest.approx([2.0, 3.0, 4.0, 5.0, 6.0]), d.port.call_log


def test_unknown_normalisation_assumes_the_pymol_default_and_says_so(session):
    port, _, _ = session()  # FakePort answers "OK", which is not a setting
    d = localres(port, "main", "locres")
    out = d.report

    assert "would not report normalize_ccp4_maps" in out
    assert d.port.calls("ramp_new")[0][0][2] == pytest.approx([-2.0, -1.0, 0.0, 1.0, 2.0])


def test_the_setting_is_read_now_not_at_load_time(session):
    """A limitation worth stating rather than papering over."""
    port, _, _ = session(get="on")
    d = localres(port, "main", "locres")
    out = d.report
    assert "read now, not as it was at load time" in out


def test_zero_rms_refuses_rather_than_dividing(session):
    port, _, _ = session(res_kw={"rms": 0.0})
    d = localres(port, "main", "locres")
    out = d.report

    assert "REFUSED" in out
    assert "sigma is undefined" in out
    assert "normalize_ccp4_maps, off" in out
    assert not d.port.calls("ramp_new"), d.port.call_log


def test_a_map_with_no_positive_values_is_not_a_resolution_field(session):
    port, _, _ = session(res_kw={"dmin": -1.0, "dmax": 0.0, "dmean": -0.5})
    d = localres(port, "main", "locres")
    out = d.report

    assert "REFUSED" in out
    assert "does not look like a resolution field" in out


def test_zero_padding_outside_the_mask_does_not_anchor_the_ramp(session):
    """Estimators write 0 outside the mask. Taking that as the best-resolved
    end would compress every real value into the top of the scale."""
    port, _, _ = session(res_kw={"dmin": 0.0})
    d = localres(port, "main", "locres")
    out = d.report

    assert "outside-the-mask padding" in out
    assert d.port.calls("ramp_new")[0][0][2][0] > -4.0, d.port.call_log


# ── breakpoints and palette given explicitly ────────────────────────────────


def test_explicit_breakpoints_are_used(session):
    port, _, _ = session()
    d = localres(port, "main", "locres", breaks=[3.0, 5.0], palette=["blue", "red"])
    assert d.port.calls("ramp_new")[0][0][2] == pytest.approx([-1.0, 1.0])


def test_descending_breakpoints_are_rejected(session):
    port, _, _ = session()
    with pytest.raises(ValueError, match="must ascend"):
        localres(port, "main", "locres", breaks=[5.0, 3.0], palette=["blue", "red"])


def test_non_positive_breakpoints_are_rejected(session):
    port, _, _ = session()
    with pytest.raises(ValueError, match="must be positive"):
        localres(port, "main", "locres", breaks=[0.0, 3.0], palette=["blue", "red"])


def test_a_single_breakpoint_is_rejected(session):
    port, _, _ = session()
    with pytest.raises(ValueError, match="at least two breakpoints"):
        localres(port, "main", "locres", breaks=[3.0], palette=["blue"])


def test_palette_must_match_the_breakpoints(session):
    port, _, _ = session()
    with pytest.raises(ValueError, match="one to one"):
        localres(port, "main", "locres", breaks=[2.0, 4.0, 6.0])


def test_explicit_breakpoints_raise_rather_than_refusing(session):
    """A bad default is the map's fault and gets a report; a bad argument is
    the caller's and gets an exception."""
    port, _, _ = session(res_kw={"dmax": 0.0, "dmin": -1.0, "dmean": -0.5})
    with pytest.raises(ValueError, match="must be positive"):
        localres(port, "main", "locres", breaks=[-1.0, 2.0], palette=["blue", "red"])


# ── the surface itself ──────────────────────────────────────────────────────


def test_the_surface_is_contoured_in_the_main_maps_sigma(session):
    port, _, _ = session()
    d = localres(port, "main", "locres", level=2.0)

    args, _ = d.port.calls("isosurface")[0]
    assert args == ("main_localres", "main", 2.0), d.port.call_log


def test_an_absolute_contour_is_converted_against_the_main_map(session):
    """main.mrc has dmean=0, rms=0.5, so 1.0 absolute is 2 sigma — and the
    resolution map's statistics must not be the ones used."""
    port, _, _ = session()
    d = localres(port, "main", "locres", level=1.0, units="absolute")
    assert d.port.calls("isosurface")[0][0][2] == pytest.approx(2.0), d.port.call_log


def test_default_contour_is_labelled_generic(session):
    port, _, _ = session()
    d = localres(port, "main", "locres")
    out = d.report

    assert d.port.calls("isosurface")[0][0][2] == DEFAULT_SIGMA
    assert "not a recommendation for this map" in out


def test_bad_units_are_rejected(session):
    port, _, _ = session()
    with pytest.raises(ValueError, match="must be 'sigma' or 'absolute'"):
        localres(port, "main", "locres", level=1.0, units="angstrom")


def test_the_ramp_is_attached_to_the_surface(session):
    port, _, _ = session()
    d = localres(port, "main", "locres")
    # The ramp is a PyMOL object with a PyMOL name; the scene only says
    # "colour this surface by that volume". Mol* has no ramp at all.
    (op,) = d.scene.of(ColorSurfaceByMap)
    assert (op.surface, op.volume) == ("main_localres", "locres")
    assert d.port.called(
        "set", "surface_color", "main_localres_ramp", "main_localres"
    ), d.port.call_log
    assert "deleting it un-colours it" in d.notes


def test_a_selection_carves_the_surface(session):
    port, _, _ = session()
    d = localres(port, "main", "locres", selection="chain A", carve=3.5)

    (op,) = d.scene.of(Isosurface)
    assert op.carve_radius == 3.5
    assert op.carve_around is not None and op.carve_around.dialects == {"pymol"}

    args, kwargs = d.port.calls("isosurface")[0]
    assert args[3] == "(chain A)"
    assert kwargs["carve"] == 3.5


def test_names_can_be_overridden(session):
    port, _, _ = session()
    d = localres(port, "main", "locres", name="surf")

    assert d.port.calls("isosurface")[0][0][0] == "surf"
    # The ramp name follows the surface. It was a parameter when the view
    # spoke PyMOL; now it is the backend's business, because no other viewer
    # has a ramp object to name.
    assert d.port.calls("ramp_new")[0][0][0] == "surf_ramp"
    assert d.port.called("set", "surface_color", "surf_ramp", "surf")


def test_validate_only_creates_nothing(session):
    port, _, _ = session()
    d = localres(port, "main", "locres", validate_only=True)
    out = d.report

    assert "Grid check passed" in out
    assert "nothing created" in out
    assert not d.port.calls("isosurface"), d.port.call_log
    assert not d.port.calls("ramp_new"), d.port.call_log


# ── what the report has to say ──────────────────────────────────────────────


def test_both_maps_get_a_provenance_banner(session):
    """I1: this renders a volume and colours it by a second one, so both
    origins have to be stated."""
    port, _, _ = session()
    d = localres(port, "main", "locres")
    out = d.report

    assert "Provenance: MEASURED" in out
    assert "Provenance: UNKNOWN" in out


def test_the_colour_direction_is_stated_in_words(session):
    """Low Å is good, so the ramp runs opposite to the usual reading."""
    port, _, _ = session()
    d = localres(port, "main", "locres")
    out = d.report

    assert "BEST-resolved" in out
    assert "Low numbers are good" in out


def test_the_two_sigma_scales_are_distinguished(session):
    port, _, _ = session()
    d = localres(port, "main", "locres")
    out = d.report
    assert "not interchangeable" in out


def test_the_legend_says_the_field_is_an_estimate(session):
    port, _, _ = session()
    d = localres(port, "main", "locres")
    out = d.report
    assert "ESTIMATE, not a measurement" in out


def test_geometry_warnings_from_both_maps_are_labelled(session):
    """Both are anisotropic in the same way, so the grids still match and the
    warnings are what is left to report."""
    port, _, _ = session(
        main_kw={"cella": (256.0, 256.0, 384.0)},
        res_kw={"cella": (256.0, 256.0, 384.0)},
    )
    d = localres(port, "main", "locres")
    out = d.report

    assert "Geometry warnings, main" in out
    assert "Geometry warnings, locres" in out
    assert "ANISOTROPIC" in out


# ── refusing to guess ───────────────────────────────────────────────────────


def test_an_unloaded_density_map_is_refused(session):
    port, _, _ = session()
    with pytest.raises(PortError, match="the density map"):
        localres(port, "never_loaded", "locres")


def test_an_unloaded_resolution_map_is_refused(session):
    port, _, _ = session()
    with pytest.raises(PortError, match="the resolution map"):
        localres(port, "main", "never_loaded")


def test_the_refusal_explains_what_the_header_is_needed_for(session):
    port, _, _ = session()
    with pytest.raises(PortError, match="share a grid"):
        localres(port, "main", "never_loaded")


# ── the normalisation answer must be one answer ─────────────────────────────


def test_the_report_and_the_ramp_cannot_disagree_about_units(session):
    """The view's report and the backend's ramp come from one reading.

    Before the seam both were computed inside the view from a single
    normalisation_state() call. Split across two readers they could
    contradict: the colour key listing sigma while the surface was ramped in
    Angstrom, so the user read resolutions off a table that did not describe
    the picture. `normalised` is now a required argument on the view and an
    argument to the backend, so a host reads once and hands the same answer to
    both — which is what `localres` here does.
    """
    port, _, _ = session(get="off")
    d = localres(port, "main", "locres")

    drawn = d.port.calls("ramp_new")[0][0][2]
    # With normalisation off the stored values are still resolutions, so the
    # ramp must carry Angstrom — and the report's table must say the same.
    assert drawn == pytest.approx([2.0, 3.0, 4.0, 5.0, 6.0]), d.port.call_log
    assert "2.00 Å" in d.report

    rows = [line for line in d.report.splitlines() if "->" in line and "Å" in line]
    assert rows, d.report
    for value, row in zip(drawn, rows, strict=True):
        assert f"{value:.6g}" in row, f"ramp says {value}, report row says {row!r}"


def test_a_contour_is_not_converted_when_the_viewer_did_not_normalise(session):
    """The breakpoints honour `normalised`; the contour level must too.

    With normalize_ccp4_maps off, PyMOL reads an isosurface level as an
    absolute map value. Converting 0.05 to 3.16 sigma and sending that gives a
    surface contoured far above dmax — empty — while the report claims 0.05.
    """
    port, _, _ = session(get="off")
    d = localres(port, "main", "locres", level=0.05, units="absolute")
    level = d.port.calls("isosurface")[0][0][2]
    assert level == pytest.approx(0.05), (
        f"sent {level} to an unnormalised map that reads levels as absolute"
    )


class TestAnUnrenderableSceneIsRefusedNotReturned:
    """A view must not hand back a scene the backend will then reject.

    Surviving an unusable rms is right for a normalised session: the caller
    asked for a sigma contour, and only the *absolute equivalent* is unknown.
    But with `normalize_ccp4_maps` off the backend has to convert that contour
    to absolute before it can contour anything, and an unusable rms makes the
    conversion impossible — so the view was returning a full report naming a
    surface that could never be built, and the failure surfaced as a bare
    ValueError from inside the backend after the user had been told it existed.
    """

    def _two_maps(self, tmp_path, main_rms):
        main = write_map(tmp_path, "main.mrc", rms=main_rms)
        res = write_map(tmp_path, "res.mrc", rms=0.5)
        port = FakePort({"get_names": ["main", "res"], "iterate_to_list": []})
        load_map(port, main, "main", provenance=Provenance.MEASURED)
        load_map(port, res, "res", provenance=Provenance.MEASURED)
        return port

    def test_an_unnormalised_session_refuses_when_the_rms_is_unusable(self, tmp_path):
        port = self._two_maps(tmp_path, main_rms=-1.0)

        report, scene = local_resolution_view("main", "res", normalised=False)

        assert "REFUSED" in report, report
        assert not list(scene), "a refusal must draw nothing"
        # And the refusal must survive contact with the backend.
        PymolBackend(port, normalised=False).render(scene)

    def test_the_refusal_explains_which_combination_is_the_problem(self, tmp_path):
        self._two_maps(tmp_path, main_rms=-1.0)

        report, _ = local_resolution_view("main", "res", normalised=False)

        assert "rms" in report.lower()
        assert "normalize_ccp4_maps" in report, (
            "the report should name the session setting that makes this "
            f"impossible, since that is what the user can change:\n{report}"
        )

    def test_a_normalised_session_still_survives_an_unusable_rms(self, tmp_path):
        """The #3 fix must not be undone: with normalisation on, the contour is
        already in the unit PyMOL wants and only the absolute equivalent is
        unknown."""
        port = self._two_maps(tmp_path, main_rms=-1.0)

        report, scene = local_resolution_view("main", "res", normalised=True)

        assert "REFUSED" not in report, report
        assert list(scene), "the view should still draw"
        PymolBackend(port, normalised=True).render(scene)

    def test_an_unnormalised_session_with_a_usable_rms_is_unaffected(self, tmp_path):
        port = self._two_maps(tmp_path, main_rms=0.5)

        report, scene = local_resolution_view("main", "res", normalised=False)

        assert "REFUSED" not in report, report
        PymolBackend(port, normalised=False).render(scene)
