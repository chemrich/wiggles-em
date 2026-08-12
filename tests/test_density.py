"""Tests for density_view — mostly about the sigma/absolute unit trap."""

from __future__ import annotations

import pytest
from conftest import render
from test_mapinfo import write_map

from wiggles_em.density import DEFAULT_SIGMA, density_view, to_absolute, to_sigma
from wiggles_em.mapinfo import read_map_header
from wiggles_em.maps import forget_map, load_map
from wiggles_em.port import FakePort, PortError
from wiggles_em.provenance import Provenance
from wiggles_em.scene import Isosurface, Unit


@pytest.fixture(autouse=True)
def _clean_maps():
    forget_map()
    yield
    forget_map()


@pytest.fixture
def loaded(tmp_path):
    """A map in the session, with known header statistics.

    write_map stamps dmin/dmax/dmean = -1/1/0 and rms = 0.5, so sigma
    conversions here are exact and easy to reason about.
    """

    def _load(filename="m.mrc", **kw):
        path = write_map(tmp_path, filename, **kw)
        obj = filename.split(".")[0]
        port = FakePort({"get_names": [obj]})
        load_map(port, path, obj, provenance=Provenance.MEASURED)
        return FakePort({"get_names": [obj]}), obj, read_map_header(path)

    return _load


# ── the conversion ──────────────────────────────────────────────────────────


def test_sigma_and_absolute_are_inverses(loaded):
    _, _, header = loaded()
    for sigma in (0.0, 1.5, 3.16, -2.0):
        assert to_sigma(header, to_absolute(header, sigma)) == pytest.approx(sigma)


def test_absolute_to_sigma_uses_mean_and_rms(loaded):
    """dmean=0, rms=0.5 in the fixture, so 1.0 absolute is 2 sigma."""
    _, _, header = loaded()
    assert to_sigma(header, 1.0) == pytest.approx(2.0)


def test_zero_rms_cannot_be_converted(tmp_path):
    """A header with rms=0 makes sigma undefined; say so rather than divide."""
    path = write_map(tmp_path, "z.mrc")
    header = read_map_header(path)
    header = type(header)(**{**header.__dict__, "rms": 0.0})
    with pytest.raises(ValueError, match="sigma is undefined"):
        to_sigma(header, 1.0)


def test_real_emdb_contour_converts_to_a_sensible_sigma(loaded):
    """Regression for the trap: EMD-30913 publishes 0.05 absolute, which is
    3.16 sigma against its real header — not 0.05 sigma."""
    _, _, header = loaded()
    header = type(header)(**{**header.__dict__, "dmean": 0.000921512, "rms": 0.0155224})
    assert to_sigma(header, 0.05) == pytest.approx(3.16, abs=0.01)


# ── the view ────────────────────────────────────────────────────────────────


def test_reports_the_level_in_both_units(loaded):
    port, obj, _ = loaded()
    d = render(density_view(obj, "chain A", level=2.0), port=port)

    assert "2 sigma" in d.report
    assert "absolute" in d.report
    assert d.port.calls("isomesh")[0][0][:3] == (f"{obj}_mesh", obj, 2.0), d.port.call_log


def test_absolute_level_is_converted_before_reaching_pymol(loaded):
    """The whole point: PyMOL contours in sigma."""
    port, obj, _ = loaded()
    d = render(density_view(obj, "chain A", level=1.0, units="absolute"), port=port)

    assert d.port.calls("isomesh")[0][0][2] == pytest.approx(2.0), d.port.call_log
    assert "would contour near zero and show mostly noise" in d.report


def test_default_level_is_labelled_as_generic(loaded):
    port, obj, _ = loaded()
    d = render(density_view(obj, "chain A"), port=port)

    assert f"{DEFAULT_SIGMA} sigma was used" in d.report
    assert "not a recommendation for this map" in d.report


def test_emdb_map_points_at_the_author_contour(tmp_path):
    """The finding this tool exists for — and it must say the level is
    absolute, or the advice would cause the very bug it warns about."""
    path = write_map(tmp_path, "emd_30913.mrc")
    port = FakePort({"get_names": ["emd_30913"]})
    load_map(port, path, "emd_30913")

    port = FakePort({"get_names": ["emd_30913"]})
    d = render(density_view("emd_30913", "chain A"), port=port)

    assert "EMD-30913 has an author-recommended contour" in d.report
    assert "ABSOLUTE value" in d.report
    assert "units='absolute'" in d.report
    assert "ebi.ac.uk/emdb/api/entry/EMD-30913" in d.report


def test_non_emdb_map_does_not_invent_an_author_contour(loaded):
    port, obj, _ = loaded("plain.mrc")
    d = render(density_view(obj, "chain A"), port=port)
    assert "author-recommended" not in d.report


def test_carries_the_provenance_banner(loaded):
    """I1: this renders a volume, so the readout must say where it came from."""
    port, obj, _ = loaded()
    d = render(density_view(obj, "chain A"), port=port)
    assert "Provenance: MEASURED" in d.report
    # I1 on the scene, not only in the prose: a view that renders a volume
    # carries a provenance legend as a value, which a substring cannot fake.
    assert any(op.provenance == obj for op in d.scene.legends), d.scene


def test_unloaded_map_is_refused_with_the_reason(loaded):
    with pytest.raises(PortError, match="not loaded through load_map"):
        density_view("never_loaded", "chain A")


def test_refusal_explains_the_unit_hazard(loaded):
    """The error should teach the trap, not just decline."""
    with pytest.raises(PortError, match="author contour gets used as a sigma level"):
        density_view("never_loaded", "chain A")


def test_bad_units_are_rejected(loaded):
    _, obj, _ = loaded()
    with pytest.raises(ValueError, match="must be 'sigma' or 'absolute'"):
        density_view(obj, "chain A", level=1.0, units="angstrom")


def test_carve_radius_is_passed_through(loaded):
    port, obj, _ = loaded()
    d = render(density_view(obj, "chain A", level=1.0, carve=3.5), port=port)
    assert d.port.calls("isomesh")[0][1]["carve"] == 3.5
    assert d.scene.of(Isosurface)[0].carve_radius == 3.5


def test_mesh_name_can_be_overridden(loaded):
    port, obj, _ = loaded()
    d = render(density_view(obj, "chain A", level=1.0, name="pocket"), port=port)
    assert d.port.calls("isomesh")[0][0][0] == "pocket"


def test_geometry_warnings_reach_the_report(loaded):
    port, obj, _ = loaded("aniso.mrc", nx=100, ny=100, nz=100, cella=(100.0, 100.0, 150.0))
    d = render(density_view(obj, "chain A", level=1.0), port=port)
    assert "ANISOTROPIC" in d.report


def test_the_scene_names_the_unit_it_states_the_level_in(loaded):
    """A bare number is the trap. The op carries which unit it means, so a
    backend converts instead of guessing — and Mol*, which takes absolute
    levels natively, is not forced through PyMOL's sigma at all."""
    port, obj, _ = loaded()
    d = render(density_view(obj, "chain A", level=2.0), port=port)
    (op,) = d.scene.of(Isosurface)
    assert op.unit is Unit.SIGMA
    assert op.level == 2.0


# ── a negative RMS is not a usable sigma scale (MCPymol PR #58) ─────────────


def test_a_negative_rms_cannot_define_a_sigma_scale(loaded):
    """MRC2014 writes rms=-1 for "statistics not computed" — what mrcfile
    leaves behind without update_header_stats().

    It divides cleanly, so `if not header.rms` lets it through: to_sigma(0.05)
    returned -2.05, and a resolution ramp turned ascending Angstrom
    breakpoints into a descending one — blue bound to the worst-resolved
    density under a legend saying blue was the best.
    """
    _, obj, _ = loaded("stale.mrc", rms=-1.0)
    with pytest.raises(ValueError, match="sigma is undefined"):
        density_view(obj, "chain A", level=0.05, units="absolute")


def test_to_sigma_and_to_absolute_both_reject_it(tmp_path):
    """to_absolute had no guard at all. `dmean + sigma * -1` returns a number
    of the right shape and the wrong sign, which is the harder failure."""
    header = read_map_header(write_map(tmp_path, "stale.mrc", rms=-1.0))

    with pytest.raises(ValueError, match="never computed"):
        to_sigma(header, 0.05)
    with pytest.raises(ValueError, match="never computed"):
        to_absolute(header, 1.5)


def test_zero_rms_is_still_rejected(tmp_path):
    header = read_map_header(write_map(tmp_path, "flat.mrc", rms=0.0))
    with pytest.raises(ValueError, match="flat map"):
        to_sigma(header, 0.05)
