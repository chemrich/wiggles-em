"""Every view called the way its own signature says it may be.

The second half of M6, and the gap MCPymol's PR #58 found in its live sweep:
*the sweep supplied its own arguments, so it never called a tool the way its
schema said it could.* ``spheroid`` passed every run on a default that always
failed.

This package has the identical shape. Every other test here passes explicit
arguments, so **no documented default is ever executed** — a default that is
mistyped, that names a palette entry which does not exist, or that fails only
when it is actually reached, would pass the whole suite.

So each view is called with its **required arguments only**, and nothing else.
Optional parameters are deliberately not passed, even where a value would be
more interesting: passing one is what the rest of the suite already does, and
it is precisely what hides this class of bug.

The assertion is deliberately shallow — the view returns a report and a Scene,
and the Scene renders through both backends without raising. Depth belongs in
the per-view modules. What this file proves is that the default *path* is
executed at all, which nothing else here does.
"""

from __future__ import annotations

import pytest

from tests.conftest import make_atoms, render
from tests.test_mapinfo import write_map
from wiggles_em.composition import composition_view
from wiggles_em.deformation import deformation_view
from wiggles_em.density import density_view
from wiggles_em.ensembles import ensemble_spread_view, morph_states
from wiggles_em.heterogeneity import forget_ensemble, load_ensemble
from wiggles_em.latent import latent_traverse_view
from wiggles_em.localres import local_resolution_view
from wiggles_em.maps import forget_map, load_map
from wiggles_em.occupancy import altloc_view, occupancy_view
from wiggles_em.port import FakePort
from wiggles_em.provenance import Provenance
from wiggles_em.qscore import qscore_view
from wiggles_em.scene import Sel

# Two residues, two altlocs, partial occupancy — enough for every atoms-first
# view to have something to say without any of them being a special case.
ROWS = [
    ("A", "1", "MET", "CA", "", 1.0, 20.0),
    ("A", "1", "MET", "CB", "A", 0.6, 22.0),
    ("A", "2", "SER", "CA", "", 1.0, 25.0),
    ("A", "2", "SER", "CB", "B", 0.4, 28.0),
]

COORDS_START = [(0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (2.0, 0.0, 0.0), (3.0, 0.0, 0.0)]
COORDS_END = [(0.5, 0.0, 0.0), (1.5, 0.0, 0.0), (2.5, 0.0, 0.0), (3.5, 0.0, 0.0)]


@pytest.fixture(autouse=True)
def _clean():
    forget_map()
    forget_ensemble()
    yield
    forget_map()
    forget_ensemble()


def _assert_renders(result, *, port=None, normalised=None):
    """A view returned a report and a Scene, and both backends took it."""
    report, scene = result
    assert isinstance(report, str) and report.strip(), "empty report"
    assert list(scene), "the view emitted no ops at all"
    render(result, ROWS, port=port, normalised=normalised)
    return report


# ── atoms-first views ───────────────────────────────────────────────────────


def test_occupancy_view_on_its_defaults():
    _assert_renders(occupancy_view(make_atoms(ROWS), "obj"))


def test_altloc_view_on_its_defaults():
    """`label=True` is the default and draws labels nothing else exercises."""
    _assert_renders(altloc_view(make_atoms(ROWS), "obj"))


def test_qscore_view_on_its_defaults(tmp_path):
    path = tmp_path / "q.xml"
    path.write_text(
        '<?xml version="1.0"?><wwPDB-validation-information>'
        '<ModelledSubgroup chain="A" resnum="1" said="A" Q_score="0.8"/>'
        '<ModelledSubgroup chain="A" resnum="2" said="A" Q_score="0.4"/>'
        "</wwPDB-validation-information>"
    )
    _assert_renders(qscore_view(make_atoms(ROWS), "obj", str(path)))


def test_ensemble_spread_view_on_its_defaults():
    """`as_putty=True` is the default; `superposed` is required and has none,
    which is the point of the `end_state`/`normalised` removals."""
    _assert_renders(
        ensemble_spread_view(
            make_atoms(ROWS),
            [COORDS_START, COORDS_END],
            "obj",
            superposed=True,
        )
    )


def test_deformation_view_on_its_defaults():
    """arrows=True, arrow_scale=1.0, max_arrows=60, as_putty=False,
    uncertainty=None — five defaults, none of them exercised elsewhere."""
    _assert_renders(
        deformation_view(
            make_atoms(ROWS),
            COORDS_START,
            COORDS_END,
            "obj",
            2,
            start_state=1,
            end_state=2,
        )
    )


def test_composition_view_on_its_defaults():
    """transparency=True, label=True, and the ('red', 'skyblue') palette —
    a palette default naming a colour PyMOL does not know would only fail
    here."""
    _assert_renders(
        composition_view(
            "obj",
            {"state_a": 0.7, "state_b": 0.3},
            lambda _sel: 10,
        )
    )


# ── map-backed views ────────────────────────────────────────────────────────


def test_density_view_on_its_defaults(tmp_path):
    """level=None, units='sigma', carve=2.0, name=None. The level default is
    the one that matters: it picks DEFAULT_SIGMA and converts it."""
    path = write_map(tmp_path, "m.mrc")
    port = FakePort({"get_names": ["m"], "iterate_to_list": []})
    load_map(port, path, "m", provenance=Provenance.MEASURED)

    _assert_renders(density_view("m", "polymer"), port=port)


def test_local_resolution_view_on_its_defaults(tmp_path):
    """Everything optional left alone: level, units, breaks, palette,
    selection, carve, name, validate_only. `normalised` is required and
    deliberately has no default — that removal was a fix in its own right."""
    main = write_map(tmp_path, "main.mrc")
    res = write_map(tmp_path, "res.mrc")
    port = FakePort({"get_names": ["main", "res"], "iterate_to_list": []})
    load_map(port, main, "main", provenance=Provenance.MEASURED)
    load_map(port, res, "res", provenance=Provenance.MEASURED)

    _assert_renders(
        local_resolution_view("main", "res", normalised=True), port=port, normalised=True
    )


def test_latent_traverse_view_on_its_defaults(tmp_path):
    """level=None, units='sigma', name=None, color='skyblue',
    build_movie=True. The movie default is what emits Frames at all."""
    root = tmp_path / "ens"
    root.mkdir()
    (root / "z.pkl").write_text("x")
    for i, rms in enumerate((0.5, 0.6, 0.8), start=1):
        write_map(root, f"vol_{i}.mrc", rms=rms)
    names = [f"ens_f{i:02d}" for i in range(1, 4)]
    port = FakePort({"get_names": names, "iterate_to_list": []})
    load_ensemble(port, root, "ens")

    _assert_renders(latent_traverse_view("ens"), port=port)


# ── the non-view entry points that also carry defaults ──────────────────────


def test_morph_states_on_its_defaults():
    """name=None (so `<obj>_morph`), steps=30, validate_only=False. The name
    default is built by string interpolation and is not exercised elsewhere."""
    report = _assert_renders(morph_states([4, 4], "obj"))

    assert "obj_morph" in report, report


def test_composition_view_default_palette_is_two_distinct_colours():
    """The default is a tuple of two; one colour for two states would render
    an indistinguishable picture under a legend claiming otherwise."""
    import inspect

    palette = inspect.signature(composition_view).parameters["palette"].default

    assert len(palette) == 2, palette
    assert palette[0] != palette[1], palette


def test_every_optional_parameter_has_been_left_at_its_default_somewhere():
    """The guard on this file itself.

    If a view grows a new optional parameter, the call above still compiles and
    still passes — the new default simply goes unexercised, which is the exact
    silence this file exists to break. So enumerate the optional parameters of
    every view and require this module to mention each view by name, forcing a
    decision when the signature changes.
    """
    import inspect
    from pathlib import Path

    views = {
        "occupancy_view": occupancy_view,
        "altloc_view": altloc_view,
        "qscore_view": qscore_view,
        "ensemble_spread_view": ensemble_spread_view,
        "deformation_view": deformation_view,
        "composition_view": composition_view,
        "density_view": density_view,
        "local_resolution_view": local_resolution_view,
        "latent_traverse_view": latent_traverse_view,
        "morph_states": morph_states,
    }
    source = Path(__file__).read_text()

    missing = []
    for name, fn in views.items():
        optional = [
            p.name
            for p in inspect.signature(fn).parameters.values()
            if p.default is not inspect.Parameter.empty
        ]
        if optional and f"def test_{name}_on_its_defaults" not in source:
            missing.append((name, optional))

    assert not missing, f"views with defaults but no defaults-only test here: {missing}"


def test_the_selection_helper_is_reachable_from_defaults():
    """`Sel` is what every view lowers through; if the import above ever stops
    resolving, the failures elsewhere would be confusing rather than direct."""
    assert Sel.obj("m").kind == "obj"
