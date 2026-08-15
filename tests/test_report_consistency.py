"""Do the report's factual claims match the ops the view actually emitted?

M3 of FIX-PLAN.md. Four of the ten findings in REVIEW-2026-08-13.md are one
shape: **the report asserts something the scene or the backend does not do.**

  #5  the report promises "a surface's number is always the frame it was made
      from"; the timeline renumbers over the survivors
  #8  the report says the level was taken against the FIRST frame; the code
      anchors on the first frame with a usable RMS
  #6  the ramp table states the breakpoints were sent unconverted; the backend
      converts them
  #4  the report describes a render that the backend then refuses outright

Testing those one at a time is what the per-module suites already do, and it
missed all four. They are a *class*, so this file tests the class: whatever a
view says in prose, the ops have to bear out. It runs over every view rather
than the ones known to have failed, because the next divergence will be in a
view nobody is currently suspicious of.

Three checks, deliberately narrow
---------------------------------
Only claims that can be read off the report unambiguously are checked. A loose
parser that guesses at prose would produce failures nobody trusts, and a check
nobody trusts gets deleted.

1. **A named object exists.** Every ``foo_surf_04``-shaped identifier the report
   names must be an object some op actually creates or references.
2. **A level is not labelled with the wrong unit.** Deliberately not "every
   number in the report is an emitted level": a report legitimately states one
   level both ways ("1.5 sigma  =  0.75 absolute"), and the first draft of this
   check fired on that correct text. It now fires only on a real contradiction
   — the same number, labelled as the unit it is not.
3. **A timeline frame shows the surface its number names.** Checked against the
   *backend*, not the prose. The first draft parsed "Skipped: frame 3" as an
   instruction to type ``frame 3`` and flagged a correct sentence; the actual
   #5 divergence is in the lowering, where ``_frames`` numbers by position.

Both of those first drafts are worth recording: a consistency check that reads
prose loosely produces false failures, and a false failure is how a check stops
being believed. Twice the code was right and the check was wrong.

Every view is also run over **degenerate input** — an ensemble with an unusable
frame — because every divergence in the review appeared exactly when something
was skipped. A suite that only sees well-formed data cannot catch this class.

Views that emit no isosurface or timeline have nothing to check for (2) and (3),
and skip explicitly rather than passing silently.
"""

from __future__ import annotations

import re

import pytest

from tests.conftest import iterate_response, make_atoms
from tests.test_mapinfo import write_map
from wiggles_em.backends.pymol import PymolBackend
from wiggles_em.composition import composition_view
from wiggles_em.deformation import deformation_view
from wiggles_em.density import density_view
from wiggles_em.ensembles import ensemble_spread_view, morph_states
from wiggles_em.heterogeneity import forget_ensemble, load_ensemble
from wiggles_em.latent import latent_traverse_view
from wiggles_em.localres import local_resolution_view
from wiggles_em.maps import forget_map, load_map
from wiggles_em.occupancy import altloc_view, occupancy_view
from wiggles_em.occupancy_states import state_occupancy_view
from wiggles_em.populations import Populations, WeightSource
from wiggles_em.port import FakePort
from wiggles_em.provenance import Provenance
from wiggles_em.qscore import qscore_view
from wiggles_em.scene import Frames, Isosurface, Scene, Unit

ROWS = [
    ("A", "1", "MET", "CA", "", 1.0, 20.0),
    ("A", "1", "MET", "CB", "A", 0.6, 22.0),
    ("A", "2", "SER", "CA", "", 1.0, 25.0),
    ("A", "2", "SER", "CB", "B", 0.4, 28.0),
]
COORDS_A = [(0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (2.0, 0.0, 0.0), (3.0, 0.0, 0.0)]
COORDS_B = [(0.5, 0.0, 0.0), (1.5, 0.0, 0.0), (2.5, 0.0, 0.0), (3.5, 0.0, 0.0)]

#: An identifier of the shape views generate for objects they create.
_NAMED_OBJECT = re.compile(r"\b([A-Za-z_][A-Za-z0-9_]*_(?:surf|morph|mesh)_?\d*)\b")

#: "1.5 sigma", "10.75 absolute" — a number with the unit it is stated in.
_STATED_UNIT = re.compile(r"(-?\d+(?:\.\d+)?)\s+(sigma|absolute)\b", re.IGNORECASE)


@pytest.fixture(autouse=True)
def _clean():
    forget_map()
    forget_ensemble()
    yield
    forget_map()
    forget_ensemble()


def _ensemble(tmp_path, rms=(0.5, 0.6, 0.8), name="ens"):
    root = tmp_path / name
    root.mkdir()
    (root / "z.pkl").write_text("x")
    for i, r in enumerate(rms, start=1):
        write_map(root, f"vol_{i}.mrc", rms=r)
    names = [f"{name}_f{i:02d}" for i in range(1, len(rms) + 1)]
    port = FakePort({"get_names": names, "iterate_to_list": []})
    load_ensemble(port, root, name)
    return port, name


def _two_maps(tmp_path):
    main = write_map(tmp_path, "main.mrc", rms=0.5)
    res = write_map(tmp_path, "res.mrc", rms=0.5)
    port = FakePort({"get_names": ["main", "res"], "iterate_to_list": [], "get": "1"})
    load_map(port, main, "main", provenance=Provenance.MEASURED)
    load_map(port, res, "res", provenance=Provenance.MEASURED)
    return port


def _all_views(tmp_path):
    """(label, report, scene, port) for every view, on a minimal real setup.

    The port comes back with the view because check 3 has to render the scene
    to see what the timeline actually does; a claim about lowering cannot be
    checked against the scene alone.
    """
    atoms = make_atoms(ROWS)
    out = []

    def atom_port():
        return FakePort(iterate_response(ROWS))

    out.append(("occupancy_view", *occupancy_view(atoms, "m"), atom_port()))
    out.append(("altloc_view", *altloc_view(atoms, "m"), atom_port()))
    out.append(
        (
            "ensemble_spread_view",
            *ensemble_spread_view(atoms, [COORDS_A, COORDS_B], "m", superposed=True),
            atom_port(),
        )
    )
    out.append(
        (
            "deformation_view",
            *deformation_view(atoms, COORDS_A, COORDS_B, "m", 2, start_state=1, end_state=2),
            atom_port(),
        )
    )
    out.append(("morph_states", *morph_states([4, 4], "m"), atom_port()))
    out.append(
        (
            "composition_view",
            *composition_view("m", {"a": 0.7, "b": 0.3}, lambda _s: 10),
            atom_port(),
        )
    )

    q = tmp_path / "q.xml"
    q.write_text(
        '<?xml version="1.0"?><wwPDB-validation-information>'
        '<ModelledSubgroup chain="A" resnum="1" said="A" Q_score="0.8"/>'
        '<ModelledSubgroup chain="A" resnum="2" said="A" Q_score="0.4"/>'
        "</wwPDB-validation-information>"
    )
    out.append(("qscore_view", *qscore_view(atoms, "m", str(q)), atom_port()))

    density_dir = tmp_path / "d"
    density_dir.mkdir()
    forget_map()
    map_path = write_map(density_dir, "m.mrc", rms=0.5)
    port = FakePort({"get_names": ["m"], "iterate_to_list": []})
    load_map(port, map_path, "m", provenance=Provenance.MEASURED)
    out.append(("density_view", *density_view("m", "polymer"), port))
    forget_map()

    localres_dir = tmp_path / "lr"
    localres_dir.mkdir()
    lr_port = _two_maps(localres_dir)
    out.append(
        ("local_resolution_view", *local_resolution_view("main", "res", normalised=True), lr_port)
    )
    forget_map()

    ensemble_dir = tmp_path / "e"
    ensemble_dir.mkdir()
    ens_port, _ = _ensemble(ensemble_dir)
    out.append(("latent_traverse_view", *latent_traverse_view("ens"), ens_port))
    forget_ensemble()

    # The same view over an ensemble with an UNUSABLE frame. Without this the
    # checks only ever see well-formed input, and every divergence in the
    # review appeared exactly when something was skipped — a frame dropped, a
    # header without statistics. A consistency suite that runs only on the
    # happy path cannot catch the class it exists for.
    gapped_dir = tmp_path / "g"
    gapped_dir.mkdir()
    gap_port, _ = _ensemble(gapped_dir, rms=(0.5, 0.6, 0.0, 0.7, 0.9), name="gap")
    out.append(("latent_traverse_view (gapped)", *latent_traverse_view("gap"), gap_port))
    forget_ensemble()

    # State occupancy, in both of its rendering modes. The two differ in what
    # they are allowed to draw, so checking only one would leave the other's
    # report unguarded — and the non-quantitative branch is the one that has to
    # say most.
    occ_dir = tmp_path / "o"
    occ_dir.mkdir()
    occ_port, _ = _ensemble(occ_dir, name="occ")
    good = Populations.declare([0.6, 0.3, 0.1], WeightSource.DECONVOLVED, temperature_k=298.15)
    out.append(("state_occupancy_view", *state_occupancy_view("occ", good), occ_port))
    counted = Populations.declare([0.6, 0.3, 0.1], WeightSource.LATENT_HISTOGRAM)
    out.append(("state_occupancy_view (counted)", *state_occupancy_view("occ", counted), occ_port))
    forget_ensemble()

    return out


VIEW_LABELS = [
    "occupancy_view",
    "altloc_view",
    "ensemble_spread_view",
    "deformation_view",
    "morph_states",
    "composition_view",
    "qscore_view",
    "density_view",
    "local_resolution_view",
    "latent_traverse_view",
    "latent_traverse_view (gapped)",
    "state_occupancy_view",
    "state_occupancy_view (counted)",
]


@pytest.fixture
def views(tmp_path):
    return {label: (report, scene, port) for label, report, scene, port in _all_views(tmp_path)}


def _op_names(scene) -> set[str]:
    """Every object name any op creates or refers to."""
    names: set[str] = set()
    for op in scene:
        for attr in ("name", "volume"):
            value = getattr(op, attr, None)
            if isinstance(value, str):
                names.add(value)
        if isinstance(op, Frames):
            names.update(op.names)
        sel = getattr(op, "sel", None)
        if sel is not None:
            names.update(str(s.value) for s in sel.walk() if s.kind == "obj")
    return names


def _missing_objects(report: str, scene) -> set[str]:
    """Object names the report claims that no op creates or references."""
    emitted = _op_names(scene)
    return {n for n in _NAMED_OBJECT.findall(report) if n not in emitted}


def _mislabelled_levels(report: str, scene) -> list[tuple[str, str, list[str]]]:
    """(number, unit the report calls it, unit the op actually carries)."""
    emitted = {(round(op.level, 6), op.unit) for op in scene if isinstance(op, Isosurface)}
    out = []
    for number, unit_word in _STATED_UNIT.findall(report):
        value = round(float(number), 6)
        said = Unit.SIGMA if unit_word.lower() == "sigma" else Unit.ABSOLUTE
        carried = {unit for level, unit in emitted if level == value}
        if carried and said not in carried:
            out.append((number, unit_word, sorted(u.value for u in carried)))
    return out


# ── check 1: a named object exists ──────────────────────────────────────────


@pytest.mark.parametrize("label", VIEW_LABELS)
def test_every_object_the_report_names_is_one_the_scene_creates(label, views):
    """A report naming `ens_surf_04` when no such surface was emitted tells the
    user to look at something that is not there — #5's shape exactly."""
    report, scene, _port = views[label]
    emitted = _op_names(scene)
    missing = _missing_objects(report, scene)

    assert not missing, (
        f"{label}'s report names {sorted(missing)}, which no op creates or "
        f"references. Emitted: {sorted(emitted)}"
    )


# ── check 2: a stated unit matches the emitted Unit ─────────────────────────


@pytest.mark.parametrize("label", VIEW_LABELS)
def test_a_level_the_report_names_is_not_labelled_with_the_wrong_unit(label, views):
    """A number the ops actually carry must not be labelled with the other unit.

    Deliberately *not* "every number in the report is an emitted level". A
    report legitimately states one level both ways — density_view writes
    "Contour: 1.5 sigma  =  0.75 absolute" while the op carries 1.5 sigma — and
    a check that demanded both appear as ops would fire on correct code. A
    check nobody trusts gets switched off, so this one only fires on a genuine
    contradiction: the same number, labelled as the unit it is not.

    That is the EMD-30913 error the Unit enum exists to make structural. Its
    published 0.05 is an absolute value and 3.16 sigma; a report calling the
    emitted 0.05 "sigma" sends the reader to contour noise.
    """
    report, scene, _port = views[label]
    if not any(isinstance(op, Isosurface) for op in scene):
        pytest.skip(f"{label} emits no isosurface, so it states no contour")
    if not _STATED_UNIT.findall(report):
        pytest.skip(f"{label} states no contour in its report")

    mislabelled = _mislabelled_levels(report, scene)

    assert not mislabelled, (
        f"{label}'s report labels a level with the unit it is not: "
        f"{mislabelled} (number, called, actually). One of the two is wrong, "
        f"and a level in the wrong unit contours somewhere else entirely."
    )


# ── check 3: a named frame steps to the surface it claims ───────────────────


@pytest.mark.parametrize("label", VIEW_LABELS)
def test_a_timeline_frame_shows_the_surface_its_number_names(label, views):
    """#5 directly, and asserted against the backend rather than the prose.

    The report's claim is explicit: "a surface's number is always the frame it
    was made from", plus "`frame N` to step". The scene honours it — the
    surfaces really do keep their original numbers. The divergence is in the
    *lowering*: `_frames` numbers the movie by position in the surviving list,
    so with a frame skipped, `frame 3` enables `..._surf_04`.

    Parsing the prose cannot see that, and a first attempt at this check tried
    to — it read "Skipped: frame 3" as an instruction to type `frame 3` and
    flagged correct text. So the claim is checked structurally: every timeline
    index must equal the number of the surface it enables.
    """
    report, scene, port = views[label]
    if not any(isinstance(op, Frames) for op in scene):
        pytest.skip(f"{label} emits no timeline")

    PymolBackend(port, normalised=True).render(scene)
    mdo = [(args[0], str(args[1])) for args, _ in port.calls("mdo")]
    if not mdo:
        pytest.skip(f"{label} emitted no movie frames")

    mismatched = []
    for index, command in mdo:
        enabled = re.search(r"enable\s+\S*?_(\d+)\b", command)
        if enabled and int(enabled.group(1)) != int(index):
            mismatched.append((index, enabled.group(1)))

    assert not mismatched, (
        f"{label}: timeline frame(s) {[i for i, _ in mismatched]} enable "
        f"surfaces {[s for _, s in mismatched]}. The report promises a "
        f"surface's number is always the frame it was made from and tells the "
        f"user `frame N` to step, so typing that number shows a different "
        f"frame's density. Report says: "
        f"{'renumbered' if 'renumber' in report else 'gaps kept'}."
    )


# ── check 4: nothing in a report is an unrendered placeholder ───────────────


@pytest.mark.parametrize("label", VIEW_LABELS)
def test_no_report_contains_an_unsubstituted_placeholder(label, views):
    """A report is prose, and prose cannot fail a test — so a broken format
    string reaches the user with the whole suite green.

    This is not hypothetical: an edit to the rigid-motion warning left a literal
    `{n}` in `ensemble_spread_view`'s output and 439 tests passed. Cheap to
    check, and it covers every view at once.
    """
    report, _scene, _port = views[label]

    leftovers = re.findall(r"\{[a-z_][a-z0-9_]*(?::[^}]*)?\}", report, re.IGNORECASE)

    assert not leftovers, (
        f"{label}'s report contains unrendered placeholder(s) {leftovers} — a "
        f"format string that lost its f-prefix, or a .format() that was never "
        f"applied."
    )


# ── the guard on this file ──────────────────────────────────────────────────


def test_the_checks_run_against_every_view():
    """A view added without a row here would be silently unchecked, which is
    the same silence the four findings lived in."""
    import inspect

    import wiggles_em

    exported = {
        name for name in wiggles_em.__all__ if name.endswith("_view") or name == "morph_states"
    }
    # A label may carry a parenthesised scenario suffix; the view it covers is
    # the part before it.
    covered = {label.split(" (")[0] for label in VIEW_LABELS}
    missing = exported - covered

    assert not missing, f"views exported but not covered by the consistency checks: {missing}"
    assert len(VIEW_LABELS) == len(set(VIEW_LABELS)), "duplicate label"
    assert inspect  # the import is the check that this module is importable


# ── the checks, checked ─────────────────────────────────────────────────────


class TestTheChecksHaveTeeth:
    """A consistency check that cannot fail is decoration.

    Check 3 proved itself by catching #5 on the gapped ensemble. Checks 1 and 2
    pass on every view today, which is either good news or a broken checker, and
    those look identical from the outside. So both are run against a report
    deliberately inconsistent with its scene.
    """

    def test_a_named_object_that_does_not_exist_is_caught(self):
        scene = Scene([Isosurface("m_mesh", "m", 1.5, Unit.SIGMA)])

        assert not _missing_objects("Mesh `m_mesh` was drawn.", scene)
        assert _missing_objects("Mesh `ghost_surf_04` was drawn.", scene) == {"ghost_surf_04"}

    def test_a_level_labelled_with_the_wrong_unit_is_caught(self):
        scene = Scene([Isosurface("m_mesh", "m", 0.05, Unit.ABSOLUTE)])

        # EMD-30913's published 0.05 is absolute; calling it sigma sends the
        # reader to contour noise.
        assert _mislabelled_levels("Contour: 0.05 absolute", scene) == []
        assert _mislabelled_levels("Contour: 0.05 sigma", scene) == [
            ("0.05", "sigma", ["absolute"])
        ]

    def test_stating_one_level_both_ways_is_not_a_contradiction(self):
        """The false positive the first draft produced, kept as a test so it
        cannot come back: density_view writes both forms of one level."""
        scene = Scene([Isosurface("m_mesh", "m", 1.5, Unit.SIGMA)])

        assert _mislabelled_levels("Contour: 1.5 sigma  =  0.75 absolute", scene) == []
