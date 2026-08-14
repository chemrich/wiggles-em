"""A Scene may not name one viewer's colour palette.

``Colour = str | tuple[float, float, float]``, and every string a view used to
emit was a *PyMOL* colour name — ``grey70``, ``skyblue``. A second viewer can
only honour those by reimplementing PyMOL's table, which is exactly what
protean had to do (``backends/molstar.py::_COLOUR_NAMES``) and what
``composition._blend`` had already done here. Three copies of one table, and
the two that existed before this module had **disjoint** name sets.

So the invariant is on the *Scene*, not on any one view: a colour that reaches
a backend is RGB, and a name is resolved before it gets there. Names remain
perfectly good **arguments** — ``composition_view(palette=("red", "skyblue"))``
still works — because a view's signature is not the seam. The Scene is.

The sweep below is annotation-driven rather than a list of ops, so an op added
later with a ``Colour`` field is covered without anyone remembering to come
here. :func:`assert_no_palette_names` is also called from
``test_defaults._assert_renders``, which means **every view is checked on its
own defaults** — the path where a hard-coded name is most likely to hide.
"""

from __future__ import annotations

import dataclasses
import shutil
import subprocess
import tempfile
from pathlib import Path
from typing import Any

import pytest

from tests.conftest import make_atoms
from wiggles_em.occupancy import altloc_view, occupancy_view
from wiggles_em.scene import (
    BLUE_WHITE_RED,
    RED_WHITE_BLUE,
    RED_YELLOW_GREEN,
    Arrow,
    ColorByScalar,
    ColorFlat,
    ColorSurfaceByMap,
    ScalarField,
    Scene,
    Sel,
    ramp,
    resolve_colour,
)

ROWS = [
    ("A", "1", "MET", "CA", "", 1.0, 20.0),
    ("A", "1", "MET", "CB", "A", 0.6, 22.0),
    ("A", "2", "SER", "CA", "", 1.0, 25.0),
    ("A", "2", "SER", "CB", "B", 0.4, 28.0),
]


def _is_rgb(value: Any) -> bool:
    """An RGB triple, as opposed to a palette of three colours."""
    return (
        isinstance(value, tuple)
        and len(value) == 3
        and all(isinstance(component, int | float) for component in value)
    )


def _colours_of_field(value: Any) -> list[Any]:
    """Split a ``Colour``-typed value into the individual colours in it."""
    if _is_rgb(value):
        return [value]
    if isinstance(value, list | tuple):
        return list(value)
    return [value]


def colours_in(op: Any) -> list[Any]:
    """Every value sitting in a ``Colour``-typed position on ``op``, nested.

    Driven by ``field.type``, which dataclasses keeps as the annotation's
    *source text* (``'Colour'``, ``'tuple[Colour, ...]'``) — so this finds a
    new carrier by its declaration rather than by a list maintained here.
    Non-Colour fields are still descended into, because :class:`Arrow` carries
    one and sits inside ``Arrows.segments``.
    """
    found: list[Any] = []
    for field in dataclasses.fields(op):
        value = getattr(op, field.name)
        if "Colour" in str(field.type):
            found.extend(_colours_of_field(value))
            continue
        for item in value if isinstance(value, list | tuple) else [value]:
            if dataclasses.is_dataclass(item) and not isinstance(item, type):
                found.extend(colours_in(item))
    return found


def palette_names_in(scene: Scene) -> list[str]:
    """Every colour in ``scene`` that is a name rather than RGB."""
    names: list[str] = []
    for op in scene:
        names += [c for c in colours_in(op) if isinstance(c, str)]
    return names


def assert_no_palette_names(scene: Scene) -> None:
    """The seam invariant. Imported by ``test_defaults`` for its whole sweep."""
    names = palette_names_in(scene)
    assert not names, (
        f"this Scene names {len(names)} colour(s) from one viewer's palette: "
        f"{sorted(set(names))}. A second viewer can only honour a name by "
        f"reimplementing PyMOL's table. Resolve it in the view."
    )


# ── the walker itself, checked against hand-built scenes ────────────────────


def test_walker_finds_a_flat_colour_name():
    """Necessity: the sweep is worthless if it cannot see a name."""
    assert palette_names_in(Scene([ColorFlat(Sel.obj("o"), "grey70")])) == ["grey70"]


def test_walker_passes_an_rgb_triple():
    """An RGB triple is three floats, not three colours."""
    assert palette_names_in(Scene([ColorFlat(Sel.obj("o"), (0.7, 0.7, 0.7))])) == []


def test_walker_finds_names_inside_a_palette():
    """``ColorSurfaceByMap.palette`` is ``tuple[Colour, ...]`` — a sequence."""
    op = ColorSurfaceByMap("s", "v", (1.0, 2.0), ("blue", (0.0, 1.0, 0.0), "red"))
    assert palette_names_in(Scene([op])) == ["blue", "red"]


def test_walker_descends_into_arrows():
    """``Arrow`` carries a Colour and is nested inside ``Arrows.segments``."""
    from wiggles_em.scene import Arrows

    arrow = Arrow((0.0, 0.0, 0.0), (1.0, 0.0, 0.0), "skyblue")
    assert palette_names_in(Scene([Arrows((arrow,), "o")])) == ["skyblue"]


# ── the views ───────────────────────────────────────────────────────────────


def test_altloc_view_emits_no_palette_names():
    _report, scene = altloc_view(make_atoms(ROWS), "obj")
    assert_no_palette_names(scene)


def test_occupancy_view_emits_no_palette_names():
    _report, scene = occupancy_view(make_atoms(ROWS), "obj")
    assert_no_palette_names(scene)


# ── the resolver ────────────────────────────────────────────────────────────


@pytest.mark.parametrize(
    ("name", "rgb"),
    [
        ("grey70", (70 / 99, 70 / 99, 70 / 99)),
        ("grey50", (50 / 99, 50 / 99, 50 / 99)),
        ("skyblue", (0.2, 0.5, 0.8)),
        ("salmon", (1.0, 0.6, 0.6)),
        ("red", (1.0, 0.0, 0.0)),
    ],
)
def test_names_resolve_to_pymols_value(name, rgb):
    """Pinned, not merely "is a triple". The conversion is only safe if the
    colour drawn is unchanged, and PyMOL renders the RGB now rather than
    looking the name up itself.

    These were wrong once. ``skyblue`` was ``(0.34, 0.63, 0.83)`` here, taken
    from two existing tables that agreed with each other and not with PyMOL.
    Agreement between transcriptions is not evidence — see
    :func:`test_the_palette_matches_pymols_own`, which asks the source.
    """
    assert resolve_colour(name) == pytest.approx(rgb)


def test_the_palette_matches_pymols_own():
    """Check every entry against a real PyMOL, when one is installed.

    The package must import without PyMOL, so this shells out to the binary
    and skips when it is absent — it will skip in CI and run on a workstation.
    That is worth having anyway: a transcribed table is exactly the artefact
    that drifts silently, and the two it replaced were both wrong.

    ``-k`` skips the user's pymolrc, so this loads no plugins and cannot touch
    a session they have open.
    """
    pymol = shutil.which("pymol")
    if pymol is None:
        pytest.skip("no pymol on PATH; the table cannot be checked against it")

    from wiggles_em.scene import _PALETTE

    names = [*_PALETTE, "grey00", "grey07", "gray25", "grey50", "grey70", "grey99"]
    script = Path(tempfile.mkdtemp()) / "dump.py"
    script.write_text(
        "from pymol import cmd\n"
        f"for n in {names!r}:\n"
        '    print("VAL", n, *cmd.get_color_tuple(n))\n'
    )
    out = subprocess.run([pymol, "-ckq", str(script)], capture_output=True, text=True, timeout=120)
    theirs = {
        parts[1]: tuple(float(c) for c in parts[2:5])
        for line in out.stdout.splitlines()
        if (parts := line.split()) and parts[0] == "VAL"
    }
    assert theirs, f"pymol produced no colours: {out.stdout[-400:]} {out.stderr[-400:]}"

    wrong = {
        name: (resolve_colour(name), theirs[name])
        for name in names
        if name in theirs and resolve_colour(name) != pytest.approx(theirs[name], abs=1e-6)
    }
    assert not wrong, "\n".join(
        f"  {name}: this package says {ours}, PyMOL says {pymols}"
        for name, (ours, pymols) in sorted(wrong.items())
    )


def test_rgb_passes_through_unchanged():
    assert resolve_colour((0.1, 0.2, 0.3)) == (0.1, 0.2, 0.3)


def test_the_grey_ramp_is_computed_not_transcribed():
    """PyMOL's greyNN is a uniform ramp, so computing it cannot drift the way a
    transcribed entry can — but only if the divisor is right.

    It is **99, not 100**: the ramp is inclusive, so ``grey99`` is white. This
    test asserted 0.99 and passed, which is what a test written from the same
    wrong assumption as the code always does.
    """
    assert resolve_colour("grey00") == (0.0, 0.0, 0.0)
    assert resolve_colour("grey99") == (1.0, 1.0, 1.0)
    assert resolve_colour("gray25") == pytest.approx((25 / 99, 25 / 99, 25 / 99))
    assert resolve_colour("grey70") == pytest.approx((0.70707, 0.70707, 0.70707), abs=1e-5)


def test_an_unknown_name_is_refused_not_approximated():
    """The whole point of the table. An approximated colour draws a plausible
    picture and returns cleanly, which is the expensive failure here."""
    with pytest.raises(ValueError) as excinfo:
        resolve_colour("chartreuse")

    message = str(excinfo.value)
    assert "chartreuse" in message
    # The remedy has to be followable, or the refusal just blocks the caller.
    assert "RGB triple" in message


def test_an_rgb_triple_on_the_wrong_scale_is_refused():
    """0-255 is the other common spelling and PyMOL's `set_color` accepts both,
    so a wrong-scale triple renders correctly in the one backend tested here
    and is broken in every other viewer."""
    with pytest.raises(ValueError, match="255"):
        resolve_colour((255, 0, 0))


def test_an_unresolvable_arrow_colour_fails_at_construction():
    """Not mid-render. CGO resolves inline while drawing, so a bad name would
    otherwise surface as a bare ValueError over a half-applied scene."""
    from wiggles_em.scene import Arrow as _Arrow

    with pytest.raises(ValueError, match="chartreuse"):
        _Arrow((0.0, 0.0, 0.0), (1.0, 0.0, 0.0), "chartreuse")


def test_a_named_arrow_colour_reaches_pymol_as_that_colour():
    """CGO carries RGB inline, so there is no `color` call to hand a name to.
    The backend used to substitute **white** for any name — invisible to the
    views, which all emit RGB, and silent in the one op whose point is
    direction. Hand-built scenes are a supported caller.
    """
    from tests.conftest import render
    from wiggles_em.scene import Arrows

    arrow = Arrow((0.0, 0.0, 0.0), (1.0, 0.0, 0.0), "red")
    d = render(("report", Scene([Arrows((arrow,), "obj")])))

    # CGO layout: opcode, start(3), shaft_end(3), radius, then r,g,b.
    buffer = d.port.calls("load_cgo")[0][0][0]
    assert tuple(buffer[8:11]) == (1.0, 0.0, 0.0), "red arrow did not reach PyMOL as red"
    assert tuple(buffer[8:11]) != (1.0, 1.0, 1.0), "still substituting white"


# ── scalar ramps ────────────────────────────────────────────────────────────


def test_ramp_reads_pymols_underscore_spelling():
    assert ramp("blue_white_red") == ((0.0, 0.0, 1.0), (1.0, 1.0, 1.0), (1.0, 0.0, 0.0))


def test_ramp_needs_something_to_interpolate_between():
    with pytest.raises(ValueError, match="at least two stops"):
        ramp("red")


def test_the_two_default_ramps_are_exact_reverses():
    """The fact that makes the reversal detectable. A viewer with one fixed
    ramp can honour exactly one of these, and until the stops were in the Scene
    it had no way to tell which it had been handed."""
    assert tuple(reversed(BLUE_WHITE_RED)) == RED_WHITE_BLUE


def test_a_palette_name_is_refused_with_the_remedy():
    """A bare string is the trap: it is iterable, so it would silently become
    one stop per *character*."""
    field = ScalarField.per_atom([(("m", "0"), 1.0)])
    with pytest.raises(ValueError, match="ramp"):
        ColorByScalar(Sel.obj("o"), field, (0.0, 1.0), palette="red_white_blue")


def test_ramp_direction_survives_the_trip_to_pymol():
    """End to end: a scene asking for blue→red must reach `spectrum` as a
    palette whose first stop is blue and last is red. Ordering is the whole
    claim — reversing it reverses the reading of every number on screen while
    the render looks entirely normal."""
    from tests.conftest import render

    field = ScalarField.per_atom([(("m", "0"), 1.0)])
    scene = Scene([ColorByScalar(Sel.obj("obj"), field, (0.0, 1.0), palette=BLUE_WHITE_RED)])
    d = render(("report", scene), [("A", "1", "MET", "CA", "", 1.0, 20.0)])

    (_expression, palette, *_rest), _kwargs = d.port.calls("spectrum")[0]
    stops = palette.split("_")
    assert len(stops) == 3, palette
    defined = {a[0]: tuple(a[1]) for a, _ in d.port.calls("set_color")}
    assert [defined[s] for s in stops] == list(BLUE_WHITE_RED), palette


def test_generated_colour_names_carry_no_underscore():
    """`spectrum` splits its palette on underscores, so a generated name with
    one in it would be read as two colours that do not exist. This is the
    constraint that makes the naming scheme what it is."""
    from tests.conftest import render

    field = ScalarField.per_atom([(("m", "0"), 1.0)])
    scene = Scene([ColorByScalar(Sel.obj("obj"), field, (0.0, 1.0), palette=RED_YELLOW_GREEN)])
    d = render(("report", scene), [("A", "1", "MET", "CA", "", 1.0, 20.0)])

    generated = [a[0] for a, _ in d.port.calls("set_color")]
    assert generated, "no colours were defined"
    assert all("_" not in name for name in generated), generated


def test_occupancy_ramp_runs_the_direction_its_report_claims():
    """The report says "red (q=0) → white (q=0.5) → blue (q=1)". Reversing the
    palette left every one of 546 tests passing while the prose described the
    opposite of the picture — found by mutation, which is the only reason this
    test exists.
    """
    report, scene = occupancy_view(make_atoms(ROWS), "obj")
    (op,) = scene.of(ColorByScalar)

    assert op.palette[0] == resolve_colour("red"), "low end is not red"
    assert op.palette[-1] == resolve_colour("blue"), "high end is not blue"
    # Both halves, so a reversal fails the first pair and a silent rewording
    # of the legend fails the second.
    assert "red (q=0)" in report
    assert "blue (q=1)" in report


def test_qscore_ramp_runs_the_direction_its_report_claims(tmp_path):
    """Low Q is bad, so this ramp carries a judgement: red is unresolvable."""
    from wiggles_em.qscore import qscore_view

    path = tmp_path / "q.xml"
    path.write_text(
        '<?xml version="1.0"?><wwPDB-validation-information>'
        '<ModelledSubgroup chain="A" resnum="1" said="A" Q_score="0.8"/>'
        '<ModelledSubgroup chain="A" resnum="2" said="A" Q_score="0.4"/>'
        "</wwPDB-validation-information>"
    )
    report, scene = qscore_view(make_atoms(ROWS), "obj", str(path))
    (op,) = scene.of(ColorByScalar)

    assert op.palette[0] == resolve_colour("red"), "Q=0 is not red"
    assert op.palette[-1] == resolve_colour("green"), "Q=1 is not green"
    assert "red (Q=0" in report
    assert "green (Q=1)" in report


def test_displacement_and_spread_ramp_low_to_high():
    """Both quantities are magnitudes where more should read as hotter, so
    blue is the low end. Reversed, a rigid core would read as the mobile part.
    """
    from wiggles_em.deformation import deformation_view
    from wiggles_em.ensembles import ensemble_spread_view

    start = [(0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (2.0, 0.0, 0.0), (3.0, 0.0, 0.0)]
    end = [(0.5, 0.0, 0.0), (1.5, 0.0, 0.0), (2.5, 0.0, 0.0), (3.5, 0.0, 0.0)]

    _r, deformation = deformation_view(
        make_atoms(ROWS), start, end, "obj", 2, start_state=1, end_state=2
    )
    _r2, spread = ensemble_spread_view(make_atoms(ROWS), [start, end], "obj", superposed=True)

    for scene in (deformation, spread):
        (op,) = scene.of(ColorByScalar)
        assert op.palette == BLUE_WHITE_RED
        assert op.palette[0] == resolve_colour("blue")
        assert op.palette[-1] == resolve_colour("red")


def test_composition_refuses_a_palette_name_it_cannot_resolve():
    """`_blend` used to fall back to red/skyblue for a name it did not know, so
    this drew the *default* ramp and reported success. Caught by mutation
    testing: fixing it without this test left the fix unguarded."""
    from wiggles_em.composition import composition_view

    with pytest.raises(ValueError, match="purple"):
        composition_view(
            "obj",
            {"state_a": 0.7, "state_b": 0.3},
            lambda _sel: 10,
            palette=("purple", "teal"),
        )
