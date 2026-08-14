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
from typing import Any

import pytest

from tests.conftest import make_atoms
from wiggles_em.occupancy import altloc_view, occupancy_view
from wiggles_em.scene import (
    Arrow,
    ColorFlat,
    ColorSurfaceByMap,
    Scene,
    Sel,
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
        ("grey70", (0.7, 0.7, 0.7)),
        ("grey50", (0.5, 0.5, 0.5)),
        ("skyblue", (0.34, 0.63, 0.83)),
        ("salmon", (1.0, 0.6, 0.6)),
        ("red", (1.0, 0.0, 0.0)),
    ],
)
def test_names_resolve_to_the_value_they_always_had(name, rgb):
    """Pinned, not merely "is a triple". The conversion is only safe if the
    colour drawn is unchanged, and PyMOL renders the RGB now rather than
    looking the name up itself."""
    assert resolve_colour(name) == rgb


def test_rgb_passes_through_unchanged():
    assert resolve_colour((0.1, 0.2, 0.3)) == (0.1, 0.2, 0.3)


def test_the_grey_ramp_is_computed_not_transcribed():
    """PyMOL's greyNN is a uniform NN/100 ramp, so computing it cannot drift
    from PyMOL the way a transcribed table entry can."""
    assert resolve_colour("grey00") == (0.0, 0.0, 0.0)
    assert resolve_colour("gray25") == (0.25, 0.25, 0.25)
    assert resolve_colour("grey99") == (0.99, 0.99, 0.99)


def test_an_unknown_name_is_refused_not_approximated():
    """The whole point of the table. An approximated colour draws a plausible
    picture and returns cleanly, which is the expensive failure here."""
    with pytest.raises(ValueError) as excinfo:
        resolve_colour("chartreuse")

    message = str(excinfo.value)
    assert "chartreuse" in message
    # The remedy has to be followable, or the refusal just blocks the caller.
    assert "RGB triple" in message


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
