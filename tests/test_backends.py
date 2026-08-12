"""Tests for the PyMOL lowering itself, rather than for any one view.

A backend is where the scene stops being a value and starts being commands, so
these assert on the commands. Everything here is about the lowering being both
*exact* — no wildcard that could touch something the view did not make — and
*proportionate*, since a selection PyMOL has to evaluate term by term is a
selection that times out on the structures this package is built for.
"""

from __future__ import annotations

from wiggles_em.backends.pymol import PymolBackend, render_selection
from wiggles_em.port import FakePort
from wiggles_em.scene import Frames, Scene, Sel


def test_a_large_residue_set_does_not_become_one_term_per_residue():
    """Sel.residues exists so a backend can emit one compact expression.

    qscore_view on 9C0K carries 1123 scored residues. Lowered as 1123
    parenthesised `or` terms the string is ~30 KB, PyMOL evaluates every term
    across the whole object, and the call runs past the 10 s port timeout — on
    a structure the tool is explicitly designed for.
    """
    residues = [("A", str(i)) for i in range(1, 1124)]
    text = render_selection(Sel.obj("m") & Sel.residues(residues))

    assert text.count(" or ") < 20, f"{text.count(' or ')} or-terms, {len(text)} chars"
    assert len(text) < 8000, f"{len(text)} chars"


def test_residues_across_chains_stay_grouped_by_chain():
    residues = [("A", "1"), ("A", "2"), ("B", "7")]
    text = render_selection(Sel.residues(residues))

    assert text.count(" or ") == 1, text
    assert "chain A" in text and "chain B" in text
    # And it must still mean the same set.
    assert "1" in text and "2" in text and "7" in text


def test_a_residue_set_is_not_reordered_into_a_different_set():
    """Compactness must not merge 5 and 7 into 5-7."""
    text = render_selection(Sel.residues([("A", "5"), ("A", "7")]))
    assert "6" not in text, text


def test_the_movie_is_entered_not_merely_built():
    """mdo attaches a command to a frame; it runs when the frame is entered.

    Without a final `frame` call the timeline exists and nothing has executed
    it, so every isosurface stays enabled and the user sees all of them
    superimposed instead of a stepped trajectory.
    """
    port = FakePort()
    PymolBackend(port).render(Scene([Frames(("s_01", "s_02", "s_03"), build_timeline=True)]))

    assert port.queried("frame"), port.call_log
    order = [name for name, _, _ in port.queries]
    assert order.index("mset") < order.index("frame")
    assert order.index("mdo") < order.index("frame")


def test_the_movie_never_disables_anything_it_did_not_make():
    """The repo's own rule: cleanup is scoped to what the view created.

    `disable v_*` also switches off an unrelated v_model or v_map every time
    the user steps a frame.
    """
    port = FakePort()
    names = ("v_01", "v_02")
    PymolBackend(port).render(Scene([Frames(names, build_timeline=True)]))

    for args, _ in port.calls("mdo"):
        command = args[1]
        assert "*" not in command, command
        for word in command.replace(";", " ").split():
            if word not in ("disable", "enable"):
                assert word in names, f"{word!r} is not one of {names}"
