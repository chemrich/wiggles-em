"""Tests for the PyMOL lowering itself, rather than for any one view.

A backend is where the scene stops being a value and starts being commands, so
these assert on the commands. Everything here is about the lowering being both
*exact* — no wildcard that could touch something the view did not make — and
*proportionate*, since a selection PyMOL has to evaluate term by term is a
selection that times out on the structures this package is built for.
"""

from __future__ import annotations

import pytest

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
    assert 'chain "A"' in text and 'chain "B"' in text
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
    frames = Frames(("s_01", "s_02", "s_03"), (1, 2, 3), build_timeline=True)
    # normalised=None: this test says nothing about volumes, and the
    # argument is required so it cannot be assumed on a test's behalf.
    PymolBackend(port, normalised=None).render(Scene([frames]))

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
    PymolBackend(port, normalised=None).render(Scene([Frames(names, (1, 2), build_timeline=True)]))

    for args, _ in port.calls("mdo"):
        command = args[1]
        assert "*" not in command, command
        for word in command.replace(";", " ").split():
            if word not in ("disable", "enable"):
                assert word in names, f"{word!r} is not one of {names}"


# ── selection quoting ───────────────────────────────────────────────────────
#
# Ported from MCPymol PR #57, which found these against PyMOL 3.1.0 rather
# than reasoning about them. Both shapes occur constantly in real depositions,
# and both make PyMOL read the selection as something other than what was
# meant — silently, and in the direction that destroys data.


def test_a_blank_chain_cannot_swallow_the_next_token():
    """A file with no chain ID leaves `chain  and resi 2`, where PyMOL takes
    `and` as the chain name and the selection stops being scoped to the object.

    Checked upstream on PyMOL 3.1.0: in a session holding a 7-atom object and
    a 10-atom one, that selection matched 10. qscore_view defaults a missing
    chain to "", so it ran `alter <every atom loaded>, b=<one Q-score>` — and
    restore_bfactors only restores the object it is given, so every other
    structure's B-factors were gone with no way back.
    """
    text = render_selection(Sel.obj("gly") & Sel.residues([("", "2")]))

    assert 'chain ""' in text, text
    assert "chain  and" not in text, text


def test_a_negative_residue_number_is_not_a_range():
    """`resi -3` means 1-3 to PyMOL, so one residue's value lands on three.

    Negative auth numbering marks expression-tag remnants and is routine in
    NMR and EM entries.

    Quoting alone does NOT fix this, though both this package and MCPymol
    claimed it did: `resi "-3"` is still read as the range. Checked against
    PyMOL 3.1.0 on an object holding residues -3, 1 and 2, where it matched all
    six atoms instead of two. The backslash is the part that does the work —
    see `test_selection_live.py`, which asserts the atom counts rather than the
    spelling.
    """
    text = render_selection(Sel.residues([("A", "-3")]))

    assert 'resi "\\-3"' in text, text
    assert 'resi "-3"' not in text, text
    assert "resi -3" not in text.replace('resi "\\-3"', ""), text


def test_an_insertion_code_survives_quoting():
    text = render_selection(Sel.residues([("A", "52A")]))
    assert 'resi "52A"' in text, text


def test_a_quote_inside_an_identifier_is_refused_not_guessed_at():
    """Quoting is only safe while the value cannot close it. Nothing in the
    PDB or mmCIF grammar puts a double quote in one of these, so the file is
    corrupt — and proceeding anyway is how the blank-chain bug behaved."""
    with pytest.raises(ValueError, match="double quote"):
        render_selection(Sel.residues([('A" or all and chain "B', "1")]))

    with pytest.raises(ValueError, match="double quote"):
        render_selection(Sel.obj('m" or all'))


def test_ordinary_residues_still_group_into_one_compact_term():
    """The grouping that keeps 1123 residues out of a 1123-term disjunction
    must survive the quoting, or the fix for one finding undoes the other."""
    text = render_selection(Sel.residues([("A", str(i)) for i in range(1, 1124)]))

    assert text.count(" or ") < 20, f"{text.count(' or ')} or-terms"
    assert '"1"+"2"+"3"' in text, text


def test_a_quoted_value_composes_into_the_plus_list():
    """The old code split plainly-numeric residues onto a `+` list and sent
    everything else the long way, because a quoted value in a `+` list was
    "grammar this package has not checked".

    It has now been checked against PyMOL 3.1.0, and the grammar accepts it:
    `resi "\\-3"+"1"+"2"` matches three residues. So there is one path rather
    than two that could drift apart — which is how half these findings started.
    """
    text = render_selection(Sel.residues([("A", "1"), ("A", "-3"), ("A", "2")]))

    # One compact term, negative value included rather than exiled to an `or`.
    assert text == '(chain "A" and resi "1"+"\\-3"+"2")', text
    assert " or " not in text, text
    # The chain is named once however the residues are written.
    assert text.count("chain") == 1, text


class TestFramesEdgeCases:
    """Guards that existed in the code the timeline rewrite replaced.

    `Frames` is public scene API, so it is reachable with shapes that
    `latent_traverse_view` never produces — it guards on `len(surfaces) > 1`
    and never emits an empty one. Rewriting a block rather than editing it is
    how both of these were lost: the invariants were in the replaced code and
    nothing named them.
    """

    def test_a_single_surface_timeline_emits_a_clean_command(self):
        """With one surface there are no siblings to disable, and the old code
        said so explicitly. The rewrite left a leading empty statement."""
        port = FakePort()
        PymolBackend(port, normalised=None).render(
            Scene([Frames(("s_01",), (1,), build_timeline=True)])
        )

        (args, _kwargs) = port.calls("mdo")[0]
        command = str(args[1])
        assert not command.strip().startswith(";"), f"leading empty statement: {command!r}"
        assert command.strip() == "enable s_01", command

    def test_an_empty_timeline_draws_nothing_rather_than_raising(self):
        """`max()` over no frames raised ValueError — an internal error escaping
        as the viewer's answer."""
        port = FakePort()
        PymolBackend(port, normalised=None).render(
            Scene([Frames((), (), build_timeline=True)])
        )

        assert not port.calls("mdo"), port.call_log
        assert not port.calls("mset"), port.call_log

    def test_names_and_numbers_must_still_correspond(self):
        """The length check is the one guard the rewrite did keep."""
        with pytest.raises(ValueError, match="frame names but"):
            Frames(("s_01", "s_02"), (1,), build_timeline=True)
