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
from wiggles_em.scene import Frames, Isosurface, ScalarField, Scene, Sel, Unit


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


def test_the_timeline_lands_on_a_frame_that_shows_something():
    """`frame 1` was hard-coded, but frame 1 is a gap whenever frame 1's header
    has no usable rms — and a gap frame disables everything, so a successful
    traversal ended up rendering blank under a report listing its surfaces.

    The old positional numbering always enabled the first surface; carrying real
    frame numbers broke that without anything noticing.
    """
    port = FakePort()
    # Frame 1 skipped: surfaces are 2 and 3.
    frames = Frames(("s_02", "s_03"), (2, 3), build_timeline=True)
    PymolBackend(port, normalised=None).render(Scene([frames]))

    (landed,), _ = port.calls("frame")[0]
    commands = {args[0]: str(args[1]) for args, _ in port.calls("mdo")}

    assert "enable" in commands[landed], (
        f"the timeline lands on frame {landed}, whose command is "
        f"{commands[landed]!r} — the session shows nothing at all"
    )
    assert landed == 2, f"should land on the first real frame, not {landed}"


def test_the_timeline_still_lands_on_frame_one_when_frame_one_exists():
    port = FakePort()
    PymolBackend(port, normalised=None).render(
        Scene([Frames(("s_01", "s_02"), (1, 2), build_timeline=True)])
    )

    assert port.calls("frame")[0][0] == (1,), port.call_log


def test_duplicate_frame_numbers_are_refused():
    """`_frames` builds dict(zip(numbers, names)), so a repeated number silently
    keeps one surface and drops the other — it is then never enabled on any
    frame, and nothing looks wrong. `ScalarField` already refuses duplicate
    keys for exactly this reason."""
    with pytest.raises(ValueError, match="duplicate"):
        Frames(("s_01", "s_02"), (1, 1), build_timeline=True)


def test_frame_numbers_below_one_are_refused():
    """The timeline runs 1..max, so a number of 0 or less is dropped without
    trace — a surface silently absent from a movie that claims to hold it."""
    with pytest.raises(ValueError, match="1 or greater"):
        Frames(("s_00",), (0,), build_timeline=True)


def test_a_carried_equivalent_is_recorded_like_any_other_conversion():
    """`converted` is documented as what the report layer reads back so it can
    state both units without recomputing. The `equivalent` short-circuit
    returned early and recorded nothing — for precisely the volumes (ensemble
    frames in an unnormalised session) the field exists to serve."""
    port = FakePort()
    backend = PymolBackend(port, normalised=False)
    backend.render(
        Scene([
            Isosurface("s", "ens_f01", level=1.5, unit=Unit.SIGMA, equivalent=0.75),
        ])
    )

    # Keyed by surface name, not volume — see TestConvertedIsKeyedTheWayItIsDocumented.
    assert backend.converted["s"] == (1.5, 0.75), backend.converted


class TestScalarFieldRefusesAKeyItCannotMatch:
    """The per-atom key is a pair — `(model, rank)` for atoms, `(chain, resi)`
    for residues. The old, now-wrong contract was a 4-tuple
    `(chain, resi, name, alt)`, and `per_atom`'s docstring spells out what
    happens to a field built that way: it matches nothing, nothing raises, and
    the structure is coloured by whatever the B-factor column already held under
    a legend naming a quantity that was never drawn.

    Documenting a silent-wrongness path is weaker than closing it, and the rest
    of this dataclass closes its hazards — length mismatch and duplicate keys
    both raise.
    """

    def test_the_old_four_part_key_is_refused(self):
        with pytest.raises(ValueError, match="two parts"):
            ScalarField.per_atom([(("A", "12", "CA", ""), 0.5)])

    def test_a_pair_is_accepted(self):
        field = ScalarField.per_atom([(("m", "1"), 0.5), (("m", "2"), 0.25)])
        assert len(field) == 2

    def test_per_residue_keys_are_checked_the_same_way(self):
        with pytest.raises(ValueError, match="two parts"):
            ScalarField.per_residue([(("A", "12", "CA"), 0.5)])

    def test_an_empty_field_is_still_allowed(self):
        assert len(ScalarField.per_atom([])) == 0


def test_the_timeline_cost_is_bounded_by_the_highest_frame_number():
    """Pinned because it is a real cost and an easy one to make worse.

    The timeline spans frame *numbers*, not surfaces, so an ensemble whose only
    usable frames are 491-500 builds 500 movie frames — 490 of them clearing the
    view. The round trips cannot be avoided without breaking what happens when a
    user jumps straight into a gap, but the identical gap command should be
    built once.
    """
    names = tuple(f"s_{i:03d}" for i in range(491, 501))
    port = FakePort()
    PymolBackend(port, normalised=None).render(
        Scene([Frames(names, tuple(range(491, 501)), build_timeline=True)])
    )

    commands = [str(args[1]) for args, _ in port.calls("mdo")]
    gaps = [c for c in commands if "enable" not in c]

    assert len(commands) == 500, "one command per frame number, by design"
    assert len(gaps) == 490
    assert len(set(gaps)) == 1, "every gap frame clears the same set of surfaces"


class TestTheArityGuardSurvivesWhatItCatches:
    """A guard that raises TypeError on the malformed input it exists to reject
    has replaced a silent wrong answer with a worse error message.

    Both shapes below are plausible: a bare `rank` is an easy slip given how
    much the docs emphasise it, and the old 4-tuple contract's `resi` was an int
    in some callers and a str in others.
    """

    def test_a_non_tuple_key_gets_the_guards_own_message(self):
        with pytest.raises(ValueError, match="two parts"):
            ScalarField.per_atom([(3, 0.5)])

    def test_keys_that_cannot_be_ordered_get_the_guards_own_message(self):
        with pytest.raises(ValueError, match="two parts"):
            ScalarField.per_atom([(("A", 12, "CA"), 0.5), (("A", "12", "CA"), 0.6)])

    def test_a_string_key_is_refused_rather_than_split_into_characters(self):
        """A str has a len() and is iterable, so it would sail past a naive
        length check and be reported as a two-part key when it is one value."""
        with pytest.raises(ValueError, match="two parts"):
            ScalarField.per_atom([("m1", 0.5)])


class TestConvertedIsKeyedTheWayItIsDocumented:
    """`converted` is documented as "by surface name" and was written by volume.

    Nothing in src/ reads it, so neither the key nor the claim was load-bearing
    yet — which is exactly when a contract is cheapest to fix and easiest to
    leave wrong. A host implementing the report layer against the docstring
    would look up a key that is never written.
    """

    def test_the_key_is_the_surface_not_the_volume(self):
        port = FakePort()
        backend = PymolBackend(port, normalised=False)
        backend.render(
            Scene([Isosurface("surf_a", "vol", level=1.5, unit=Unit.SIGMA, equivalent=0.75)])
        )

        assert "surf_a" in backend.converted, backend.converted
        assert "vol" not in backend.converted, backend.converted

    def test_two_contours_of_one_volume_do_not_collide(self):
        """The reason the surface is the right key: a volume can carry several
        surfaces at different levels, and under a volume key the second silently
        overwrites the first."""
        port = FakePort()
        backend = PymolBackend(port, normalised=False)
        backend.render(
            Scene([
                Isosurface("low", "vol", level=1.0, unit=Unit.SIGMA, equivalent=0.5),
                Isosurface("high", "vol", level=3.0, unit=Unit.SIGMA, equivalent=1.5),
            ])
        )

        assert backend.converted["low"] == (1.0, 0.5)
        assert backend.converted["high"] == (3.0, 1.5)
