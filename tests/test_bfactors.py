"""Tests for B-factor save/restore around destructive views."""

from __future__ import annotations

import json

import pytest

from wiggles_em.atoms import Atom
from wiggles_em.bfactors import (
    bfactors_destroyed,
    clear_stash,
    destroyed_note,
    has_stash,
    mark_bfactors_destroyed,
    preservation_note,
    restore_bfactors,
    stash_bfactors,
)
from wiggles_em.port import FakePort, PortError

ATOMS = [
    Atom("A", "1", "MET", "CA", "", 1.0, 20.5, "m", 1),
    Atom("A", "1", "MET", "CB", "", 1.0, 22.0, "m", 2),
    Atom("A", "2", "SER", "CA", "A", 0.6, 31.25, "m", 3),
]


def test_stash_records_every_atom():
    assert stash_bfactors("obj", ATOMS) == 3
    assert has_stash("obj")


def test_stash_is_per_object():
    stash_bfactors("a", ATOMS)
    assert has_stash("a")
    assert not has_stash("b")


def test_clear_one_and_clear_all():
    stash_bfactors("a", ATOMS)
    stash_bfactors("b", ATOMS)
    clear_stash("a")
    assert not has_stash("a") and has_stash("b")
    clear_stash()
    assert not has_stash("b")


def test_restore_pushes_values_back_in_one_batch():
    """One alter, not one per atom — a 10k-atom structure would otherwise
    mean 10k round trips."""
    stash_bfactors("obj", ATOMS)
    port = FakePort()

    out = restore_bfactors(port, "obj")

    assert "Restored 3 B-factors" in out
    assert len(port.calls("alter")) == 1, port.call_log

    payload = next(c for c in port.commands if "stored.wiggles_b" in c)
    values = json.loads(payload.split("= ", 1)[1])
    # Keyed on model|index. Chain + residue + name + altloc collided on PDB
    # insertion codes and across a selection spanning two models, and a
    # collision restored one atom's B-factor onto another.
    assert values["m|1"] == 20.5
    assert values["m|2"] == 22.0
    assert values["m|3"] == 31.25


def test_restore_key_separates_atoms_that_share_everything_but_identity():
    """Altloc pairs, insertion codes and cross-model selections all collided
    on the old chain/residue/name/altloc key. Identity is (model, index)."""
    atoms = [
        Atom("A", "2", "SER", "CA", "A", 0.6, 10.0, "m", 4),
        Atom("A", "2", "SER", "CA", "B", 0.4, 90.0, "m", 5),
    ]
    stash_bfactors("obj", atoms)
    port = FakePort()
    restore_bfactors(port, "obj")

    payload = next(c for c in port.commands if "stored.wiggles_b" in c)
    values = json.loads(payload.split("= ", 1)[1])
    assert len(values) == 2
    assert values["m|4"] == 10.0
    assert values["m|5"] == 90.0


def test_the_old_key_would_have_collided_where_this_one_does_not():
    """Guards the guard: the atoms above differ only by altloc, and two
    residues differing only by an insertion code differ in nothing the old key
    looked at. Identity has to come from somewhere the model cannot repeat."""
    insertion = [
        Atom("A", "100", "SER", "CA", "", 0.3, 11.0, "m", 7),
        Atom("A", "100A", "SER", "CA", "", 0.9, 22.0, "m", 8),
    ]
    old_keys = {(a.chain, a.resi, a.name, a.alt) for a in insertion}
    assert len({a.key for a in insertion}) == 2
    # The old key does separate these two only because PyMOL happens to fold
    # the insertion code into `resi`; across two models it does not.
    two_models = [
        Atom("A", "1", "MET", "CA", "", 0.4, 11.0, "first", 1),
        Atom("A", "1", "MET", "CA", "", 0.8, 22.0, "second", 1),
    ]
    assert len({(a.chain, a.resi, a.name, a.alt) for a in two_models}) == 1, old_keys
    assert len({a.key for a in two_models}) == 2


def test_restore_falls_back_to_current_value_for_unknown_atoms():
    """The alter expression must default to `b`, not 0 — an atom added since
    the stash was taken should keep its value, not be zeroed."""
    stash_bfactors("obj", ATOMS)
    port = FakePort()
    restore_bfactors(port, "obj")

    (_, expression), _ = port.calls("alter")[0]
    assert expression.rstrip().endswith(", b)"), expression


def test_restore_without_a_stash_is_an_error():
    """Restoring from an empty stash would silently zero the column."""
    with pytest.raises(PortError, match="no saved B-factors"):
        restore_bfactors(FakePort(), "obj")


def test_restore_after_clear_is_an_error():
    stash_bfactors("obj", ATOMS)
    clear_stash("obj")
    with pytest.raises(PortError, match="no saved B-factors"):
        restore_bfactors(FakePort(), "obj")


def test_preservation_note_is_actionable():
    """A claim the user cannot act on is not worth making."""
    note = preservation_note("myobj", 42)
    assert "42" in note
    assert "restore_bfactors" in note
    assert "myobj" in note


# ── the first stash wins (MCPymol PR #57) ──────────────────────────────────


def test_a_second_stash_does_not_overwrite_the_first():
    """Every caller reads `b` *after* an earlier view may have overwritten it.

    So occupancy_view followed by qscore_view saved the occupancies as though
    they were the user's crystallographic data, and restore_bfactors then wrote
    those back and reported success for having destroyed the values it exists
    to protect. The stash is the only copy — nothing else holds them.
    """
    original = [Atom("A", "1", "MET", "CA", "", 1.0, 20.5, "m", 1)]
    overwritten = [Atom("A", "1", "MET", "CA", "", 1.0, 0.42, "m", 1)]

    assert stash_bfactors("obj", original) == 1
    assert stash_bfactors("obj", overwritten) == 1, "the count must stay truthful"

    port = FakePort()
    restore_bfactors(port, "obj")
    payload = next(c for c in port.commands if "stored.wiggles_b" in c)
    assert json.loads(payload.split("= ", 1)[1]) == {"m|1": 20.5}


def test_restoring_clears_the_stash_so_a_later_view_can_take_a_baseline():
    """Without this, first-stash-wins pins an object to B-factors that are no
    longer the ones worth saving."""
    stash_bfactors("obj", [Atom("A", "1", "MET", "CA", "", 1.0, 20.5, "m", 1)])
    restore_bfactors(FakePort(), "obj")
    assert not has_stash("obj")

    stash_bfactors("obj", [Atom("A", "1", "MET", "CA", "", 1.0, 99.0, "m", 1)])
    port = FakePort()
    restore_bfactors(port, "obj")
    payload = next(c for c in port.commands if "stored.wiggles_b" in c)
    assert json.loads(payload.split("= ", 1)[1]) == {"m|1": 99.0}


class TestADestroyedColumnIsNotStashed:
    """REVIEW #1. `preserve_bfactors=False` left _STASH empty, and
    first-stash-wins protects nothing when there is no first stash.

    "Nothing saved" and "nothing left to save" look identical in a dict and
    mean opposite things: the first invites a stash, the second forbids one.
    So the loss is recorded separately.
    """

    def test_a_stash_on_a_destroyed_object_is_refused(self):
        mark_bfactors_destroyed("obj")

        assert stash_bfactors("obj", ATOMS) == 0
        assert not has_stash("obj")

    def test_a_refused_stash_leaves_nothing_for_restore_to_write_back(self):
        """The whole point: restore must not put a view's scalar into the
        B-factor column and report success."""
        mark_bfactors_destroyed("obj")
        stash_bfactors("obj", ATOMS)

        with pytest.raises(PortError, match="no saved B-factors"):
            restore_bfactors(FakePort(), "obj")

    def test_an_undestroyed_object_still_stashes_normally(self):
        """The guard must not stop the ordinary path, which is every other
        view in the package."""
        assert stash_bfactors("obj", ATOMS) == len(ATOMS)
        assert has_stash("obj")

    def test_the_mark_is_per_object(self):
        """One destroyed object must not block preserving a different one."""
        mark_bfactors_destroyed("gone")

        assert stash_bfactors("gone", ATOMS) == 0
        assert stash_bfactors("kept", ATOMS) == len(ATOMS)

    def test_clearing_forgets_the_mark_too(self):
        """A reloaded structure's column is genuinely the user's again, so a
        stale mark would refuse to preserve data that is worth preserving."""
        mark_bfactors_destroyed("obj")
        clear_stash("obj")

        assert not bfactors_destroyed("obj")
        assert stash_bfactors("obj", ATOMS) == len(ATOMS)

    def test_the_two_notes_make_opposite_claims(self):
        """A user who saw a preservation note on an earlier view must not read
        the second one as the same guarantee."""
        preserved = preservation_note("obj", 2)
        destroyed = destroyed_note("obj")

        assert "are held" in preserved and "restore_bfactors" in preserved
        assert "are held" not in destroyed and "restore_bfactors" not in destroyed
        # The claim, not its capitalisation: the values are gone, and there is
        # a recovery route. Asserting the literal "Reload" made this fail when
        # the remedy was corrected to name clear_stash as well.
        assert "gone" in destroyed
        assert "reload" in destroyed.lower() and "clear_stash" in destroyed


class TestTheRemedyTheNotePrintsActuallyWorks:
    """`destroyed_note` tells the user how to recover. If that instruction does
    not restore the ability to preserve, the note is worse than silence: the
    user follows it, believes the session is healthy, and the *reloaded*
    crystallographic values are then overwritten with nothing saved.

    This package never loads structures — only volumes, via `load_map` and
    `load_ensemble` — so it cannot observe a reload and clear the mark itself.
    The instruction therefore has to name something the user can actually run.
    """

    def _remedy_calls(self, note: str) -> list[str]:
        """Function names the note tells the user to call."""
        import re

        return re.findall(r"\b([a-z_]+)\(", note)

    def test_the_note_names_a_call_that_clears_the_mark(self):
        mark_bfactors_destroyed("obj")
        note = destroyed_note("obj")

        named = self._remedy_calls(note)
        assert named, f"the note gives no runnable remedy at all: {note}"

        # Every function it names must exist and be importable from the package.
        import wiggles_em

        for name in named:
            assert hasattr(wiggles_em, name), (
                f"the note tells the user to call {name}(), which the package "
                f"does not export"
            )

    def test_following_the_remedy_restores_the_ability_to_preserve(self):
        """The end-to-end claim, not the wording."""
        mark_bfactors_destroyed("obj")
        assert stash_bfactors("obj", ATOMS) == 0, "precondition: stash refused"

        # What the note tells the user to do.
        clear_stash("obj")

        assert not bfactors_destroyed("obj")
        assert stash_bfactors("obj", ATOMS) == len(ATOMS), (
            "after following the note's instruction the object still cannot be "
            "preserved, so the remedy is inert"
        )

    def test_the_note_does_not_promise_a_reload_alone_is_enough(self):
        """A reload restores the column but not this package's belief about it,
        and the package has no way to notice — so 'reload' on its own is a
        promise it cannot keep."""
        note = destroyed_note("obj").lower()

        if "reload" in note:
            assert "clear_stash" in note, (
                "the note tells the user to reload, which this package cannot "
                f"observe, without naming the call that makes it true: {note}"
            )


def test_a_real_stash_survives_the_object_being_marked_destroyed():
    """`stash_bfactors` documents itself as "a no-op returning the size of the
    stash already held, so a caller's 'N values preserved' line stays true".

    The destroyed check ran first, so an object that had a genuine stash and was
    later marked destroyed reported 0 — and a direct caller (a host, not the
    backend, which is shielded by its own has_stash check) would print "The
    original 0 values are held" while restore_bfactors would in fact restore
    them correctly.
    """
    assert stash_bfactors("obj", ATOMS) == len(ATOMS)
    mark_bfactors_destroyed("obj")

    assert stash_bfactors("obj", ATOMS) == len(ATOMS), (
        "an existing stash was reported as empty because the object was later "
        "marked destroyed — but the stash is still there and still correct"
    )
    assert has_stash("obj")
    assert "Restored 3" in restore_bfactors(FakePort(), "obj")


def test_a_zero_return_means_nothing_held_not_object_destroyed():
    """The docstring invited `stash_bfactors(...) == 0` as a destruction probe,
    and after G6 that reading is wrong: a destroyed object with a real stash
    returns the stash size. `bfactors_destroyed` is the question to ask."""
    stash_bfactors("obj", ATOMS)
    mark_bfactors_destroyed("obj")

    assert stash_bfactors("obj", ATOMS) != 0
    assert bfactors_destroyed("obj"), "the destruction probe that still works"
