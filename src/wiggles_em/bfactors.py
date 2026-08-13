"""Saving and restoring B-factors around views that overwrite them.

Several views colour by pushing a scalar into the B-factor column, because
that is what PyMOL's ``spectrum`` reads. It is destructive, so anything doing
it owes the user a way back.

**The values are held here, in Python, not in PyMOL.** They come free: every
view already calls :func:`wiggles_em.atoms.fetch_atoms`, and that reads ``b``.
Keeping them on this side means a restore does not depend on some PyMOL field
surviving whatever the user does to the session in between.

Restoring pushes them back in **one** command via PyMOL's ``stored`` namespace
and a single ``alter``, rather than one ``alter`` per atom — a 10,000-atom
structure would otherwise mean 10,000 round trips. This mirrors what MCPymol's
``conservation_view`` already does in production, which is the reason to
prefer it over the ``custom`` field: ``custom`` is a *string* property, so a
float pushed into it would stringify, and a view claiming "B-factors preserved"
on that basis would be asserting something it had not established.
"""

from __future__ import annotations

import json

from wiggles_em.atoms import Atom
from wiggles_em.port import PortError, PymolPort, call

# obj name -> {(model, rank): b}  — Atom.key; rank, because index is
# renumbered by a removal and the stash would then repoint at another atom.
# Process-global. Correct for a single-session MCP server, which is what
# MCPymol is today; wrong the moment one process serves two sessions. Keyed
# by object name, so two sessions with the same object name would collide.
# Flagged rather than solved — see MOVING.md.
_STASH: dict[str, dict[tuple[str, str], float]] = {}

# Objects whose B-factor column was overwritten with **no stash taken**, so the
# originals are gone from this session. Tracked separately from _STASH because
# "nothing saved" and "nothing left to save" look identical in a dict and mean
# opposite things: the first invites a stash, the second forbids one.
#
# Without this, a view run with preserve_bfactors=False left _STASH empty, and
# first-stash-wins then offered no protection at all — the *next* view read a
# column already holding the first view's scalars, saved that as though it were
# the user's data, and restore_bfactors wrote it back reporting success.
_DESTROYED: set[str] = set()


def mark_bfactors_destroyed(obj: str) -> None:
    """Record that ``obj``'s B-factors were overwritten with nothing saved.

    Called by a backend that is about to write over the column while the caller
    has asked for no preservation. From this point a stash on ``obj`` is
    refused, because everything readable off the session is now some view's
    output rather than the user's data.
    """
    _DESTROYED.add(obj)


def bfactors_destroyed(obj: str) -> bool:
    """Were ``obj``'s B-factors overwritten with no stash taken?"""
    return obj in _DESTROYED


def _key(atom: Atom) -> tuple[str, str]:
    """Atom identity for restore purposes — see :attr:`Atom.key`.

    Was chain + residue + name + altloc, which collided on PDB insertion codes
    and across a selection spanning two models. Both collisions restored one
    atom's B-factor onto another, which is worse than losing it.
    """
    return atom.key


def stash_bfactors(obj: str, atoms: list[Atom]) -> int:
    """Record the B-factors of ``atoms`` under ``obj``. Returns how many.

    Takes atoms already fetched rather than querying again — the caller has
    them, and a second read would be a second round trip for data we hold.

    **The first stash wins.** Every caller reads ``b`` *after* some earlier
    view may already have overwritten it, so a second stash would save the
    first view's output as though it were the user's data — and
    :func:`restore_bfactors` would then write occupancies, or Q-scores, into
    the B-factor column and report success for having destroyed the
    crystallographic values. What is held here is the only copy. Re-stashing is
    a no-op returning the size of the stash already held, so a caller's "N
    values preserved" line stays true.

    **A destroyed column is not stashed, unless a real stash already exists.**
    If an earlier view overwrote ``obj``'s B-factors with no stash taken (see
    :func:`mark_bfactors_destroyed`), then what ``atoms`` carries is that view's
    scalar, not the user's data. Saving it would make
    :func:`restore_bfactors` write a Q-score or an occupancy into the B-factor
    column and report success — the exact outcome first-stash-wins exists to
    prevent, arrived at by the other door. Returns 0 in that case, and the
    caller says so rather than claiming a preservation it does not have.

    The two states are not exclusive: an object can be stashed by one view and
    then destroyed by a later one running with ``preserve_bfactors=False``. The
    stash is still the user's data and :func:`restore_bfactors` still puts it
    back correctly, so it wins and the true size is returned.

    A zero return therefore means "nothing is held", not "this object is
    destroyed" — use :func:`bfactors_destroyed` for that question.

    :func:`restore_bfactors` clears the stash, so a view run after a restore
    takes a fresh baseline.
    """
    # Order matters: an object can have a real stash AND be marked destroyed,
    # if a later view ran with preserve_bfactors=False. The stash is still the
    # user's data and restore_bfactors still puts it back correctly, so report
    # its true size — checking destroyed first said 0 and made the caller's
    # "N values preserved" line false.
    existing = _STASH.get(obj)
    if existing is not None:
        return len(existing)
    if obj in _DESTROYED:
        return 0
    _STASH[obj] = {_key(a): a.b for a in atoms}
    return len(_STASH[obj])


def has_stash(obj: str) -> bool:
    """Is there anything saved for ``obj``?"""
    return obj in _STASH


def stashed_count(obj: str) -> int:
    """How many values are held for ``obj``. 0 when nothing is.

    Exists so a caller can report what is held without re-stashing to find out
    — `stash_bfactors` returns the count, but calling it for the number alone
    would fetch atoms this caller does not need.
    """
    return len(_STASH.get(obj, ()))


def clear_stash(obj: str | None = None) -> None:
    """Forget the stash for ``obj``, or all of them.

    Clears the destroyed mark too: this means "forget what this session knows
    about that object's B-factors", and a stale mark would refuse to preserve a
    freshly reloaded structure whose column is genuinely the user's again.
    """
    if obj is None:
        _STASH.clear()
        _DESTROYED.clear()
    else:
        _STASH.pop(obj, None)
        _DESTROYED.discard(obj)


def restore_bfactors(port: PymolPort, obj: str) -> str:
    """Put ``obj``'s original B-factors back.

    Raises:
        PortError: nothing was stashed for ``obj`` — restoring from an empty
            stash would silently zero the column, which is worse than failing.
    """
    saved = _STASH.get(obj)
    if not saved:
        raise PortError(
            f"no saved B-factors for {obj!r}. Either no view has overwritten "
            f"them in this session, or the stash was cleared."
        )

    # JSON round-trip so the dict survives as a PyMOL command-line literal.
    # Tuple keys are not JSON-able, so key on a joined string and rebuild the
    # same string in the alter expression.
    flat = {"|".join(k): v for k, v in saved.items()}
    port.do(f"stored.wiggles_b = {json.dumps(flat)}")
    call(port, "alter", obj, "b=stored.wiggles_b.get('|'.join((model, str(rank))), b)")
    call(port, "rebuild", obj)
    # Drop the stash now the column holds those values again. Without this,
    # first-stash-wins would pin the object to B-factors that are no longer
    # the ones worth saving, and a view run after a restore could never record
    # a new baseline.
    clear_stash(obj)
    return (
        f"Restored {len(saved)} B-factors on {obj}. Atoms not present when the "
        f"stash was taken keep their current value."
    )


def preservation_note(obj: str, stashed: int) -> str:
    """The line a view prints about what it did to the B-factor column.

    Deliberately states where the values are and how to get them back, rather
    than the bare word "preserved" — a claim the user cannot act on is not
    worth making.
    """
    return (
        f"  B-factor column overwritten. The original {stashed} values are held "
        f"in this session;\n  call restore_bfactors(port, '{obj}') to put them "
        f"back."
    )


def destroyed_note(obj: str) -> str:
    """The line a view prints when the originals were already lost.

    Says what is actually true rather than staying silent. The alternative —
    printing nothing — leaves a user who has seen a preservation note on an
    earlier view assuming the same guarantee holds here.

    **The remedy has to be one the user can actually run.** An earlier version
    said only "reload the object". Reloading does restore the column, but this
    package never loads structures — only volumes, through :func:`load_map` and
    :func:`load_ensemble` — so it cannot observe the reload, the destroyed mark
    survives it, and the next view refuses to preserve a column that is now
    genuinely the user's again. Following that advice left someone worse off
    than ignoring it: they would believe the session was healthy while the
    freshly reloaded crystallographic values were overwritten with nothing
    saved. So the note names :func:`clear_stash`, which is what actually
    clears the mark.
    """
    return (
        f"  WARNING: {obj}'s B-factor column was already overwritten by an "
        f"earlier view\n  that preserved nothing, so the original values are "
        f"gone from this session.\n  Nothing was saved, because what is in the "
        f"column now is that view's output.\n  To recover: reload {obj}, then "
        f"call clear_stash('{obj}') so a later\n  view can preserve it again — "
        f"reloading alone does not tell this session\n  that the column is "
        f"trustworthy."
    )
