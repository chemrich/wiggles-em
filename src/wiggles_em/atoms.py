"""Reading per-atom properties out of PyMOL, and the residue key they group by.

Shared by every tier-1 view. Kept separate from the views so that the one
place which knows how an atom is fetched is not also the place that decides
how it is coloured.
"""

from __future__ import annotations

from dataclasses import dataclass

from wiggles_em.port import ITERATE_TO_LIST, PortError, PymolPort

# The iterate expression tier-1 views read. Order matters: it defines the
# tuple layout that _to_atom unpacks.
#
# ``rank``, not ``index``: see :attr:`Atom.key`. PyMOL renumbers ``index`` when
# atoms are removed, which silently repoints every stashed value past the
# deletion at a different atom.
ATOM_EXPR = "chain, resi, resn, name, alt, q, b, model, rank"


@dataclass(frozen=True)
class Atom:
    """One atom's identity and the two scalars tier 1 cares about."""

    chain: str
    resi: str
    resn: str
    name: str
    alt: str
    q: float  # occupancy — SENSE 1, per-atom. Never compositional.
    b: float
    # No defaults. A default would make every atom built without them share
    # one key, which is the collision this field pair exists to remove.
    model: str
    rank: int

    @property
    def residue(self) -> tuple[str, str]:
        """Chain + residue number. The grouping key for per-residue views."""
        return (self.chain, self.resi)

    @property
    def key(self) -> tuple[str, str]:
        """A unique identity for this atom: its model and its rank within it.

        Chain + residue + name + altloc is **not** unique, and both ways it
        fails are ordinary rather than exotic:

        - a PDB insertion code (residues 100 and 100A in one chain) — antibody
          models are full of them;
        - a selection spanning two loaded structures, where chain, residue and
          atom name repeat by construction.

        Either collision silently hands one atom another atom's value, which a
        fixed colour ramp then renders as a perfectly ordinary result.

        **``rank``, not ``index``.** Both are unique within an object, so both
        remove the collision — but ``index`` is *renumbered when atoms are
        removed*, and a stash keyed on it silently repoints at a different atom
        the moment anyone tidies a structure. ``rank`` is the atom's position in
        the file as loaded and does not move.

        Checked against a live PyMOL 3.1.0 rather than taken from the
        documentation, because the first attempt at this used ``index`` on the
        strength of its uniqueness alone:

        - ``remove hydro`` renumbered ``index`` (6→4, 7→5, 8→6) and left every
          ``rank`` alone;
        - so did removing a residue from the middle (5→3, 6→4);
        - ``remove solvent`` renumbered *nothing*, because PyMOL sorts solvent
          to the end of index order — which is why the finding's own repro did
          not reproduce, and why this was checked by observation instead;
        - ``rank`` came back unique across the object, and stays
          non-contiguous after a removal (0, 2, 4, 5…), which is exactly what
          makes it stable.
        """
        return (self.model, str(self.rank))


def _to_atom(row: object) -> Atom:
    if not isinstance(row, (list, tuple)) or len(row) != 9:
        raise PortError(
            f"malformed atom row {row!r}: expected 9 fields ({ATOM_EXPR}), "
            f"got {len(row) if isinstance(row, (list, tuple)) else type(row).__name__}"
        )
    try:
        return Atom(
            chain=str(row[0]),
            resi=str(row[1]),
            resn=str(row[2]),
            name=str(row[3]),
            alt=str(row[4]),
            q=float(row[5]),
            b=float(row[6]),
            model=str(row[7]),
            rank=int(row[8]),
        )
    except (TypeError, ValueError) as exc:
        raise PortError(f"malformed atom row {row!r}: {exc}") from exc


def fetch_atoms(port: PymolPort, selection: str) -> list[Atom]:
    """Read every atom in ``selection``.

    Raises:
        PortError: the selection is empty, or the port returned something
            that is not a list of atom rows.

    An empty selection is an error rather than an empty list on purpose. Every
    caller here is about to render a view, and a view of nothing is worse than
    a message saying the selection matched nothing — see MCPymol issue #15,
    where an empty result was reported as success.
    """
    raw = port.query(ITERATE_TO_LIST, selection, ATOM_EXPR)
    if raw is None:
        raise PortError(f"selection {selection!r} returned nothing")
    if not isinstance(raw, (list, tuple)):
        raise PortError(f"expected a list of atom rows, got {type(raw).__name__}")
    if not raw:
        raise PortError(f"selection {selection!r} matched no atoms")
    return [_to_atom(row) for row in raw]


def count_states(port: PymolPort, obj: str) -> int:
    """How many states (models) the object has."""
    raw = port.query("count_states", obj)
    try:
        return int(raw)  # type: ignore[arg-type]
    except (TypeError, ValueError) as exc:
        raise PortError(f"count_states({obj!r}) returned {raw!r}") from exc


def fetch_state_coords(
    port: PymolPort, selection: str, state: int
) -> list[tuple[float, float, float]]:
    """Coordinates of ``selection`` in one state, in atom order."""
    raw = port.query("get_coords", selection, state)
    if not isinstance(raw, (list, tuple)) or not raw:
        raise PortError(f"no coordinates for {selection!r} in state {state}")
    out: list[tuple[float, float, float]] = []
    for point in raw:
        try:
            x, y, z = point  # type: ignore[misc]
        except (TypeError, ValueError) as exc:
            raise PortError(f"malformed coordinate {point!r}") from exc
        out.append((float(x), float(y), float(z)))
    return out


def group_by_residue(atoms: list[Atom]) -> dict[tuple[str, str], list[Atom]]:
    """Group atoms by (chain, resi), preserving first-seen order."""
    out: dict[tuple[str, str], list[Atom]] = {}
    for atom in atoms:
        out.setdefault(atom.residue, []).append(atom)
    return out


def altloc_groups(atoms: list[Atom]) -> list[str]:
    """Distinct non-blank altloc identifiers, sorted.

    A blank altloc means "no alternate here" and is deliberately excluded —
    it is not a group, it is the absence of grouping.
    """
    return sorted({a.alt for a in atoms if a.alt.strip()})
