"""Occupancy views — **sense 1 only**, per-atom crystallographic ``q``.

These tools read the ``occupancy`` column of a coordinate file and nothing
else. They never infer, estimate or display compositional occupancy — what
fraction of imaged particles contain a subunit. That is sense 2, it lives in
method output rather than the model, and it belongs to ``composition_view`` in
tier 3.

The distinction is the single most important design decision in SPEC.md, and
the reason is that the failure is invisible: a model can be ``q = 1.0``
everywhere while its subunit is present in half the particles. Both statements
are true, they answer different questions, and a render that conflates them
looks entirely normal. So every legend these functions emit names the sense.

See the Wiggles compendium entry `occupancy-two-senses`.
"""

from __future__ import annotations

from wiggles_em.atoms import Atom, altloc_groups, group_by_residue
from wiggles_em.scene import (
    ColorByScalar,
    ColorFlat,
    Hide,
    Label,
    Legend,
    Rep,
    ScalarField,
    Scene,
    Sel,
    Sense,
    Show,
)

# Occupancy at or above this is "full" for reporting purposes. Occupancies are
# stored as float32 and refined, so exact 1.0 comparison mislabels ordinary
# fully-occupied atoms.
FULL_OCCUPANCY = 0.999

# qFit removes conformers below 10% occupancy during refinement, so a
# multiconformer model's conformer count is partly a threshold artifact. Worth
# saying when low-occupancy atoms are present.
QFIT_FLOOR = 0.10

# Distinct colours for altloc groups, in assignment order.
_ALTLOC_COLOURS = (
    "skyblue",
    "salmon",
    "palegreen",
    "wheat",
    "lightpink",
    "paleyellow",
    "lightblue",
    "lightorange",
)

SENSE_1_LEGEND = (
    "Occupancy shown is SENSE 1: per-atom crystallographic occupancy (q) read "
    "from the coordinate file. It is NOT compositional occupancy — the fraction "
    "of imaged particles containing this subunit. A model can be fully occupied "
    "here and half-present in that sense, simultaneously and correctly."
)


def _stats(atoms: list[Atom]) -> dict[str, float | int]:
    qs = [a.q for a in atoms]
    partial = [a for a in atoms if a.q < FULL_OCCUPANCY]
    return {
        "n_atoms": len(atoms),
        "n_partial": len(partial),
        "n_below_qfit_floor": sum(1 for a in atoms if 0.0 < a.q < QFIT_FLOOR),
        "q_min": min(qs),
        "q_max": max(qs),
    }


def occupancy_view(atoms: list[Atom], obj: str) -> tuple[str, Scene]:
    """Colour ``obj`` by per-atom occupancy, de-emphasising partial atoms.

    If nothing in the object is partially occupied, this says so plainly and
    draws nothing. A spectrum across a constant value is a rainbow that means
    nothing, and a viewer that renders one is inviting a false reading.

    The domain is a fixed ``(0.0, 1.0)`` and not the observed range, which is
    the whole reason :class:`~wiggles_em.scene.ColorByScalar` demands one: a
    model sitting between 0.95 and 1.0 stretched over its own range becomes a
    full rainbow implying variation that is not there.

    How the scalar reaches the viewer is a backend's business. Both route
    per-atom values through the B-factor column, but PyMOL has one copy of the
    object and must stash the originals, while protean builds a display copy —
    so the "your B-factors were overwritten" caveat comes from the backend,
    not from here.

    Args:
        atoms: Every atom in ``obj``, already read.
        obj: Object or selection name, for the scene's selections and the report.

    Returns:
        A report ending with the sense-1 legend, and the scene to draw.
    """
    st = _stats(atoms)
    lines = [f"occupancy_view({obj})", ""]
    legend = Legend(SENSE_1_LEGEND, sense=Sense.ATOM_OCCUPANCY)

    if st["n_partial"] == 0:
        lines += [
            f"  All {st['n_atoms']} atoms are fully occupied (q >= {FULL_OCCUPANCY}).",
            "  Nothing to show: no colouring applied.",
            "",
            "  A spectrum over a constant value would be a rainbow that means",
            "  nothing. If you expected partial occupancy here, check that the",
            "  model actually carries alternate conformations — see altloc_view.",
            "",
            SENSE_1_LEGEND,
        ]
        # A scene of nothing but a legend does not draw — see Scene.draws. The
        # refusal is the result, and it has to be distinguishable from a view
        # that drew something.
        return "\n".join(lines), Scene([legend])

    target = Sel.obj(obj)
    field = ScalarField.per_atom([(a.key, a.q) for a in atoms])
    scene = Scene(
        [
            ColorByScalar(target, field, domain=(0.0, 1.0), palette="red_white_blue"),
            # Partial atoms shown as sticks so they read as "modelled alternative"
            # rather than as the single truth. Note this must NOT exclude protein:
            # side-chain alternates are the main thing this view exists to show.
            Show(target & Sel.lt("q", FULL_OCCUPANCY), Rep.STICKS),
            legend,
        ]
    )

    lines += [
        f"  {st['n_partial']} of {st['n_atoms']} atoms partially occupied (q < {FULL_OCCUPANCY}).",
        f"  q range: {st['q_min']:.3f} – {st['q_max']:.3f}",
        "  Coloured red (q=0) → white (q=0.5) → blue (q=1), fixed scale.",
        "",
    ]
    if st["n_below_qfit_floor"]:
        lines += [
            f"  {st['n_below_qfit_floor']} atoms sit below q={QFIT_FLOOR:g}. qFit removes",
            "  conformers under that during refinement, so their presence here",
            "  suggests a different refinement protocol — worth checking.",
            "",
        ]

    lines.append(SENSE_1_LEGEND)
    return "\n".join(lines), scene


def altloc_view(atoms: list[Atom], obj: str, *, label: bool = True) -> tuple[str, Scene]:
    """Show every alternate conformation at once, one colour per altloc group.

    Occupancies go into the labels rather than being implied by the colouring,
    because the occupancy *is* the measurement — drawing all altlocs in equal
    weight discards it, which is the commonest way this view is got wrong.

    Args:
        atoms: Every atom in ``obj``, already read.
        obj: Object or selection name.
        label: Label one atom per alternate residue with its altloc and occupancy.

    Returns:
        A report ending with the sense-1 legend, and the scene to draw.
    """
    groups = altloc_groups(atoms)
    lines = [f"altloc_view({obj})", ""]
    legend = Legend(SENSE_1_LEGEND, sense=Sense.ATOM_OCCUPANCY)

    if not groups:
        lines += [
            f"  No alternate conformations: all {len(atoms)} atoms have a blank altloc.",
            "  Nothing to show.",
            "",
            "  This is an ordinary single-conformer model. It says nothing about",
            "  whether the molecule is heterogeneous — only that this file does",
            "  not model it.",
            "",
            SENSE_1_LEGEND,
        ]
        return "\n".join(lines), Scene([legend])

    target = Sel.obj(obj)
    ops = [
        Hide(target, Rep.EVERYTHING),
        Show(target, Rep.CARTOON),
        ColorFlat(target, "grey70"),
    ]

    per_group: list[str] = []
    for i, alt in enumerate(groups):
        colour = _ALTLOC_COLOURS[i % len(_ALTLOC_COLOURS)]
        sel = target & Sel.prop("alt", alt)
        ops += [Show(sel, Rep.STICKS), ColorFlat(sel, colour)]
        members = [a for a in atoms if a.alt == alt]
        occ = sorted({round(a.q, 3) for a in members})
        occ_text = ", ".join(f"{q:g}" for q in occ[:4]) + (" …" if len(occ) > 4 else "")
        per_group.append(
            f"  alt {alt}: {colour:<12} {len(members):>5} atoms   occupancy {occ_text}"
        )

    if label:
        # One label per alternate residue — labelling every atom is unreadable.
        for alt in groups:
            ops.append(
                Label(
                    target & Sel.prop("alt", alt) & Sel.prop("name", "CA"),
                    f"{alt} %.2f",
                    fields=("q",),
                )
            )

    ops.append(legend)
    scene = Scene(ops)
    residues = group_by_residue([a for a in atoms if a.alt.strip()])
    lines += [
        f"  {len(groups)} altloc group(s) across {len(residues)} residue(s).",
        "",
        *per_group,
        "",
        "  Backbone in grey; alternates as coloured sticks.",
    ]
    if label:
        lines.append("  CA atoms labelled with altloc and occupancy.")
    lines += [
        "",
        "  Occupancies are shown rather than encoded in the colouring: drawing",
        "  every alternate at equal weight would discard the measurement.",
        "",
        SENSE_1_LEGEND,
    ]
    return "\n".join(lines), scene
