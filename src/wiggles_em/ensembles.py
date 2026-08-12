"""Multi-state ensembles: spread, and morphing that refuses when ill-posed.

An ensemble and a multiconformer model both say "not one state", and they look
similar on disk, but they encode different claims. Altloc says *at this site,
these discrete alternatives, in these proportions*. Multi-model says *here are
N draws from a distribution* — the spread is the variation, no member is
privileged, and there are no occupancies at all.

These two views are for the multi-model case. Occupancy lives in
:mod:`wiggles_em.occupancy`. See the Wiggles compendium entry `multiconformer`.

No numpy: the arithmetic is a mean and a root-mean-square, and keeping tier 1
dependency-free is worth more than the vectorisation.
"""

from __future__ import annotations

import math

from wiggles_em.atoms import Atom
from wiggles_em.port import PortError
from wiggles_em.scene import (
    ColorByScalar,
    Legend,
    Morph,
    ScalarField,
    Scene,
    SceneOp,
    Sel,
    SizeByScalar,
)

#: How far the positional spread may exceed the rigid-invariant radial spread
#: before the difference is called rigid motion.
#:
#: High, and the reason is worth keeping. A first attempt at 3.0 fired on a
#: genuine hinge: any *asymmetric* conformational change moves the centroid,
#: so the positional spread picks up that shift on its own and runs a few
#: times the radial figure without a gram of rigid-body motion. Truly rigid
#: states have a radial spread of exactly zero — translation and rotation both
#: leave every atom's distance to its own centroid untouched — so the rigid
#: case sits at infinity and a wide threshold separates the two cleanly.
#:
#: This detects rigid motion that *dominates*, not rigid motion that is
#: present. A modest drift under a lot of real flexing will not trip it, and
#: should not: the report already says whether the states were fitted.
RIGID_RATIO = 10.0

SPREAD_LEGEND = (
    "Spread is the RMS deviation of each atom's position across states. It is a "
    "description of how much the deposited members differ, NOT a calibrated "
    "uncertainty and NOT an error bar. The number of members and the refinement "
    "protocol both shape it, so a wider tube means these models disagree more — "
    "it does not mean the true position is less well determined."
)


def _rmsf(points: list[tuple[float, float, float]]) -> float:
    """RMS deviation of a set of positions about their centroid."""
    n = len(points)
    cx = sum(p[0] for p in points) / n
    cy = sum(p[1] for p in points) / n
    cz = sum(p[2] for p in points) / n
    total = sum((p[0] - cx) ** 2 + (p[1] - cy) ** 2 + (p[2] - cz) ** 2 for p in points)
    return math.sqrt(total / n)


def per_atom_spread(
    coords_by_state: list[list[tuple[float, float, float]]],
) -> list[float]:
    """RMSF per atom across states.

    Raises:
        PortError: states disagree on atom count, which makes the whole
            calculation meaningless rather than merely approximate.
    """
    if len(coords_by_state) < 2:
        raise PortError("need at least two states to compute spread")
    counts = {len(c) for c in coords_by_state}
    if len(counts) != 1:
        raise PortError(
            f"states have differing atom counts {sorted(counts)} — cannot pair "
            f"atoms across states, so spread is undefined"
        )
    n_atoms = counts.pop()
    return [_rmsf([state[i] for state in coords_by_state]) for i in range(n_atoms)]


def per_residue_spread(atoms: list[Atom], atom_spread: list[float]) -> dict[tuple[str, str], float]:
    """Mean atom spread within each residue, keyed by (chain, resi)."""
    if len(atoms) != len(atom_spread):
        raise PortError(
            f"{len(atoms)} atoms but {len(atom_spread)} spread values — the "
            f"iterate and coordinate orders have diverged"
        )
    acc: dict[tuple[str, str], list[float]] = {}
    for atom, value in zip(atoms, atom_spread, strict=True):
        acc.setdefault(atom.residue, []).append(value)
    return {key: sum(v) / len(v) for key, v in acc.items()}


def radial_spread(coords_by_state: list[list[tuple[float, float, float]]]) -> list[float]:
    """Per-atom spread of the distance to that state's own centroid.

    Invariant under both translation and rotation, which is what makes it a
    check on :func:`per_atom_spread` rather than another version of it. An
    internally rigid molecule tumbling across its box has a large positional
    spread and a radial spread of zero; a genuinely flexing one has both.

    It is deliberately *not* offered as a flexibility measure: it is blind to
    tangential motion, so a residue swinging around the centre at constant
    radius scores zero. It answers one question — is the positional number
    dominated by rigid motion — and nothing else.
    """
    radii_by_state = []
    for state in coords_by_state:
        n = len(state)
        cx = sum(p[0] for p in state) / n
        cy = sum(p[1] for p in state) / n
        cz = sum(p[2] for p in state) / n
        radii_by_state.append(
            [math.dist(p, (cx, cy, cz)) for p in state]
        )
    return [
        _rmsf([(radii[i], 0.0, 0.0) for radii in radii_by_state])
        for i in range(len(radii_by_state[0]))
    ]


def ensemble_spread_view(
    atoms: list[Atom],
    coords: list[list[tuple[float, float, float]]],
    obj: str,
    *,
    as_putty: bool = True,
    superposed: bool,
) -> tuple[str, Scene]:
    """Colour and thicken ``obj`` by how much its states disagree, per residue.

    Unlike occupancy, this domain **is** data-derived — ``(0, max spread)``,
    scaled to this object — and that is why the domain is a parameter rather
    than a convention. Spread has no absolute scale: 2 Å is enormous in a
    rigid core and unremarkable in a flexible loop, so a fixed domain would
    make most ensembles look uniformly blue. The report states the maximum so
    the scaling is never implicit.

    **Spread is only a conformational quantity once the states share a frame.**
    It is the RMS of each atom about its own position across states, which
    measures whatever separates those states — including a rigid-body offset
    that says nothing about flexibility. A viewer loads a multi-model PDB's
    models as states without aligning them, so an MD trajectory that drifted
    across its box, or two independently refined ensemble members 2 Å apart,
    would report that offset at *every* residue and paint an internally rigid
    molecule red end to end.

    Fitting belongs to the host, which has the session and a superposition
    routine already — PyMOL's ``intra_fit``, biotite's ``superimpose``. This
    view takes ``superposed`` as a statement of what was done.

    But a flag is a claim, and a wrong claim is silent, so the numbers are
    checked against it: :func:`radial_spread` is invariant under translation
    and rotation, and a positional spread that dwarfs it means rigid motion
    whatever the caller said. The report leads with that.

    Args:
        atoms: Every atom in ``obj``, already read.
        coords: Coordinates per state, outer list one entry per state.
        obj: The object name, for the scene's selections and the report.
        as_putty: Also vary tube width by spread, not only colour.
        superposed: Whether the host fitted the states onto a common frame
            before reading these coordinates. Required, with no default: a
            default would be this view guessing at the one thing that decides
            whether its number means anything.

    Returns:
        A report ending with the spread legend, and the scene to draw.
    """
    n_states = len(coords)
    if n_states < 2:
        return (
            f"ensemble_spread_view({obj})\n\n"
            f"  Object has {n_states} state; spread needs at least 2.\n"
            f"  Nothing to show.\n\n"
            f"  If this is a multiconformer model rather than an ensemble, the\n"
            f"  heterogeneity is in altloc groups, not states — use altloc_view."
        ), Scene([Legend(SPREAD_LEGEND)])

    atom_spread = per_atom_spread(coords)
    residue_spread = per_residue_spread(atoms, atom_spread)

    # The check on the caller's claim. Comparing means rather than maxima
    # because a single flexible loop should not mask a whole-body drift.
    positional = sum(atom_spread) / len(atom_spread)
    radial = sum(radial_spread(coords)) / len(atom_spread)
    rigid_dominated = positional > RIGID_RATIO * max(radial, 1e-9)

    values = sorted(residue_spread.values())
    lo, hi = values[0], values[-1]
    median = values[len(values) // 2]

    target = Sel.obj(obj)
    field = ScalarField.per_residue(list(residue_spread.items()))
    domain = (0.0, round(hi, 4))
    ops: list[SceneOp] = [ColorByScalar(target, field, domain=domain, palette="blue_white_red")]
    if as_putty:
        ops.append(SizeByScalar(target, field, domain=domain))
    ops.append(Legend(SPREAD_LEGEND))

    widest = sorted(residue_spread.items(), key=lambda kv: -kv[1])[:5]
    widest_text = ", ".join(f"{c}/{r} {v:.2f}Å" for (c, r), v in widest)

    return "\n".join(
        [
            f"ensemble_spread_view({obj})",
            "",
            f"  {n_states} states, {len(atoms)} atoms, {len(residue_spread)} residues.",
            f"  Per-residue spread: {lo:.2f} – {hi:.2f} Å (median {median:.2f} Å)",
            f"  Most variable: {widest_text}",
            "",
            f"  Coloured blue (rigid) → red ({hi:.2f} Å), scaled to this object.",
            (
                "  States were superposed before measuring."
                if superposed
                else "  States were NOT superposed — spread includes any rigid-body"
                "\n  offset between them."
            ),
            *(
                [
                    "",
                    "  ! RIGID-BODY MOTION DOMINATES. Each atom's distance to its own",
                    f"  state's centroid varies by {radial:.2f} Å on average, while its",
                    f"  position varies by {positional:.2f} Å. That gap is a rigid",
                    "  offset or rotation between the states, not flexibility — the",
                    "  radial figure is invariant to both. Fit the states onto a",
                    "  common frame and measure again; the number above is not a",
                    "  conformational quantity.",
                ]
                if rigid_dominated
                else []
            ),
            "" if not as_putty else "  Tube width also tracks spread (putty).",
            "",
            SPREAD_LEGEND,
        ]
    ), Scene(ops)


def morph_states(
    counts: list[int],
    obj: str,
    *,
    name: str | None = None,
    steps: int = 30,
    validate_only: bool = False,
) -> tuple[str, Scene]:
    """Interpolate between states, refusing when interpolation is ill-posed.

    The value here is the refusal, not the interpolation — PyMOL morphs
    natively. Morphing is only meaningful when all states share a topology, so
    that atom *i* in state 1 and atom *i* in state 2 are the same atom. That
    holds for a deformation-model ensemble (3DFlex, DynaMight), where every
    state is the same reference warped. It does **not** hold for volumes
    reconstructed independently and modelled separately, and a morph across
    those is an animation of a correspondence that was never established.

    The topology check is the whole judgement, and it is arithmetic over the
    per-state atom counts — no viewer required. Whether the viewer can then
    perform the interpolation is a separate question the backend answers;
    ``cmd.morph`` is Incentive-only, so on open-source PyMOL the request is
    sound and unhonourable at once.

    Args:
        counts: Atom count per state, one entry per state.
        obj: A multi-state object.
        name: Name for the morph object. Defaults to ``<obj>_morph``.
        steps: Frames to interpolate.
        validate_only: Check topology and report without creating the morph.

    Returns:
        A report and a scene. Refusals explain what differs rather than just
        declining, and return a scene that draws nothing.
    """
    n_states = len(counts)
    if n_states < 2:
        return (
            f"morph_states({obj})\n\n"
            f"  REFUSED: object has {n_states} state; morphing needs at least 2."
        ), Scene()

    if len(set(counts)) != 1:
        detail = ", ".join(f"state {i + 1}: {c}" for i, c in enumerate(counts))
        return "\n".join(
            [
                f"morph_states({obj})",
                "",
                "  REFUSED: states differ in atom count, so atoms cannot be paired",
                "  across them and interpolation has no defined meaning.",
                "",
                f"  {detail}",
                "",
                "  A morph here would animate a correspondence that was never",
                "  established. If these came from independently reconstructed",
                "  volumes, that is expected — the ensemble is a set of separate",
                "  models, not one model deformed. Compare them with",
                "  ensemble_spread_view instead, which needs no pairing.",
            ]
        ), Scene()

    lines = [
        f"morph_states({obj})",
        "",
        f"  Topology check passed: {n_states} states, {counts[0]} atoms each.",
        "  Atoms pair across states, so interpolation is well posed.",
        "",
    ]
    if validate_only:
        lines.append("  validate_only=True — no morph created.")
        return "\n".join(lines), Scene()

    morph_name = name or f"{obj}_morph"
    lines += [
        f"  Requested `{morph_name}` with {steps} interpolated frames.",
        "",
        "  Interpolated frames are generated, not observed. Only the states you",
        "  started from correspond to anything that was reconstructed; every",
        "  frame between them is this tool's invention.",
    ]
    return "\n".join(lines), Scene([Morph(morph_name, obj, steps=steps)])
