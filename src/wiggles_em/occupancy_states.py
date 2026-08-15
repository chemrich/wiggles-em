"""Showing how much of each conformational state there is — sense 3.

This is the view the package was reorganised around. It exists because state
occupancies **are** recoverable from cryo-EM: what fails is one family of
estimator, not the measurement. See :mod:`wiggles_em.populations` for the
evidence and the taxonomy.

**A third sense of "occupancy", and it stays separate from the other two
forever.** ``occupancy_view`` reads per-atom crystallographic ``q``:
*how much of this atom is there*. ``composition_view`` reads compositional
occupancy: *is this subunit present in the particle*. This reads state
occupancy: *which conformation is the particle in*. A complex can be fully
intact — sense 2 saturated — and still be split evenly between two
conformations. Different question, different unit of analysis, different
estimator. Every legend here names :attr:`~wiggles_em.scene.Sense.STATE_OCCUPANCY`
so the three can never be read as each other.

## Two rendering decisions, and the first one is the important one

**Occupancy is never encoded as visibility.** The obvious rendering is to fade
each state in proportion to how populated it is, and it is exactly wrong here.
The sparsely populated intermediate is the state these methods have historically
*lost* — it is the entire reason this part of the field exists — and drawing it
at 5% opacity reproduces that loss in the picture. A viewer would look at a
correct result and see the same nothing that a broken estimator gave. So every
state is drawn equally visible, and the number is carried by colour and by the
report instead.

**A ramp is only drawn for weights that are quantitative.** Colour along a scale
asserts that the differences between the numbers mean something. When the
weights came from counting particles into classes or histogramming a latent
space, they do not — Evans *et al.* show that estimator returning a nearly flat
distribution for a genuinely trimodal one. In that case every state gets one
flat colour, the reported values still appear in the text as *what the method
said*, and the legend says why they are not on a scale. The check is
:attr:`~wiggles_em.populations.Populations.is_quantitative`, which is decided by
the declared source and never by inspecting the values.
"""

from __future__ import annotations

from wiggles_em.density import DEFAULT_SIGMA, MapStats, to_absolute, usable_rms
from wiggles_em.heterogeneity import Method, loaded_ensemble
from wiggles_em.latent import FRAME_QUALITY_LEGEND, frame_levels
from wiggles_em.populations import Populations
from wiggles_em.port import PortError
from wiggles_em.provenance import provenance_banner
from wiggles_em.scene import (
    BLUE_WHITE_RED,
    ColorFlat,
    Colour,
    Isosurface,
    Legend,
    Rep,
    Scene,
    Sel,
    Sense,
    Unit,
    resolve_colour,
)

SENSE_3_LEGEND = (
    "STATE OCCUPANCY (sense 3): the fraction of imaged particles in each "
    "CONFORMATION. Not per-atom occupancy `q`, which says how much of an atom is "
    "there, and not compositional occupancy, which says whether a subunit is "
    "present at all. A complex can be 100% intact and still be split evenly "
    "between two conformations — different question, different answer."
)

VISIBILITY_LEGEND = (
    "EVERY STATE IS DRAWN EQUALLY VISIBLE, deliberately. Occupancy is shown by "
    "colour and stated in the report; it is NOT shown by transparency. Fading a "
    "rare state in proportion to its population would make the sparsely "
    "populated intermediate almost invisible — and that state is precisely the "
    "one these methods have historically lost, so drawing it that way would "
    "reproduce the failure in the picture rather than reveal it."
)

FLAT_COLOUR_LEGEND = (
    "NO COLOUR SCALE IS DRAWN, because these weights are not quantitative. A "
    "ramp would assert that the differences between them mean something. The "
    "reported values are shown as what the method returned, not as measurements. "
    "Recomputing them by deconvolution or ensemble reweighting is what would put "
    "them on a scale."
)


def _ramp_colour(fraction: float) -> tuple[float, float, float]:
    """Interpolate the blue-white-red stops at ``fraction`` in 0..1.

    Blue is the least populated end and red the most, which is the direction
    ``BLUE_WHITE_RED`` already carries as a value — the backend does not have
    to know the convention.
    """
    stops = BLUE_WHITE_RED
    if len(stops) == 1:
        return stops[0]
    position = max(0.0, min(1.0, fraction)) * (len(stops) - 1)
    low = min(int(position), len(stops) - 2)
    t = position - low
    a, b = stops[low], stops[low + 1]
    return tuple(x + (y - x) * t for x, y in zip(a, b, strict=True))  # type: ignore[return-value]


def state_occupancy_view(
    ensemble_name: str,
    populations: Populations,
    *,
    level: float | None = None,
    units: str = "sigma",
    name: str | None = None,
    flat_color: Colour = "skyblue",
) -> tuple[str, Scene]:
    """Draw each conformational state, with how much of it there is.

    Args:
        ensemble_name: An ensemble loaded through ``load_ensemble``. Its volumes
            are the states, in order.
        populations: One weight per state, carrying how it was computed. Build
            it with :meth:`~wiggles_em.populations.Populations.declare` or read
            it with :func:`~wiggles_em.recovar.read_deconvolved_weights`.
        level: Contour level, held **absolute** across states so the density
            differences between them survive per-map normalisation. ``None``
            takes :data:`~wiggles_em.density.DEFAULT_SIGMA` against the first
            state with usable statistics.
        units: ``"sigma"`` against that anchor state, or ``"absolute"``.
        name: Prefix for the surfaces. Defaults to ``<ensemble>_state``.
        flat_color: The single colour used when the weights are not
            quantitative and no ramp may be drawn.

    Returns:
        A report and a Scene. The report carries the occupancy of each state,
        the source's caveat, and relative free energies when a temperature was
        declared.

    Raises:
        PortError: no such ensemble, its method was never identified (I2), or
            it has no state with usable map statistics.
        ValueError: ``units`` is unrecognised, or the number of weights does
            not match the number of states.
    """
    if units not in ("sigma", "absolute"):
        raise ValueError(f"units must be 'sigma' or 'absolute', got {units!r}")

    ensemble = loaded_ensemble(ensemble_name)
    if ensemble is None:
        raise PortError(
            f"{ensemble_name!r} was not loaded through load_ensemble, so its states, "
            f"their order and their map statistics are unavailable — and the order is "
            f"what pairs a weight with a state."
        )

    # I2 applies here as much as to a traversal: this is a rendering of a
    # method's output and it carries that method's interpretive caveat.
    if ensemble.method is Method.UNKNOWN:
        return "\n".join(
            [
                f"state_occupancy_view({ensemble_name})",
                "",
                "  REFUSED: this ensemble's generating method was never identified, so",
                "  no interpretive caveat can be attached to the occupancies drawn.",
                "",
                "  SPEC invariant I2 — no unlabelled latent plot. An occupancy is a",
                "  stronger claim than a shape, not a weaker one, so the rule applies",
                "  with more force here rather than less.",
                "",
                "  " + "\n  ".join(ensemble.evidence),
                "",
                f"      load_ensemble(path, {ensemble_name!r}, method='recovar')",
            ]
        ), Scene()

    if populations.n_states != ensemble.n_frames:
        raise ValueError(
            f"{populations.n_states} weights for {ensemble.n_frames} states. These are "
            f"paired by position, so a mismatch does not mean 'some states are "
            f"unweighted' — it means every pairing after the first difference is "
            f"wrong, and the picture would look entirely normal."
        )

    anchor_at = next((i for i, h in enumerate(ensemble.headers) if usable_rms(h)), None)
    if anchor_at is None and (level is None or units == "sigma"):
        raise PortError(
            f"no state of {ensemble_name!r} has a usable RMS in its header, so a sigma "
            f"level cannot be converted. Pass units='absolute' to contour directly."
        )
    if level is None:
        absolute = to_absolute(MapStats.stated(ensemble.headers[anchor_at]), DEFAULT_SIGMA)  # type: ignore[index]
    elif units == "absolute":
        absolute = float(level)
    else:
        absolute = to_absolute(MapStats.stated(ensemble.headers[anchor_at]), float(level))  # type: ignore[index]

    levels = frame_levels(ensemble, absolute)
    prefix = name or f"{ensemble_name}_state"
    width = max(2, len(str(ensemble.n_frames)))

    quantitative = populations.is_quantitative
    biggest = max(populations.probabilities)
    ops: list = []
    drawn: list[int] = []
    skipped: list[int] = []

    for index, (obj, sigma) in enumerate(zip(ensemble.objects, levels, strict=True)):
        number = index + 1
        if sigma is None:
            skipped.append(number)
            continue
        surface = f"{prefix}_{number:0{width}d}"
        ops.append(
            Isosurface(
                surface, obj, level=sigma, unit=Unit.SIGMA, style=Rep.SURFACE, equivalent=absolute
            )
        )
        if quantitative:
            # Scaled against the most populated state so the ramp spans the
            # data. Never against a fixed 0..1: with three states near a third
            # each, a fixed domain would render them indistinguishable.
            share = populations.probabilities[index] / biggest if biggest > 0 else 0.0
            ops.append(ColorFlat(Sel.obj(surface), _ramp_colour(share)))
        else:
            ops.append(ColorFlat(Sel.obj(surface), resolve_colour(flat_color)))
        drawn.append(number)

    if not drawn:
        raise PortError(
            f"no state of {ensemble_name!r} could be contoured: every header reports an "
            f"rms that cannot define a sigma scale."
        )

    ops.append(Legend(SENSE_3_LEGEND, sense=Sense.STATE_OCCUPANCY))
    ops.append(Legend(VISIBILITY_LEGEND, sense=Sense.STATE_OCCUPANCY))
    if not quantitative:
        ops.append(Legend(FLAT_COLOUR_LEGEND, sense=Sense.STATE_OCCUPANCY))
    ops.append(Legend(FRAME_QUALITY_LEGEND))

    return _report(
        ensemble_name,
        ensemble,
        populations,
        prefix,
        drawn,
        skipped,
        absolute,
        quantitative,
    ), Scene(ops)


def _report(
    ensemble_name: str,
    ensemble,
    populations: Populations,
    prefix: str,
    drawn: list[int],
    skipped: list[int],
    absolute: float,
    quantitative: bool,
) -> str:
    lines = [
        f"state_occupancy_view({ensemble_name})",
        "",
        f"  Method: {ensemble.method.label}",
        f"  {ensemble.method.caveat}",
        "",
        f"  {len(drawn)} state(s) as `{prefix}_*`, all drawn equally visible.",
        f"  Contour: {absolute:.6g} absolute, held constant across states.",
        "",
        "  State occupancies:",
    ]
    for number in drawn:
        share = populations.probabilities[number - 1]
        bar = "#" * max(1, round(share * 40)) if share > 0 else ""
        lines.append(f"    state {number:>3}  {share * 100:6.2f}%  {bar}")

    if skipped:
        which = ", ".join(str(n) for n in skipped)
        lines += [
            "",
            f"  {len(skipped)} state(s) not contoured — headers report an rms that",
            f"  cannot define a sigma scale: {which}. Their weights are still listed",
            "  above only for the states that were drawn, so the percentages shown do",
            "  not sum to 100%.",
        ]

    lines += ["", populations.banner()]

    if populations.temperature_k is not None:
        lines += [
            "",
            f"  Relative free energy at {populations.temperature_k:g} K, kJ/mol, "
            f"against state {drawn[0]}:",
        ]
        for number, (delta_g, minus, plus) in zip(
            range(1, populations.n_states + 1),
            populations.relative_free_energy(reference=drawn[0] - 1),
            strict=True,
        ):
            if number not in drawn:
                continue
            if minus is None and plus is None:
                lines.append(f"    state {number:>3}  {delta_g:+8.2f}")
            elif plus is None:
                lines.append(f"    state {number:>3}  {delta_g:+8.2f}  -{minus:.2f} / +unbounded")
            else:
                lines.append(f"    state {number:>3}  {delta_g:+8.2f}  -{minus:.2f} / +{plus:.2f}")
        lines += [
            "  An unbounded upper error means the state's uncertainty reaches its own",
            "  population: it may be arbitrarily unfavourable. That is a property of",
            "  taking a logarithm near zero, not a failure of the fit.",
        ]

    lines += [
        "",
        provenance_banner(ensemble.objects[0]),
        "",
        "  " + SENSE_3_LEGEND,
        "",
        "  " + VISIBILITY_LEGEND,
    ]
    if not quantitative:
        lines += ["", "  " + FLAT_COLOUR_LEGEND]
    return "\n".join(lines)
