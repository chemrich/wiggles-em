"""Stepping through a latent trajectory, without claiming what it cannot.

This is the tool the compendium is most careful about, and the care is not
decorative. The Flatiron Institute blind challenge put 41 submissions on the
same thyroglobulin data and found the field recovers **motion** well and
**populations** badly: the simulated ground truth had three modes, most
submissions recovered two, and only three of the 41 found all three. The mode
they missed was the middle one. `limits` gives the information-theoretic reason
— measurement noise sets a floor on how finely conformation space can be
resolved, and a state sandwiched between two others is what that floor eats
first.

Two SPEC invariants fall out of that, and both are enforced here.

**I2 — no unlabelled latent plot.** Every rendering names its generating method
and carries that method's caveat. The caveats are not interchangeable: cryoDRGN
latent density can bear no relation to the truth, while 3DVA frames are linear
interpolations no particle ever occupied. When
:func:`~wiggles_em.heterogeneity.load_ensemble` could not identify the
method, this **refuses to render** rather than draw an unlabelled trajectory.

**I3 — a gap is not an absence.** Nothing here may state or imply that an
unpopulated region of latent space is unvisited. :data:`ABSENCE_CLAIMS` lists
the phrasings that would breach that, and a test asserts no output contains
one. There is deliberately **no scatter plot and no density estimate** in this
tool: the single most requested latent rendering is the one the evidence says
is most often wrong, and the honest version of it is not a better-drawn plot,
it is not drawing it.

**The trap this tool exists to close.** PyMOL normalises each map independently
on load, so contouring every frame at "1.5 sigma" contours each frame at *its
own* sigma. A traversal is precisely the case where the density genuinely
changes between frames — that is the signal — and per-frame normalisation
rescales it away. Frames drawn that way look reassuringly similar no matter
what the ensemble does. So the level is fixed in **absolute** terms and
converted to each frame's own sigma, and the spread of those sigmas is reported,
because that spread *is* the density change.
"""

from __future__ import annotations

from wiggles_em.density import (
    DEFAULT_SIGMA,
    MapStats,
    rms_meaning,
    to_absolute,
    to_sigma,
    usable_rms,
)
from wiggles_em.heterogeneity import Ensemble, Method, loaded_ensemble
from wiggles_em.port import PortError
from wiggles_em.provenance import provenance_banner
from wiggles_em.scene import (
    ColorFlat,
    Colour,
    Frames,
    Isosurface,
    Legend,
    Rep,
    Scene,
    Sel,
    Unit,
    resolve_colour,
)

#: Phrasings that would turn a gap into an absence claim. I3 forbids all of
#: them, in any output this package produces about latent space. Checked by
#: :func:`contains_absence_claim` and asserted over every tier-3 report.
ABSENCE_CLAIMS: tuple[str, ...] = (
    "unvisited",
    "never visited",
    "unpopulated",
    "not populated",
    "no molecules",
    "no particles here",
    "forbidden state",
    "forbidden region",
    "does not occur",
    "never occupied",
    "empty region means",
    "absent conformation",
)

GAP_LEGEND = (
    "A GAP IS NOT AN ABSENCE. An empty region between frames is a region of "
    "UNKNOWN OCCUPANCY, not a state the molecule avoids. In a blind challenge on "
    "41 submissions the commonest failure was missing an intermediate state that "
    "was genuinely present — usually the middle one of three — and there is an "
    "information-theoretic floor underneath that result, not just method error. "
    "Read this trajectory as: these conformations are supported. Read nothing "
    "into what lies between them."
)

#: What the viewer cannot tell the user, said out loud.
#:
#: The Flatiron challenge measured per-frame resolution against ground truth and
#: found it ranged from about 9 Å to about 18.5 Å — twofold across submissions,
#: and "a large range of values" *within* a single submission. This tool draws
#: every frame at one level in one colour, so a frame built from less data is
#: rendered exactly like a frame built from more.
#:
#: That is a viewer problem specifically: blur reads as floppiness, not as data
#: quality, so a soft frame looks like a conformational finding rather than a
#: reconstruction artefact.
#:
#: There is deliberately no number attached. Resolution is not in an MRC header;
#: the FSC that would measure it needs a reference no viewer session has; and a
#: score derived from the header statistics would look measured without being
#: measured, which is what `MapStats.stated` vs `.measured` and
#: `ensemble_spread_view`'s refusal to call spread an uncertainty both exist to
#: prevent. So this states the limit of the tool rather than a claim about the
#: caller's ensemble, and points at the job directory, where the answer does
#: live if it lives anywhere.
FRAME_QUALITY_LEGEND = (
    "THESE FRAMES ARE NOT NECESSARILY OF EQUAL QUALITY, and nothing in this render "
    "distinguishes them. In a community challenge, frame resolution within a single "
    "submission ranged from about 9 Å to about 18.5 Å, and across submissions it "
    "varied about twofold. This viewer cannot measure which frame is which: "
    "resolution is not recorded in a map header, and the Fourier shell correlation "
    "that would establish it needs a reference this session does not have. Read a "
    "smooth frame and a detailed frame as equally uncertain unless you know "
    "otherwise from the job that produced them."
)

POPULATION_LEGEND = (
    "NO POPULATION IS SHOWN, deliberately. Cryo-EM heterogeneity methods recover "
    "motion reliably and relative populations unreliably; those are different "
    "claims and a viewer conflates them constantly. This tool renders the motion. "
    "It draws no latent scatter and estimates no density, because the rendering "
    "users most want from latent space is the one the evidence says is most often "
    "wrong."
)


def contains_absence_claim(text: str) -> str | None:
    """The first I3-breaching phrase in ``text``, or None.

    Public so the invariant can be asserted over any report this package grows
    later, not only the ones that exist today.
    """
    lowered = text.lower()
    return next((phrase for phrase in ABSENCE_CLAIMS if phrase in lowered), None)


def frame_levels(ensemble: Ensemble, absolute: float) -> list[float | None]:
    """One sigma level per frame for a shared absolute contour.

    ``None`` for a frame whose header cannot define a sigma scale — a zero
    RMS, or MRC's ``rms=-1`` "statistics not computed" marker. Reported rather
    than silently substituted, because a substituted level would draw a frame
    at a contour nobody chose, and the negative case divides cleanly enough to
    produce a plausible wrong one.
    """
    levels: list[float | None] = []
    for header in ensemble.headers:
        try:
            levels.append(to_sigma(MapStats.stated(header), absolute))
        except ValueError:
            levels.append(None)
    return levels


def latent_traverse_view(
    ensemble_name: str,
    *,
    level: float | None = None,
    units: str = "sigma",
    name: str | None = None,
    color: Colour = "skyblue",
    build_movie: bool = True,
) -> tuple[str, Scene]:
    """Build a steppable isosurface trajectory through an ensemble's frames.

    Args:
        ensemble_name: An ensemble loaded through ``load_ensemble``.
        level: Contour level. ``None`` uses :data:`DEFAULT_SIGMA` against the
            first frame **whose header carries a usable RMS**, then holds that
            **absolute** value across all frames. The report names which frame
            that was; it is frame 1 unless frame 1's statistics are missing.
        units: ``"sigma"`` interpreted against that same anchor frame, or
            ``"absolute"``, which needs no anchor at all.
        name: Prefix for the isosurfaces. Defaults to ``<ensemble>_surf``.
        color: Colour for every frame. One colour on purpose — a per-frame
            spectrum would encode frame index as if it were a measured quantity.
        build_movie: Wire the frames to PyMOL's movie timeline so they can be
            stepped and played.

    Returns:
        A report: the method and its caveat, the shared absolute level with the
        per-frame sigmas it corresponds to, and the I2/I3 legends.

    Raises:
        PortError: no such ensemble, or its method was never identified (I2).
        ValueError: ``units`` is not recognised, the level cannot be
            converted, or ``color`` is a name
            :func:`~wiggles_em.scene.resolve_colour` has no value for.
    """
    if units not in ("sigma", "absolute"):
        raise ValueError(f"units must be 'sigma' or 'absolute', got {units!r}")

    ensemble = loaded_ensemble(ensemble_name)
    if ensemble is None:
        raise PortError(
            f"{ensemble_name!r} was not loaded through load_ensemble, so its method, "
            f"frame order and headers are unavailable. Every one of those is needed "
            f"here: the method to satisfy invariant I2, the order to make the "
            f"trajectory mean anything, and the headers to hold one contour level "
            f"across frames PyMOL normalised separately."
        )

    # I2, the refusal. An unlabelled latent rendering is the thing the
    # invariant exists to prevent, so this is the tool working, not failing.
    if ensemble.method is Method.UNKNOWN:
        return "\n".join(
            [
                f"latent_traverse_view({ensemble_name})",
                "",
                "  REFUSED: this ensemble's generating method was never identified,",
                "  so no interpretive caveat can be attached to what would be drawn.",
                "",
                "  SPEC invariant I2 — no unlabelled latent plot. The caveats are not",
                "  interchangeable and they are not decoration: cryoDRGN's latent",
                "  density can bear no relation to the truth, while 3DVA's frames are",
                "  linear interpolations no particle ever occupied. Rendering these",
                "  identically would assert something for one that holds only for the",
                "  other, and the picture looks equally convincing either way.",
                "",
                "  " + "\n  ".join(ensemble.evidence),
                "",
                "  If you know what produced this directory, say so and reload:",
                f"      load_ensemble(path, {ensemble_name!r}, method='cryodrgn')",
                "  Accepted: " + ", ".join(m.value for m in Method if m is not Method.UNKNOWN),
            ]
        ), Scene()

    no_usable_rms = PortError(
        f"no frame of {ensemble_name!r} has a usable RMS in its header, so an "
        f"absolute contour level cannot be converted to the sigma a viewer "
        f"contours in. The headers may be stale or the maps unnormalised."
    )

    # Anchor the absolute level on the first header that can define a sigma
    # scale, not simply on frame 0. A header carrying rms=0, or MRC's rms=-1
    # "statistics not computed" marker, converts nothing — and taking it
    # anyway yielded an anchor of dmean, a number unrelated to the level asked
    # for, which was then applied to every other frame. A level already given
    # in absolute units needs no anchor at all.
    # Which frame this is has to travel with it. The report used to say "the
    # FIRST frame" while the anchor had moved to the first frame with a usable
    # RMS — so a reader checking the contour against frame 1's header, as the
    # report told them to, got a different number and concluded the tool had
    # converted wrongly.
    anchor_position = next((i for i, h in enumerate(ensemble.headers) if usable_rms(h)), None)
    anchor = ensemble.headers[anchor_position] if anchor_position is not None else None
    anchor_frame = None if anchor_position is None else anchor_position + 1
    # The anchor is used whenever the level is in sigma — whether it was
    # defaulted or supplied. Attributing it only on the defaulted path left a
    # caller-supplied sigma silently interpreted against a frame they never
    # named, which is the same defect the FIRST-frame claim was.
    anchor_used = level is None or units == "sigma"

    if level is None:
        if anchor is None:
            raise no_usable_rms
        absolute = to_absolute(MapStats.stated(anchor), DEFAULT_SIGMA)
        used_default = True
    elif units == "absolute":
        absolute = float(level)
        used_default = False
    else:
        if anchor is None:
            raise no_usable_rms
        absolute = to_absolute(MapStats.stated(anchor), float(level))
        used_default = False

    levels = frame_levels(ensemble, absolute)
    # Carry the frame's own 1-based number, not its position among the
    # survivors. Numbering over the filtered list makes `_03` hold frame 4's
    # density the moment a frame is dropped for rms=0, and nothing downstream
    # can tell — which is the off-by-one heterogeneity._natural_key exists to
    # prevent, reintroduced one layer up.
    usable = [
        (number, obj, sigma)
        for number, (obj, sigma) in enumerate(zip(ensemble.objects, levels, strict=True), start=1)
        if sigma is not None
    ]
    skipped_frames = [number for number, sigma in enumerate(levels, start=1) if sigma is None]
    if not usable:
        raise PortError(
            f"no frame of {ensemble_name!r} has a usable RMS in its header, so an "
            f"absolute contour level cannot be converted to the sigma PyMOL "
            f"contours in. The headers may be stale or the maps unnormalised."
        )

    prefix = name or f"{ensemble_name}_surf"
    width = max(2, len(str(ensemble.n_frames)))
    surfaces: list[str] = []
    numbers: list[int] = []
    ops: list = []
    for number, obj, sigma in usable:
        surface = f"{prefix}_{number:0{width}d}"
        # Each frame carries its OWN sigma, converted from one absolute level
        # against that frame's header. A fixed sigma across frames contours
        # each against its own normalisation, which flattens away the density
        # change the traversal exists to show.
        # `equivalent` carries the shared absolute level. A host with
        # normalize_ccp4_maps off needs the contour in absolute units, and an
        # ensemble frame has no load_map record for the backend to convert
        # against — so without this the whole traversal refused to render.
        ops.append(
            Isosurface(
                surface,
                obj,
                level=sigma,
                unit=Unit.SIGMA,
                style=Rep.SURFACE,
                equivalent=absolute,
            )
        )
        # One colour for every frame, on purpose: a per-frame spectrum would
        # encode frame index as if it were a measured quantity.
        ops.append(ColorFlat(Sel.obj(surface), resolve_colour(color)))
        surfaces.append(surface)
        # Carried, not re-derived. The report promises "a surface's number is
        # always the frame it was made from"; that promise is only keepable if
        # the number travels with the surface.
        numbers.append(number)

    if len(surfaces) > 1:
        ops.append(Frames(tuple(surfaces), tuple(numbers), build_timeline=build_movie))
    ops.append(Legend(POPULATION_LEGEND))
    ops.append(Legend(GAP_LEGEND))
    ops.append(Legend(FRAME_QUALITY_LEGEND))

    return _report(
        ensemble,
        prefix,
        surfaces,
        absolute,
        levels,
        used_default,
        build_movie,
        skipped_frames,
        anchor_frame,
        anchor_used,
    ), Scene(ops)


def _skipped_reasons(ensemble: Ensemble, skipped_frames: list[int]) -> list[str]:
    """One line per distinct rms among the skipped frames, saying what it means.

    Grouped by value rather than listed per frame: twenty frames from one
    unprocessed run all carry the same sentinel, and twenty identical lines
    bury the one that differs.
    """
    by_value: dict[float, list[int]] = {}
    for number in skipped_frames:
        header = ensemble.headers[number - 1]
        by_value.setdefault(header.rms, []).append(number)
    return [
        f"frame{'s' if len(frames) > 1 else ''} "
        f"{', '.join(str(n) for n in frames)}: rms={value:g} "
        f"({rms_meaning(value, brief=True)})"
        for value, frames in sorted(by_value.items())
    ]


def _report(
    ensemble: Ensemble,
    prefix: str,
    surfaces: list[str],
    absolute: float,
    levels: list[float | None],
    used_default: bool,
    build_movie: bool,
    skipped_frames: list[int],
    anchor_frame: int | None,
    anchor_used: bool,
) -> str:
    known = [s for s in levels if s is not None]
    lo, hi = min(known), max(known)
    skipped = len(levels) - len(known)

    lines = [
        f"latent_traverse_view({ensemble.name})",
        "",
        f"  Method: {ensemble.method.label}",
        f"  {ensemble.method.caveat}",
        "",
        f"  {len(surfaces)} frames as `{surfaces[0]}` … `{surfaces[-1]}`, one colour.",
    ]
    if build_movie and len(surfaces) > 1:
        lines += [
            "  Wired to the movie timeline: `mplay` to run it, `frame N` to step.",
        ]

    lines += [
        "",
        f"  Contour: {absolute:.6g} absolute, held constant across every frame.",
        f"  That is {lo:.3g}–{hi:.3g} sigma depending on the frame.",
    ]

    spread = (hi - lo) / max(abs(hi), abs(lo), 1e-12)
    if spread > 0.05:
        lines += [
            "",
            f"  Those sigmas differ by {spread * 100:.0f}% across the trajectory, and",
            "  that difference IS the density change this traversal shows. PyMOL",
            "  normalises every map independently on load, so contouring each frame",
            "  at a fixed sigma would have contoured each at its own scale and",
            "  flattened exactly this signal — the frames would look reassuringly",
            "  alike whatever the ensemble does.",
        ]
    else:
        lines += [
            "",
            "  Those sigmas barely differ, so the frames carry near-identical",
            "  density statistics. The level is still held in absolute terms, which",
            "  is what makes that a finding rather than an artefact of normalising",
            "  each frame separately.",
        ]

    if anchor_used:
        lines += [
            "",
            *(
                [
                    f"  No level given, so {DEFAULT_SIGMA} sigma against frame {anchor_frame}",
                    "  was used and converted to an absolute value. A generic starting",
                    "  point, not a recommendation for this ensemble.",
                ]
                if used_default
                else [
                    f"  The level you gave is in sigma, read against frame {anchor_frame}",
                    "  and converted to the absolute value above.",
                ]
            ),
            # Which frame the sigma was read against is a fact the reader needs
            # to check the number, and it is not frame 1 whenever frame 1's
            # header carries no usable rms.
            *(
                []
                if anchor_frame == 1
                else [
                    f"  Frame {anchor_frame}, not frame 1: every earlier frame reports an",
                    "  rms that cannot define a sigma scale, so no sigma level exists",
                    "  against them to convert.",
                ]
            ),
        ]

    if skipped:
        # Name them. "1 frame(s) were skipped" plus a gap in the numbering is
        # not enough for a reader to work out which density belongs to which
        # latent coordinate, and guessing is how the traversal gets misread.
        which = ", ".join(f"frame {n}" for n in skipped_frames)
        lines += [
            "",
            f"  {skipped} frame(s) were skipped: their headers report an rms that",
            "  cannot define a sigma scale, so no level exists to draw them at. They",
            "  are loaded but not contoured.",
            # What the rms actually is, and what that value means. "rms=0" was
            # hard-coded here, so a frame whose statistics were never computed
            # (MRC writes -1) was reported as a flat map.
            *[f"    {reason}" for reason in _skipped_reasons(ensemble, skipped_frames)],
            f"  Skipped: {which}. The surfaces keep their own frame numbers, so the",
            "  sequence has gaps rather than being renumbered — a surface's number",
            "  is always the frame it was made from.",
        ]

    lines += [
        "",
        provenance_banner(ensemble.objects[0]),
        "",
        "  " + POPULATION_LEGEND,
        "",
        "  " + GAP_LEGEND,
        "",
        "  " + FRAME_QUALITY_LEGEND,
    ]
    return "\n".join(lines)
