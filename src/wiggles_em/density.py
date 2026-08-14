"""Showing density around a selection, with the contour level stated honestly.

Every viewer makes users hunt for a contour threshold. EMDB has been publishing
the depositor's own recommended level all along — ``map.contour_list[].level``
with ``source: AUTHOR`` — and nothing surfaces it. See
``REPOSITORIES.md``.

**The unit trap this module exists to close.** PyMOL normalises MRC/CCP4 maps
on load (``normalize_ccp4_maps``, on by default), so an ``isomesh`` level is in
**sigma**. EMDB's author contour is an **absolute** map value. They are not
interchangeable, and confusing them is not a subtle error:

======================  =================  ==================================
Entry                   Author contour     As a PyMOL level if used directly
======================  =================  ==================================
EMD-30913               0.05 absolute      0.05 sigma — contours the noise
EMD-11638               0.116 absolute     0.116 sigma — contours the noise
======================  =================  ==================================

Converted properly those are **3.16 sigma** and **4.42 sigma**. The conversion
needs the map's mean and RMS, which the MRC header carries, which is why
:func:`density_view` requires the map to have been loaded through
:func:`wiggles_em.maps.load_map`: without the header it cannot state a level in
both units, and stating one without saying which is how the mistake happens.

**The header's statistics are a claim, not a measurement**, and the converters
take a :class:`MapStats` rather than a header so a host that has walked the
voxels can convert against what it measured. `density_view` here uses
``MapStats.stated(header)`` — the file's own numbers — because that is all it
has; the alternative exists for hosts that have more.
"""

from __future__ import annotations

from dataclasses import dataclass
from enum import Enum

from wiggles_em.mapinfo import MapHeader
from wiggles_em.maps import loaded_map
from wiggles_em.port import PortError
from wiggles_em.provenance import provenance_banner
from wiggles_em.scene import Isosurface, Legend, Rep, Scene, Sel, Unit

#: Used when no level is given and no author contour is available. A generic
#: starting point, not a recommendation for any particular map — the report
#: says so rather than letting it read as considered.
DEFAULT_SIGMA = 1.5

#: Default carve radius (Å) around the selection.
DEFAULT_CARVE = 2.0


def usable_rms(header: MapHeader) -> bool:
    """Can this header's RMS define a sigma scale?

    Only a strictly positive RMS can, and two non-positive values occur in
    real files. Zero means a flat or unnormalised map. **Negative is the
    dangerous one:** MRC2014 writes ``rms=-1`` for "statistics not computed",
    which is what mrcfile leaves behind without ``update_header_stats()``, and
    it divides perfectly cleanly. ``to_sigma(0.05)`` returned ``-2.05``, and
    ``local_resolution_view`` turned ascending Ångström breakpoints into a
    descending ramp — blue bound to the worst-resolved density, under a legend
    saying blue was the best.

    Testing ``if not header.rms`` catches only the first of the two.

    Kept for callers that hold a header and nothing else — it is the same rule
    as :attr:`MapStats.usable` and delegates to it, so the two cannot drift
    into disagreeing about what a usable RMS is.
    """
    return MapStats.stated(header).usable


def rms_meaning(rms: float, *, brief: bool = False) -> str:
    """What a non-positive RMS in a real header usually means.

    Public because three modules print it. It was private and `latent` imported
    it anyway, which is the sort of thing that gets duplicated rather than
    shared the next time — and a second copy of this is a second place for
    "rms=0 means flat" to be wrong.

    ``brief`` gives the parenthetical form, for listing several frames at once
    where the full sentence would run four times the width of the surrounding
    report. Same knowledge, one definition, two lengths.
    """
    if rms < 0:
        return (
            "statistics never computed"
            if brief
            else "a negative RMS marks statistics that were never computed, which "
            "is what mrcfile leaves behind without update_header_stats()"
        )
    return (
        "flat map, or stale header statistics"
        if brief
        else "a zero RMS means a flat map, or one whose header statistics are stale"
    )


class StatsSource(str, Enum):
    """Where a sigma scale's two numbers came from.

    Not decoration. A viewer that reports header fields as though it had
    measured them is the failure this exists to keep visible — see
    :class:`MapStats`.
    """

    #: The file's own header said so. A claim, and claims go stale.
    STATED = "stated"
    #: Computed from the voxels.
    MEASURED = "measured"


@dataclass(frozen=True)
class MapStats:
    """The mean and RMS a sigma scale needs, and where they came from.

    These functions used to take a :class:`MapHeader`, which meant they could
    only ever convert against the *file's claim* about itself. That is usually
    fine and occasionally very wrong: a map that has been cropped or rescaled
    keeps whatever header nobody updated, and the arithmetic proceeds cleanly
    over numbers that no longer describe the data.

    **The finding that forced this**, from protean: Mol\\*'s ``grid.stats`` are,
    for CCP4/MRC, stored header fields passed straight through unexamined. A
    test fixture written with deliberately false statistics failed reporting the
    header's fake minimum **with the dimensions correct** — so the volume had
    genuinely parsed, and every number describing it was the file's claim. A
    host that walks the voxels holds better numbers and, before this, had
    nowhere to put them.

    How far it reaches: Mol\\*'s default isosurface is ``relative: 2``, computed
    as ``relativeValue * sigma + mean`` against those same statistics. A stale
    header therefore misplaces the surface in *any* viewer that trusts it, and
    the result looks entirely ordinary.

    So the scale carries its :class:`StatsSource`, for the same reason
    :class:`~wiggles_em.scene.Unit` exists: a level may not be a bare number,
    and a scale may not be two bare numbers.
    """

    mean: float
    rms: float
    source: StatsSource

    @classmethod
    def stated(cls, header: MapHeader) -> MapStats:
        """The file's own claim. Named rather than implicit, so that taking it
        is a visible choice — it is the right one whenever nothing has measured
        the voxels, which is most of the time."""
        return cls(mean=header.dmean, rms=header.rms, source=StatsSource.STATED)

    @classmethod
    def measured(cls, *, mean: float, rms: float) -> MapStats:
        """Statistics computed from the data. Keyword-only: ``mean`` and
        ``rms`` are both floats in the same range and transposing them is
        silent."""
        return cls(mean=mean, rms=rms, source=StatsSource.MEASURED)

    @property
    def usable(self) -> bool:
        """Can these numbers define a sigma scale? See :func:`usable_rms` for
        why a *negative* RMS is the dangerous case rather than zero."""
        return self.rms > 0


def to_sigma(stats: MapStats, absolute: float) -> float:
    """Convert an absolute map value to sigma units.

    Takes a :class:`MapStats` rather than a :class:`MapHeader` so that a host
    which has measured the voxels can convert against what it measured. Pass
    ``MapStats.stated(header)`` to use the file's own claim.

    Raises:
        ValueError: the RMS cannot define a sigma scale — see :func:`usable_rms`.
    """
    if not stats.usable:
        raise ValueError(
            f"{stats.source.value} statistics report rms={stats.rms:g}, so sigma "
            f"is undefined and an absolute level cannot be converted: "
            f"{rms_meaning(stats.rms)}. Give the level in absolute units "
            f"instead, or supply measured statistics if you have them."
        )
    return (absolute - stats.mean) / stats.rms


def to_absolute(stats: MapStats, sigma: float) -> float:
    """Convert a sigma level back to an absolute map value.

    Raises:
        ValueError: the RMS cannot define a sigma scale, so there is no sigma
            to convert *from* — see :func:`usable_rms`. Returning
            ``mean + sigma * rms`` on a negative RMS hands back a number of the
            right shape and the wrong sign, which is the harder failure.
    """
    if not stats.usable:
        raise ValueError(
            f"{stats.source.value} statistics report rms={stats.rms:g}, so sigma "
            f"is undefined and a sigma level has no absolute equivalent: "
            f"{rms_meaning(stats.rms)}."
        )
    return stats.mean + sigma * stats.rms


def density_view(
    map_obj: str,
    selection: str,
    *,
    level: float | None = None,
    units: str = "sigma",
    carve: float = DEFAULT_CARVE,
    name: str | None = None,
) -> tuple[str, Scene]:
    """Draw an isomesh around ``selection``, reporting the level in both units.

    The scene carries the level in the unit the *caller* gave, tagged with
    which one it is, and the backend converts if its viewer needs the other.
    Both units go in the report regardless, because that is the only way a
    reader can tell an EMDB author contour from a sigma level — the confusion
    this view exists to prevent.

    Args:
        map_obj: A volume loaded through :func:`wiggles_em.maps.load_map`.
        selection: What to carve the mesh around.
        level: Contour level. ``None`` uses :data:`DEFAULT_SIGMA` and says so.
        units: ``"sigma"`` (PyMOL's units) or ``"absolute"`` (EMDB's). An
            absolute level is converted, and both values are reported.
        carve: Carve radius in Å around the selection.
        name: Mesh object name. Defaults to ``<map_obj>_mesh``.

    Returns:
        A report: the level in sigma *and* absolute units, the sigma value
        itself, the provenance banner, and — where the map is an EMDB
        deposition and no level was given — a pointer to the author contour.

    Raises:
        PortError: ``map_obj`` was not loaded through ``load_map``, so its
            header is unavailable and the level cannot be stated in both units.
        ValueError: ``units`` is not recognised, or conversion is impossible.
    """
    if units not in ("sigma", "absolute"):
        raise ValueError(f"units must be 'sigma' or 'absolute', got {units!r}")

    record = loaded_map(map_obj)
    if record is None:
        raise PortError(
            f"{map_obj!r} was not loaded through load_map, so its header is "
            f"unavailable. Without it a contour level cannot be stated in both "
            f"sigma and absolute units, and stating one without saying which is "
            f"exactly how an EMDB author contour gets used as a sigma level. "
            f"Load it with load_map first."
        )

    header = record.header
    used_default = level is None

    if used_default:
        sigma = DEFAULT_SIGMA
        absolute: float | None = None
    elif units == "absolute":
        sigma = to_sigma(MapStats.stated(header), float(level))  # type: ignore[arg-type]
        absolute = float(level)  # type: ignore[arg-type]
    else:
        sigma = float(level)  # type: ignore[arg-type]
        absolute = None

    if absolute is None and usable_rms(header):
        absolute = to_absolute(MapStats.stated(header), sigma)

    mesh = name or f"{map_obj}_mesh"
    # Emitted in sigma because that is what was resolved above for the report;
    # tagging the unit is what stops a backend guessing. A viewer that wants
    # absolute converts, rather than every view remembering to.
    scene = Scene(
        [
            Isosurface(
                mesh,
                map_obj,
                level=sigma,
                unit=Unit.SIGMA,
                style=Rep.MESH,
                carve_around=Sel.raw(selection),
                carve_radius=carve,
            ),
            Legend(provenance_banner(map_obj), provenance=map_obj),
        ]
    )

    # Reports the rms it actually saw, and what that value means. Hard-coding
    # "rms=0" told a reader with MRC's "statistics not computed" sentinel that
    # their map was flat — a different problem with a different remedy.
    absolute_text = (
        f"{absolute:.6g}"
        if absolute is not None
        else f"unknown (rms={header.rms:g}: {rms_meaning(header.rms)})"
    )
    lines = [
        f"density_view({map_obj} around {selection})",
        "",
        f"  Contour: {sigma:.3g} sigma  =  {absolute_text} absolute",
        f"  Map sigma (header rms): {header.rms:.6g}   mean: {header.dmean:.6g}",
        f"  Mesh `{mesh}`, carved {carve:g} Å around the selection.",
        # The registry is keyed by object name, so a map deleted and reloaded
        # under the same name keeps the old header and nothing can detect it —
        # a viewer does not expose the file an object came from. Naming the
        # path is what makes a substitution visible to a reader.
        f"  Header read from: {record.path}",
        "",
    ]

    if used_default:
        lines += [
            f"  No level given, so {DEFAULT_SIGMA} sigma was used. That is a generic",
            "  starting point, not a recommendation for this map.",
        ]
        accession = record.evidence.emdb_accession
        if accession:
            lines += [
                "",
                f"  {accession} has an author-recommended contour published by the",
                "  depositor. It is an ABSOLUTE value, so pass it as:",
                f"      density_view({map_obj!r}, {selection!r}, level=<value>, units='absolute')",
                "  Retrieve it from:",
                f"      https://www.ebi.ac.uk/emdb/api/entry/{accession}"
                "  ->  map.contour_list.contour[].level",
            ]
        lines.append("")
    elif units == "absolute":
        lines += [
            "  Level was given in absolute units and converted to sigma, which is",
            "  what PyMOL contours in. Passing an absolute EMDB contour straight",
            "  to PyMOL would contour near zero and show mostly noise.",
            "",
        ]

    lines.append(provenance_banner(map_obj))

    warnings = header.warnings()
    if warnings:
        lines += ["", f"  Geometry warnings ({len(warnings)})"]
        lines += [f"    ! {w}" for w in warnings]

    return "\n".join(lines), scene
