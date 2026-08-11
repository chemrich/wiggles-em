"""Lower a :class:`~wiggles_em.scene.Scene` onto PyMOL.

This is where the σ conversion lives, and it is the reason the seam moved here.
PyMOL normalises MRC/CCP4 maps on load (``normalize_ccp4_maps``, on by
default), so an ``isomesh`` level is in **sigma** — while EMDB's published
author contour is an **absolute** map value. EMD-30913 publishes ``0.05``,
which is 3.16 σ for that map; used directly it contours noise.

Views state levels in the unit the data is in and this module converts, against
**that map's own header**. Two maps loaded in one session are normalised
independently, so a level converted against the wrong header is wrong by
whatever the two maps' statistics differ by — which is a plausible number, not
an obvious one.

Everything else here is naming. ``spectrum``, ``cartoon putty``, ``load_cgo``
and ``mset`` are PyMOL's words for things the scene describes without them.
"""

from __future__ import annotations

import json
from typing import cast

from wiggles_em.atoms import fetch_atoms
from wiggles_em.bfactors import has_stash, preservation_note, stash_bfactors
from wiggles_em.density import to_sigma
from wiggles_em.maps import loaded_map
from wiggles_em.port import PortError, PymolPort, call
from wiggles_em.scene import (
    Arrows,
    ColorByScalar,
    ColorFlat,
    ColorSurfaceByMap,
    Colour,
    Delete,
    Frames,
    Granularity,
    Hide,
    Isosurface,
    Label,
    Legend,
    Opacity,
    Refused,
    Rep,
    ScalarField,
    Scatter,
    Scene,
    SceneOp,
    Sel,
    Show,
    SizeByScalar,
    Unit,
)

#: CGO opcodes. PyMOL's own constants, inlined so this module imports nothing
#: from PyMOL — the package must stay installable without it.
_CGO_CYLINDER = 9.0
_CGO_CONE = 27.0

#: Where per-atom scalars are parked before an ``alter`` reads them back.
#: A dict lookup in one ``alter`` beats one ``alter`` per residue: the round
#: trips stop scaling with the size of the structure.
_STORED = "stored.wiggles_scalar"

#: Representations PyMOL draws. ``MESH`` is not one of them — a mesh is a
#: property of an isosurface object, not of a selection.
_REPS = {
    Rep.CARTOON: "cartoon",
    Rep.STICKS: "sticks",
    Rep.SPHERES: "spheres",
    Rep.SURFACE: "surface",
    Rep.LINES: "lines",
    Rep.EVERYTHING: "everything",
}


def render_selection(sel: Sel) -> str:
    """Lower a :class:`Sel` into a PyMOL selection expression."""
    if sel.kind == "obj":
        return f"({sel.value})"
    if sel.kind == "all":
        return "all"
    if sel.kind == "prop":
        return f"{sel.key} {sel.value}"
    if sel.kind == "lt":
        return f"{sel.key}<{sel.value}"
    if sel.kind == "residues":
        residues = cast("tuple[tuple[str, str], ...]", sel.value)
        if not residues:
            return "none"
        return "(" + " or ".join(f"(chain {c} and resi {r})" for c, r in residues) + ")"
    if sel.kind == "raw":
        if sel.key != "pymol":
            raise Refused(
                f"selection text is in the {sel.key!r} dialect and this backend "
                f"speaks PyMOL. Passing it through might parse into something "
                f"else entirely, which would select the wrong atoms silently."
            )
        return f"({sel.value})"
    if sel.kind == "and":
        return "(" + " and ".join(render_selection(p) for p in sel.parts) + ")"
    if sel.kind == "or":
        return "(" + " or ".join(render_selection(p) for p in sel.parts) + ")"
    if sel.kind == "not":
        return f"(not {render_selection(sel.parts[0])})"
    raise Refused(f"unknown selection kind {sel.kind!r}")


class PymolBackend:
    """Draws a Scene through a :class:`~wiggles_em.port.PymolPort`."""

    def __init__(self, port: PymolPort, preserve_bfactors: bool = True) -> None:
        self.port = port
        self.preserve_bfactors = preserve_bfactors
        #: Absolute levels this backend converted, by surface name. The report
        #: layer reads these back so it can state both units without
        #: recomputing a conversion that might disagree.
        self.converted: dict[str, tuple[float, float]] = {}
        #: Caveats that are true of PyMOL and of no other viewer. The host
        #: appends these to the view's report.
        #:
        #: A view cannot write them: "your B-factors were overwritten, call
        #: restore_bfactors" is a fact about a viewer with one copy of each
        #: object, and protean — which builds a display copy and re-sends it —
        #: has nothing to restore. Putting the note where the behaviour is
        #: keeps the two from drifting apart.
        self.notes: list[str] = []

    # -- entry point -------------------------------------------------------

    def render(self, scene: Scene) -> None:
        """Draw every op, in order. Raises on the first one it cannot honour."""
        for op in scene:
            self.render_op(op)

    def render_op(self, op: SceneOp) -> None:
        handler = getattr(self, f"_{type(op).__name__.lower()}", None)
        if handler is None:
            raise Refused(f"{type(op).__name__} has no PyMOL lowering")
        handler(op)

    # -- scalars -----------------------------------------------------------

    def _push(self, field: ScalarField) -> str:
        """Park a scalar field in PyMOL's ``stored`` namespace.

        Returns the ``alter`` expression that reads it back. Tuple keys are not
        JSON-able, so they are joined on ``|`` and the same join is rebuilt on
        the PyMOL side — the identical trick ``restore_bfactors`` uses, and
        deliberately the same so one atom has one key everywhere.
        """
        flat = {"|".join(str(part) for part in key): value for key, value in zip(field.keys, field.values, strict=True)}
        # do(), not a structured call: this assigns into a namespace rather
        # than invoking a cmd function, which is one of the few things that has
        # no structured equivalent.
        self.port.do(f"{_STORED} = {json.dumps(flat)}")
        if field.granularity is Granularity.ATOM:
            key_expr = "'|'.join((chain, resi, name, alt))"
        else:
            key_expr = "'|'.join((chain, resi))"
        return f"b={_STORED}.get({key_expr}, b)"

    def _stash(self, sel: Sel) -> None:
        """Save B-factors before overwriting them, once per object.

        PyMOL has one copy of an object, so routing a scalar through ``b``
        destroys the crystallographic values. The originals are read back
        rather than taken from the view, because the view's atoms may be a
        subset of what the ``alter`` is about to touch, and a partial stash
        restores a column that is right in places.
        """
        if not self.preserve_bfactors:
            self.notes.append("  WARNING: B-factors overwritten and not preserved.")
            return
        obj = next((s.value for s in sel.walk() if s.kind == "obj"), None)
        if obj is None or has_stash(str(obj)):
            return
        atoms = fetch_atoms(self.port, str(obj))
        stashed = stash_bfactors(str(obj), atoms)
        self.notes.append(preservation_note(str(obj), stashed))

    def _colorbyscalar(self, op: ColorByScalar) -> None:
        sel = render_selection(op.sel)
        self._stash(op.sel)
        call(self.port, "alter", sel, self._push(op.field))
        lo, hi = op.domain
        call(self.port, "spectrum", "b", op.palette, sel, minimum=lo, maximum=hi)

    def _sizebyscalar(self, op: SizeByScalar) -> None:
        sel = render_selection(op.sel)
        self._stash(op.sel)
        call(self.port, "alter", sel, self._push(op.field))
        call(self.port, "show", "cartoon", sel)
        call(self.port, "set", "cartoon_putty_scale_min", op.scale_min)
        call(self.port, "set", "cartoon_putty_scale_max", op.scale_max)
        call(self.port, "cartoon", "putty", sel)

    # -- colour, visibility, labels ---------------------------------------

    def _colour_name(self, colour: Colour) -> str:
        """A PyMOL colour name, defining one for an RGB triple if needed."""
        if isinstance(colour, str):
            return colour
        r, g, b = colour
        name = f"wgf_{int(r * 255):02x}{int(g * 255):02x}{int(b * 255):02x}"
        call(self.port, "set_color", name, [float(r), float(g), float(b)])
        return name

    def _colorflat(self, op: ColorFlat) -> None:
        call(self.port, "color", self._colour_name(op.colour), render_selection(op.sel))

    def _show(self, op: Show) -> None:
        call(self.port, "show", self._rep(op.rep), render_selection(op.sel))

    def _hide(self, op: Hide) -> None:
        call(self.port, "hide", self._rep(op.rep), render_selection(op.sel))

    def _rep(self, rep: Rep) -> str:
        try:
            return _REPS[rep]
        except KeyError:
            raise Refused(
                f"PyMOL has no {rep.value!r} representation for a selection"
            ) from None

    def _label(self, op: Label) -> None:
        # PyMOL evaluates a Python expression against each atom, so field names
        # go in verbatim and the text becomes a format string. A backend that
        # cannot evaluate an expression per atom has to interpolate the values
        # itself — which is exactly why Label carries fields rather than text
        # with the expression already baked into it.
        expr = f'"{op.text}" % ({", ".join(op.fields)},)' if op.fields else f'"{op.text}"'
        call(self.port, "label", render_selection(op.sel), expr)

    def _opacity(self, op: Opacity) -> None:
        sel = render_selection(op.sel)
        transparency = round(1.0 - op.value, 3)
        call(self.port, "set", "cartoon_transparency", transparency, sel)
        call(self.port, "set", "transparency", transparency, sel)

    # -- volumes -----------------------------------------------------------

    def _sigma_for(self, volume: str, level: float, unit: Unit) -> float:
        """The level PyMOL wants, in σ, converted against *this* map's header."""
        if unit is Unit.SIGMA:
            return level
        entry = loaded_map(volume)
        if entry is None:
            raise PortError(
                f"an absolute level was given for {volume!r}, but that volume "
                f"was not loaded through load_map, so its header is unknown "
                f"and the conversion to sigma cannot be made. PyMOL contours "
                f"in sigma; passing an absolute value through would contour at "
                f"the wrong level without failing."
            )
        sigma = to_sigma(entry.header, level)
        self.converted[volume] = (level, sigma)
        return sigma

    def _isosurface(self, op: Isosurface) -> None:
        sigma = self._sigma_for(op.volume, op.level, op.unit)
        action = "isomesh" if op.style is Rep.MESH else "isosurface"
        if op.carve_around is not None:
            if op.carve_radius is None:
                raise Refused("carve_around given without a carve_radius")
            call(
                self.port,
                action,
                op.name,
                op.volume,
                sigma,
                render_selection(op.carve_around),
                carve=op.carve_radius,
            )
        else:
            call(self.port, action, op.name, op.volume, sigma)

    def _colorsurfacebymap(self, op: ColorSurfaceByMap) -> None:
        entry = loaded_map(op.volume)
        if entry is None:
            raise PortError(
                f"{op.volume!r} was not loaded through load_map, so its header "
                f"is unknown and its breakpoints cannot be converted."
            )
        # Breakpoints are in the second volume's own units and convert against
        # ITS header — a different sigma scale from the density map's contour
        # level. This is the half of the sigma trap that is easiest to miss.
        sigmas = [to_sigma(entry.header, point) for point in op.breakpoints]
        ramp = f"{op.surface}_ramp"
        colours = [self._colour_name(c) for c in op.palette]
        call(self.port, "ramp_new", ramp, op.volume, sigmas, colours)
        call(self.port, "set", "surface_color", ramp, op.surface)

    # -- geometry, sequences, lifecycle -----------------------------------

    def _arrows(self, op: Arrows) -> None:
        buffer: list[float] = []
        for start, end, colour in op.segments:
            r, g, b = colour if not isinstance(colour, str) else (1.0, 1.0, 1.0)
            shaft_end = tuple(s + (e - s) * 0.75 for s, e in zip(start, end, strict=True))
            buffer += [
                _CGO_CYLINDER, *start, *shaft_end, op.radius, r, g, b, r, g, b,
                _CGO_CONE, *shaft_end, *end, op.radius * 2.2, 0.0, r, g, b, r, g, b, 1.0, 1.0,
            ]  # fmt: skip
        call(self.port, "load_cgo", buffer, "wgf_arrows")

    def _frames(self, op: Frames) -> None:
        if not op.build_timeline:
            return
        call(self.port, "mset", f"1 x{len(op.names)}")
        prefix = op.names[0].rsplit("_", 1)[0] if op.names else ""
        for index, name in enumerate(op.names, start=1):
            # mdo attaches a command line to a frame; there is no cmd
            # equivalent that takes the command as data.
            call(self.port, "mdo", index, f"disable {prefix}_*; enable {name}")

    def _delete(self, op: Delete) -> None:
        for name in op.names:
            call(self.port, "delete", name)

    def _legend(self, op: Legend) -> None:
        """Draws nothing. Legends are report text; they are in the scene so
        invariants can be checked on the value rather than on a string."""

    def _scatter(self, op: Scatter) -> None:
        raise Refused(
            "Scatter is defined so invariant I2 can name what it forbids, and "
            "no view emits it. Motion is recoverable from a heterogeneity "
            "method and populations are not, so a latent scatter renders a "
            "claim the data does not support."
        )


def draw(port: PymolPort, scene: Scene) -> PymolBackend:
    """Render ``scene`` through ``port``. Returns the backend, for its record."""
    backend = PymolBackend(port)
    backend.render(scene)
    return backend


__all__ = ["PymolBackend", "draw", "render_selection"]
