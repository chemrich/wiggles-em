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
from wiggles_em.bfactors import (
    bfactors_destroyed,
    destroyed_note,
    has_stash,
    mark_bfactors_destroyed,
    preservation_note,
    stash_bfactors,
    stashed_count,
)
from wiggles_em.density import to_absolute, to_sigma
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
    Morph,
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


def _no_quote(value: str) -> str:
    """Return ``value``, or raise if it carries a double quote."""
    if '"' in value:
        raise ValueError(
            f"identifier {value!r} contains a double quote, which cannot appear "
            f"in a chain, residue or object identifier — the file is probably "
            f"corrupt, and guessing at it is how the blank-chain bug behaved"
        )
    return value


def quote(value: str) -> str:
    """Quote an identifier read off an atom, for use in a PyMOL selection.

    PyMOL parses bare identifiers, so interpolating a raw ``chain`` or ``resi``
    is unsafe for the two shapes that occur constantly in real depositions:

    * A **blank chain** — every atom in a file with no chain ID — leaves
      ``chain  and resi 2``, where ``and`` is consumed as the chain name. The
      selection stops being scoped to the object and matches the whole session.
      Quoting fixes this one: ``chain ""`` stays scoped.
    * A **negative residue number**, an expression-tag remnant in most NMR and
      EM entries, is read as a *range*: ``resi -3`` means 1-3, so one residue's
      value gets written across three.

    **Quoting does not fix the second one**, though this package and MCPymol
    both claimed it did. ``resi "-3"`` is still the range — checked against
    PyMOL 3.1.0 on an object holding residues -3, 1 and 2, where it matched all
    six atoms rather than residue -3's two. The escape the grammar honours is a
    **backslash**, and it composes both with the quoting (``resi "\\-3"``
    matches exactly residue -3) and with ``+`` lists (``resi "\\-3"+"1"+"2"``
    matches three residues, not the range).

    So a value is quoted *and* a leading minus escaped. The escape is a no-op
    for every other shape — plain numbers, zero, insertion codes like ``52A``
    and blank, lowercase or numeric chains were all checked against a live
    session, and a chain literally named ``-`` matches under both forms.

    Raises:
        ValueError: ``value`` contains a double quote, which would close the
            quoting and hand the rest of the string to the parser. Nothing in
            the PDB or mmCIF grammar puts one in a chain or residue
            identifier, so this is a corrupt file rather than a case to
            support — and guessing at it is how the blank-chain bug behaved.
    """
    value = _no_quote(value)
    # A leading '-' is the range operator to PyMOL's parser even inside quotes.
    escaped = f"\\{value}" if value.startswith("-") else value
    return f'"{escaped}"'


def render_selection(sel: Sel) -> str:
    """Lower a :class:`Sel` into a PyMOL selection expression."""
    if sel.kind == "obj":
        # Parenthesised rather than quoted: PyMOL takes an object name bare.
        # Still checked for a quote, because an object name reaches here from
        # a tool argument — `m" or all` would parenthesise into
        # `(m" or all)` and select the session.
        return f"({_no_quote(str(sel.value))})"
    if sel.kind == "all":
        return "all"
    if sel.kind == "prop":
        return f"{sel.key} {quote(str(sel.value))}"
    if sel.kind == "lt":
        return f"{sel.key}<{sel.value}"
    if sel.kind == "residues":
        residues = cast("tuple[tuple[str, str], ...]", sel.value)
        if not residues:
            return "none"
        # One term per *chain*, not per residue: 1123 scored residues in one
        # chain become one term rather than 1123 parenthesised `or`s that
        # PyMOL evaluates individually across the whole object.
        #
        # Numbers are never collapsed into ranges — turning 5 and 7 into 5-7
        # would silently add residue 6. Every value goes through `quote`, which
        # also escapes a leading minus, and quoted values compose in a `+` list
        # (`resi "\-3"+"1"+"2"` matches three residues, checked against PyMOL
        # 3.1.0). So there is one path rather than a plain-digit path and an
        # everything-else path that could drift apart.
        by_chain: dict[str, list[str]] = {}
        for chain, resi in residues:
            by_chain.setdefault(chain, []).append(resi)

        terms = []
        for chain, numbers in by_chain.items():
            body = f"resi {'+'.join(quote(n) for n in numbers)}"
            terms.append(f"(chain {quote(chain)} and {body})")
        return terms[0] if len(terms) == 1 else "(" + " or ".join(terms) + ")"
    if sel.kind == "first":
        return f"(first ({render_selection(sel.parts[0])}))"
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


def normalisation_state(port: PymolPort) -> bool | None:
    """Is ``normalize_ccp4_maps`` on? None when PyMOL will not say.

    Lives here rather than with the views because it is a question only PyMOL
    can be asked — it is PyMOL's per-map normalisation that puts a resolution
    ramp's breakpoints in sigma in the first place. ``local_resolution_view``
    takes the answer as an argument, so a host that does not normalise passes
    ``False`` and the arithmetic comes out right for it too.

    Read now, not at load time, which is the honest limitation: a session that
    had it off when the map was loaded and on now would report the wrong thing.

    **What a real PyMOL returns, checked rather than assumed:** ``'1'`` — the
    *string* ``'1'``, not ``'on'`` and not the integer ``1`` (open-source PyMOL,
    2026-08-09, via ``tools/livefire.py``). The parse stays tolerant of the
    other spellings because that answer is one build's, and an older plugin may
    not expose ``get`` at all, which is unknown rather than an error.
    """
    try:
        raw = port.query("get", "normalize_ccp4_maps")
    except PortError:
        return None
    text = str(raw).strip().lower()
    if text in ("on", "1", "1.0", "true", "yes"):
        return True
    if text in ("off", "0", "0.0", "false", "no"):
        return False
    return None


class PymolBackend:
    """Draws a Scene through a :class:`~wiggles_em.port.PymolPort`."""

    def __init__(
        self,
        port: PymolPort,
        *,
        preserve_bfactors: bool = True,
        normalised: bool | None,
    ) -> None:
        self.port = port
        self.preserve_bfactors = preserve_bfactors
        #: Whether PyMOL normalised the volumes on load, as the host read it.
        #:
        #: Taken as an argument rather than queried here, because the *view*
        #: needs the same answer for its report. When each read the session
        #: separately they could disagree, and the failure was a colour key
        #: describing units the surface was not drawn in. One read, one answer:
        #: the host calls :func:`normalisation_state` once and passes it here
        #: and to the view.
        #:
        #: **No default.** It had one — ``None`` — and the public :func:`draw`
        #: helper never passed it, so a host that had correctly read
        #: ``normalize_ccp4_maps off`` and told the view so got a backend that
        #: silently assumed the opposite: Ångström breakpoints converted to
        #: sigma against a volume still holding Ångströms, a surface flat in one
        #: extreme colour, under a report stating they were sent unconverted.
        #: A default is what lets two parts of the system hold different
        #: beliefs, so ``None`` must now be chosen rather than fallen into.
        #:
        #: ``None`` still means "the session would not say", which
        #: :func:`normalisation_state` returns when ``get`` is unavailable. It
        #: is treated as normalised, because that is PyMOL's own default — but
        #: it is now an assertion the caller makes, not one made for them.
        self.normalised = normalised
        #: Levels this backend converted, keyed by **surface name** —
        #: ``{surface: (as given, as sent)}``. Keyed by surface rather than by
        #: volume because one volume routinely carries several surfaces at
        #: different levels, and a volume key silently keeps only the last.
        #:
        #: Written for a report layer that states both units without recomputing
        #: a conversion that might disagree. **Nothing in this package reads it
        #: yet**, so the contract is stated here rather than demonstrated by a
        #: caller — which is why it drifted from the code once already.
        self.converted: dict[str, tuple[float, float]] = {}
        #: Objects this render has already spoken about, so a view emitting
        #: two scalar ops against one object (ColorByScalar + SizeByScalar,
        #: which is every putty view) says it once rather than twice.
        self._noted: set[str] = set()
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
            key_expr = "'|'.join((model, str(rank)))"
        else:
            key_expr = "'|'.join((chain, resi))"
        return f"b={_STORED}.get({key_expr}, b)"

    def _note_once(self, obj: str, note: str) -> None:
        """Append a B-factor note for ``obj``, at most once per render."""
        if obj in self._noted:
            return
        self._noted.add(obj)
        self.notes.append(note)

    def _stash(self, sel: Sel) -> None:
        """Save B-factors before overwriting them, once per object.

        PyMOL has one copy of an object, so routing a scalar through ``b``
        destroys the crystallographic values. The originals are read back
        rather than taken from the view, because the view's atoms may be a
        subset of what the ``alter`` is about to touch, and a partial stash
        restores a column that is right in places.
        """
        obj = next((s.value for s in sel.walk() if s.kind == "obj"), None)
        if obj is None:
            return
        # An existing stash outranks everything below, including
        # preserve_bfactors=False. The originals are already held and
        # `restore_bfactors` still puts them back, so there is nothing to record
        # as lost and nothing to warn about — warning here told the user their
        # crystallographic values were gone when one call would return them.
        #
        # Same ordering `stash_bfactors` gets right; this is the second place
        # that had to, and it was missed when the first was fixed.
        if has_stash(str(obj)):
            self._note_once(str(obj), preservation_note(str(obj), stashed_count(str(obj))))
            return
        if not self.preserve_bfactors:
            # Record the loss, or the next view stashes this view's scalar as
            # though it were the user's data and restore_bfactors writes it
            # back reporting success.
            mark_bfactors_destroyed(str(obj))
            self._note_once(str(obj), "  WARNING: B-factors overwritten and not preserved.")
            return
        if bfactors_destroyed(str(obj)):
            self._note_once(str(obj), destroyed_note(str(obj)))
            return
        atoms = fetch_atoms(self.port, str(obj))
        stashed = stash_bfactors(str(obj), atoms)
        self._note_once(str(obj), preservation_note(str(obj), stashed))

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

    def _level_for(
        self,
        surface: str,
        volume: str,
        level: float,
        unit: Unit,
        equivalent: float | None = None,
    ) -> float:
        """The number PyMOL wants, converted against *this* map's header.

        Which unit that is depends on the session, not on the op. PyMOL
        contours in sigma **because** it normalised the map on load; with
        ``normalize_ccp4_maps`` off the values are stored as written and an
        isosurface level is read as an absolute map value. Sending a sigma
        number to an unnormalised map contours far above its dmax and yields an
        empty surface, while the report cheerfully states the absolute level it
        thought it had asked for.
        """
        wanted = Unit.SIGMA if self.normalised is not False else Unit.ABSOLUTE
        if unit is wanted:
            return level

        # The view may already hold this contour in the other unit. An ensemble
        # frame has no load_map record, so the header lookup below cannot
        # succeed for it and the remedy the error suggests is impossible —
        # nothing was drawn at all. Where the view could compute the equivalent
        # it carries it, and it is the only copy in existence.
        if equivalent is not None:
            # Recorded like any other conversion. Returning early skipped it,
            # leaving nothing recorded for exactly the surfaces this field was
            # added to make renderable.
            self.converted[surface] = (level, equivalent)
            return equivalent

        entry = loaded_map(volume, self.port)
        if entry is None:
            raise PortError(
                f"a level in {unit.value} was given for {volume!r}, but this "
                f"session needs it in {wanted.value} and that volume was not "
                f"loaded through load_map, so its header is unknown and the "
                f"conversion cannot be made. Passing the number through would "
                f"contour at the wrong level without failing."
            )
        converted = (
            to_sigma(entry.header, level)
            if wanted is Unit.SIGMA
            else to_absolute(entry.header, level)
        )
        self.converted[surface] = (level, converted)
        return converted

    def _isosurface(self, op: Isosurface) -> None:
        level = self._level_for(op.name, op.volume, op.level, op.unit, op.equivalent)
        action = "isomesh" if op.style is Rep.MESH else "isosurface"
        if op.carve_around is not None:
            if op.carve_radius is None:
                raise Refused("carve_around given without a carve_radius")
            call(
                self.port,
                action,
                op.name,
                op.volume,
                level,
                render_selection(op.carve_around),
                carve=op.carve_radius,
            )
        else:
            call(self.port, action, op.name, op.volume, level)

    def _colorsurfacebymap(self, op: ColorSurfaceByMap) -> None:
        # With the port: a record whose object has left the session is evicted
        # rather than used, because converting against a volume that is no
        # longer loaded gives a wrong number and no error.
        entry = loaded_map(op.volume, self.port)
        if entry is None:
            raise PortError(
                f"{op.volume!r} was not loaded through load_map, or the object "
                f"has since left the session, so its header is unknown and its "
                f"breakpoints cannot be converted."
            )
        # Breakpoints are in the second volume's own units (Angstrom) and
        # convert against ITS header — a different sigma scale from the density
        # map's contour level. This is the half of the sigma trap that is
        # easiest to miss.
        #
        # Unless PyMOL was told not to normalise, in which case the stored
        # values are still resolutions and converting would be the error. An
        # unanswerable question is the unknown case, and the PyMOL default is
        # on, so unknown converts.
        if self.normalised is False:
            sigmas = list(op.breakpoints)
        else:
            sigmas = [to_sigma(entry.header, point) for point in op.breakpoints]
        ramp = f"{op.surface}_ramp"
        colours = [self._colour_name(c) for c in op.palette]
        call(self.port, "ramp_new", ramp, op.volume, sigmas, colours)
        call(self.port, "set", "surface_color", ramp, op.surface)
        self.notes.append(
            f"  Coloured through PyMOL ramp `{ramp}`. The ramp must outlive the "
            f"surface — deleting it un-colours it, because the colour is not "
            f"baked into the mesh."
        )

    # -- geometry, sequences, lifecycle -----------------------------------

    def _arrows(self, op: Arrows) -> None:
        buffer: list[float] = []
        for arrow in op.segments:
            colour = arrow.colour
            r, g, b = colour if not isinstance(colour, str) else (1.0, 1.0, 1.0)
            shaft_end = tuple(
                s + (e - s) * 0.75 for s, e in zip(arrow.start, arrow.end, strict=True)
            )
            buffer += [
                _CGO_CYLINDER, *arrow.start, *shaft_end, arrow.radius, r, g, b, r, g, b,
                _CGO_CONE, *shaft_end, *arrow.end, arrow.radius * 2.2, 0.0,
                r, g, b, r, g, b, 1.0, 1.0,
            ]  # fmt: skip
        if not buffer:
            return
        # Delete first: load_cgo onto an existing name appends rather than
        # replaces, so a second call would leave the previous arrows behind and
        # the picture would show two displacements at once.
        call(self.port, "delete", op.name)
        call(self.port, "load_cgo", buffer, op.name)

    def _frames(self, op: Frames) -> None:
        if not op.build_timeline:
            return
        # The timeline is numbered by the frames the surfaces were MADE from,
        # not by their position in the surviving list. `enumerate(..., start=1)`
        # renumbered over the survivors, so with frame 3 skipped, `frame 3`
        # showed frame 4's density and `frame 5` showed nothing — while the
        # report promised "a surface's number is always the frame it was made
        # from" and told the reader to type `frame N` to step.
        by_number = dict(zip(op.numbers, op.names, strict=True))
        if not by_number:
            # `latent_traverse_view` guards on len(surfaces) > 1 and never emits
            # this, but Frames is public scene API. `max()` over nothing raised
            # ValueError, which is an internal error escaping as the viewer's
            # answer; an empty timeline is simply nothing to build.
            return
        span = max(by_number)
        all_disabled = "; ".join(f"disable {n}" for n in op.names)
        call(self.port, "mset", f"1 x{span}")
        for number in range(1, span + 1):
            name = by_number.get(number)
            # Every sibling named explicitly. The previous version recovered a
            # prefix with rsplit and disabled `prefix_*`, which also switched
            # off an unrelated `v_model` in the session on every frame step —
            # the wildcard the Delete op's contract exists to forbid. Frames
            # already carries the exact names, so no string surgery is needed.
            #
            # Quadratic in the number of frames, in string length. A traversal
            # is normally tens of frames and these are the cheapest commands
            # PyMOL has — but the cost is O(max frame number), not O(surfaces),
            # and those diverge in exactly the gap case this numbering exists
            # for: 10 usable frames numbered 491-500 build a 500-frame timeline,
            # 490 of them "disable everything". The round trips are inherent
            # (PyMOL runs an mdo when its frame is *entered*, so jumping into
            # the middle of a gap run must still clear the view), but the
            # identical gap command is built once rather than 490 times.
            others = all_disabled if name is None else "; ".join(
                f"disable {other}" for other in op.names if other != name
            )
            # A skipped frame disables everything rather than being left out of
            # the timeline. Omitting it would leave the previous frame's surface
            # on screen while the report says that frame was never contoured —
            # showing one frame's density under another's number.
            #
            # Both halves can be empty: a lone surface has no siblings to
            # disable, and a skipped frame enables nothing. Joining them
            # unconditionally left a leading `; ` — the guard that was in the
            # code this loop replaced.
            command = "; ".join(part for part in (others, f"enable {name}" if name else "") if part)
            # mdo attaches a command line to a frame; there is no cmd
            # equivalent that takes the command as data.
            call(self.port, "mdo", number, command)
        # mdo commands run when a frame is *entered*. Without this the timeline
        # is built and never executed, so every isosurface stays enabled and
        # the traversal renders as one superimposed blob.
        #
        # The FIRST REAL frame, not frame 1. Frame 1 is a gap whenever frame 1's
        # header carries no usable rms, and a gap frame disables everything — so
        # hard-coding 1 ended a perfectly successful traversal on a blank
        # viewport, under a report listing the surfaces it had just built. The
        # old positional numbering always enabled the first surface; carrying
        # real frame numbers took that away silently.
        call(self.port, "frame", min(by_number))

    def _morph(self, op: Morph) -> None:
        try:
            call(self.port, "morph", op.name, op.obj, refinement=0, steps=op.steps)
        except PortError as exc:
            if "incentive-only" not in str(exc).lower():
                raise
            # cmd.morph is Incentive-only. Confirmed against open-source PyMOL
            # on 2026-08-08, which is what most users run — so for most users
            # the interpolation half of morph_states does not exist here. The
            # judgement half is the topology check, it already ran in the view,
            # and it does not depend on a licence. Reporting beats raising.
            #
            # None of this is true of protean, which interpolates natively, so
            # it belongs in this backend rather than in the view's report.
            self.notes += [
                "  Morph NOT created: cmd.morph is Incentive-only and this is",
                "  open-source PyMOL.",
                "",
                "  The topology check is the part of this tool that carries a",
                "  judgement, and it passed — these states can be interpolated",
                "  meaningfully. To see the motion without Incentive PyMOL, step",
                "  the states directly (they are observed, unlike interpolated",
                "  frames):",
                "      set all_states, on      # all of them at once, or",
                "      mset 1 x<n_states>      # play them as frames",
                "",
                "  ensemble_spread_view shows the same variation as a static image.",
            ]

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


def draw(port: PymolPort, scene: Scene, *, normalised: bool | None) -> PymolBackend:
    """Render ``scene`` through ``port``. Returns the backend, for its record.

    ``normalised`` is required and passed straight through: this helper used to
    construct the backend with the default, so every host that rendered through
    it got a backend assuming a normalised session however carefully it had
    read the real answer and told the view.

    Reading the session here instead was the alternative, and it is the one
    thing this must not do — the view has already been given an answer, and a
    second independent read is exactly how the two came to disagree. The host
    calls :func:`normalisation_state` once and hands the result to both.
    """
    backend = PymolBackend(port, normalised=normalised)
    backend.render(scene)
    return backend


__all__ = ["PymolBackend", "draw", "normalisation_state", "render_selection"]
