"""What a view wants drawn, as a value rather than a sequence of calls.

A view computes and returns. It never calls a viewer. Everything above
:mod:`wiggles_em.backends` is pure: no sockets, no command strings, no
viewer-specific names.

**Why this exists rather than a wider port protocol.** The old ``PymolPort``
abstracted the *transport*, and did that well — protean's viewer bridge is the
same ``{action, args}`` shape over a WebSocket. What does not port is the
vocabulary: ``spectrum b, red_white_blue, obj, minimum=0, maximum=1`` names a
``pymol.cmd`` function and Mol\\* has never heard of it.

Three things fall out of moving the seam up a layer, in the order they matter:

1. **Normalisation becomes a backend's problem.** Every level and domain here
   is stated in the units the *data* is in — an absolute contour level, a
   fixed ``(0.0, 1.0)`` occupancy domain. PyMOL contours in σ and normalises
   each map independently, so its backend converts; Mol\\* carries
   ``Volume.IsoValue.absolute`` natively and passes straight through. This is
   the trap that has cost this project the most time, and a view can no longer
   get it wrong because a view no longer decides.

2. **Invariants become assertions on a value.** "Draws no latent scatter" was
   a substring search over a recorded command log. Here it is
   ``not any(isinstance(op, Scatter) for op in scene.ops)`` — a claim about
   what the view returned, which cannot pass for the wrong reason.

3. **A view is testable with no viewer and no fake at all.** Call it, look at
   what came back.

The cost is that a view cannot read state part-way through drawing. Checked
against all fourteen: the *views* read up front and draw at the end, so nothing
is lost. The **loaders** (``load_map``, ``load_ensemble``) genuinely do need to
interleave — they issue a load and then query back to confirm the object
arrived, which is the MCPymol issue #15 discipline — so they keep taking a port
and are not views. That distinction is real and is the reason this module does
not try to describe loading.
"""

from __future__ import annotations

from collections.abc import Iterable
from dataclasses import dataclass
from enum import Enum

# --------------------------------------------------------------------------
# Selections
# --------------------------------------------------------------------------


class Op(str, Enum):
    """How a compound selection combines its parts."""

    AND = "and"
    OR = "or"
    NOT = "not"


@dataclass(frozen=True)
class Sel:
    """A selection, as a value a backend lowers into its own dialect.

    Scene ops cannot carry PyMOL selection strings: protean parses a subset
    with real gaps, and a string is not something a backend can inspect before
    deciding whether it can honour it.

    Five shapes cover every selection the fourteen views build:

    ``Sel.obj("6xyz")``                 an object or a named handle
    ``Sel.prop("name", "CA")``          a property equals a value
    ``Sel.lt("q", 0.999)``              a numeric comparison
    ``Sel.residues([("A", "12"), …])``  a set of residues, compactly
    ``Sel.raw("byres polymer")``        text the caller supplied; see below

    Combine with ``&``, ``|`` and ``~``.
    """

    kind: str
    key: str = ""
    value: object = None
    parts: tuple[Sel, ...] = ()

    # -- constructors ------------------------------------------------------

    @classmethod
    def obj(cls, name: str) -> Sel:
        """An object or named handle."""
        return cls("obj", value=name)

    @classmethod
    def prop(cls, key: str, value: object) -> Sel:
        """A property equal to a value — ``chain A``, ``name CA``, ``alt B``."""
        return cls("prop", key=key, value=value)

    @classmethod
    def lt(cls, key: str, value: float) -> Sel:
        """A numeric property strictly below ``value`` — ``q<0.999``."""
        return cls("lt", key=key, value=value)

    @classmethod
    def residues(cls, residues: list[tuple[str, str]]) -> Sel:
        """A set of ``(chain, resi)`` pairs.

        A first-class shape rather than a chain of ``|``, because Q-score views
        build these with a thousand-plus entries and a backend wants to emit
        one compact expression rather than a thousand-term disjunction.
        """
        return cls("residues", value=tuple(residues))

    @classmethod
    def raw(cls, text: str, dialect: str = "pymol") -> Sel:
        """Selection text the *caller* supplied, in a named dialect.

        ``composition_view`` takes a selection from the user and scopes it to
        an object. That text cannot be interpreted here, and guessing at it is
        exactly how a selection silently matches the wrong atoms.

        Marking it raw makes it auditable: a backend that does not speak the
        dialect **refuses**, rather than passing a string that might happen to
        parse into something else entirely.
        """
        return cls("raw", key=dialect, value=text)

    @classmethod
    def first(cls, inner: Sel) -> Sel:
        """One atom from ``inner`` — whichever the viewer considers first.

        Labelling wants exactly one atom per part; labelling every atom is
        unreadable. An index-based narrowing cannot do it, because an atom
        index is numbered over the whole object and ANDing ``rank 0`` with a
        per-part selection matches nothing for every part but one. PyMOL and
        protean both have a first-of operator, which is the honest primitive.
        """
        return cls("first", parts=(inner,))

    @classmethod
    def all(cls) -> Sel:
        """Everything. Used by ``Hide`` to clear a slate."""
        return cls("all")

    # -- algebra -----------------------------------------------------------

    def __and__(self, other: Sel) -> Sel:
        return Sel(Op.AND.value, parts=(self, other))

    def __or__(self, other: Sel) -> Sel:
        return Sel(Op.OR.value, parts=(self, other))

    def __invert__(self) -> Sel:
        return Sel(Op.NOT.value, parts=(self,))

    def walk(self):
        """Yield this selection and every one nested inside it."""
        yield self
        for part in self.parts:
            yield from part.walk()

    @property
    def dialects(self) -> set[str]:
        """Every raw dialect anywhere in this selection. Empty if none."""
        return {s.key for s in self.walk() if s.kind == "raw"}


# --------------------------------------------------------------------------
# Vocabulary
# --------------------------------------------------------------------------


class Rep(str, Enum):
    """A representation. What a backend actually calls it is its business."""

    CARTOON = "cartoon"
    STICKS = "sticks"
    SPHERES = "spheres"
    SURFACE = "surface"
    MESH = "mesh"
    LINES = "lines"
    EVERYTHING = "everything"


class Unit(str, Enum):
    """The unit a contour level is stated in.

    Never defaulted anywhere in this package. EMDB publishes author-recommended
    contour levels as ABSOLUTE values while most viewers contour in SIGMA;
    EMD-30913 publishes 0.05, which is 3.16 σ for that map, and used as σ it
    contours noise. A bare number is unusable and the enum makes that structural.
    """

    ABSOLUTE = "absolute"
    SIGMA = "sigma"


class Sense(str, Enum):
    """Which quantity called "occupancy" a view is showing.

    The two are incompatible and a render that conflates them looks entirely
    normal, so every view that touches either declares which one it means.
    """

    #: Per-atom crystallographic ``q``, read from the coordinate file.
    ATOM_OCCUPANCY = "atom-occupancy"
    #: Fraction of imaged *particles* containing a subunit.
    PARTICLE_COMPOSITION = "particle-composition"


Colour = str | tuple[float, float, float]
"""A named colour, or RGB in 0–1. Backends that need a name define one.

**A name is an argument, not a scene value.** Every name below is one of
*PyMOL's*, so a Scene carrying one can only be honoured by a second viewer that
reimplements PyMOL's table — which is a viewer-neutral value naming a viewer.
Views therefore call :func:`resolve_colour` on whatever a caller passed and put
RGB in the Scene. Callers keep naming colours, because a view's signature is
not the seam; the Scene is.
"""


#: PyMOL's RGB for the names this package uses or accepts. Assembled from the
#: two copies that already existed — ``composition._blend``'s table and
#: protean's ``backends/molstar.py::_COLOUR_NAMES`` — which had **disjoint**
#: name sets and agreed on the one name they shared (``skyblue``). That
#: divergence is why this lives in one place now.
#:
#: Deliberately not exhaustive. A name whose value this package cannot state
#: confidently is refused rather than approximated, because a wrong RGB draws a
#: plausible picture and reports success. Callers wanting anything else pass an
#: RGB triple, which needs no table at all.
_PALETTE: dict[str, tuple[float, float, float]] = {
    # Unambiguous in any viewer.
    "red": (1.0, 0.0, 0.0),
    "green": (0.0, 1.0, 0.0),
    "blue": (0.0, 0.0, 1.0),
    "yellow": (1.0, 1.0, 0.0),
    "cyan": (0.0, 1.0, 1.0),
    "magenta": (1.0, 0.0, 1.0),
    "orange": (1.0, 0.5, 0.0),
    "white": (1.0, 1.0, 1.0),
    "black": (0.0, 0.0, 0.0),
    # PyMOL's pastels, used for altloc groups and defaults.
    "skyblue": (0.34, 0.63, 0.83),
    "salmon": (1.0, 0.6, 0.6),
    "palegreen": (0.65, 0.9, 0.65),
    "wheat": (0.99, 0.82, 0.65),
    "lightpink": (1.0, 0.75, 0.87),
    "paleyellow": (1.0, 1.0, 0.5),
    "lightblue": (0.75, 1.0, 1.0),
    "lightorange": (1.0, 0.8, 0.5),
}

#: ``grey``/``gray`` with no suffix. The suffixed family is computed — see
#: :func:`resolve_colour`.
_GREY = (0.5, 0.5, 0.5)


def resolve_colour(colour: Colour) -> tuple[float, float, float]:
    """RGB in 0–1 for ``colour``, which may already be RGB.

    Raises:
        ValueError: for a name with no recorded value. **Refusing beats
            guessing**: an approximated colour renders a plausible picture and
            returns cleanly, which is the failure mode that costs most here.
            The message names the escape hatch, since a caller who wants a
            colour this package does not know can always pass its RGB.
    """
    if not isinstance(colour, str):
        r, g, b = colour
        return (float(r), float(g), float(b))

    name = colour.strip().lower()
    if name in _PALETTE:
        return _PALETTE[name]
    if name in ("grey", "gray"):
        return _GREY
    # PyMOL defines grey00-grey99 (and the `gray` spelling) as a uniform ramp
    # at NN/100, which is why `grey50` and `grey70` are not table entries: the
    # rule is exact, so computing it cannot drift from PyMOL the way a
    # transcribed value can.
    if name[:4] in ("grey", "gray") and name[4:].isdigit() and len(name[4:]) == 2:
        level = int(name[4:]) / 100.0
        return (level, level, level)

    raise ValueError(
        f"{colour!r} is not a colour this package has a value for. Known names: "
        f"{', '.join(sorted(_PALETTE))}, and grey00-grey99. Pass an RGB triple "
        f"in 0-1 for anything else — a Scene carries RGB, so a triple needs no "
        f"table and works in every viewer. Guessing at the name was rejected "
        f"deliberately: a wrong colour draws a picture that looks fine."
    )


def ramp(spec: str | Iterable[Colour]) -> tuple[tuple[float, float, float], ...]:
    """Colour stops for a scalar ramp, low value first.

    Accepts PyMOL's underscore spelling — ``"blue_white_red"`` — because that
    is how these ramps have always been written down here and it stays the
    readable way to say one. It is resolved *at construction*, so what the
    Scene carries is stops: a second viewer needs no palette vocabulary, only
    the ability to interpolate between colours, and the **direction** becomes
    a value a backend can inspect instead of a convention it has to know.

    That direction is the whole point. ``red_white_blue`` and
    ``blue_white_red`` differ only in order, a viewer with one fixed ramp can
    honour exactly one of them, and drawing the other reverses the reading of
    every number on screen while looking entirely normal.

    Raises:
        ValueError: fewer than two stops, or a name that is not resolvable.
    """
    parts: list[Colour] = list(spec.split("_")) if isinstance(spec, str) else list(spec)
    if len(parts) < 2:
        raise ValueError(
            f"a ramp needs at least two stops to interpolate between; got "
            f"{len(parts)} from {spec!r}"
        )
    return tuple(resolve_colour(part) for part in parts)


#: Low value red, high value blue. Occupancy's default: q runs 0-1 and the
#: legend states the direction, because red-means-more is the commoner reading
#: and this is deliberately the other one.
RED_WHITE_BLUE = ramp("red_white_blue")

#: Low value blue, high value red — the everyday "more is hotter" ramp, used
#: for displacement and spread.
BLUE_WHITE_RED = ramp("blue_white_red")

#: Low value red, high green. For scores where low is bad rather than merely
#: small, so the ramp carries a judgement the other two do not.
RED_YELLOW_GREEN = ramp("red_yellow_green")


class Granularity(str, Enum):
    """What a scalar field's keys identify."""

    #: Keys are ``(model, rank)`` — :attr:`wiggles_em.atoms.Atom.key`, the same
    #: key the B-factor stash uses, so an atom is identified identically
    #: wherever it appears. Build them with ``atom.key``, never by hand:
    #: ``(chain, resi, name, alt)`` collides on insertion codes and across two
    #: loaded structures, and a field built that way matches nothing at all,
    #: which draws an ordinary-looking ramp of the *previous* column.
    ATOM = "atom"
    #: Keys are ``(chain, resi)``.
    RESIDUE = "residue"


@dataclass(frozen=True)
class ScalarField:
    """Per-atom or per-residue values, carried with the keys they belong to.

    Keyed rather than positional. A bare parallel array assumes the backend
    reads atoms in the same order the view did, which is true for one transport
    and an unstated assumption for any other — and when it breaks, it does not
    fail, it colours the wrong atoms.

    Granularity is explicit because the two are drawn differently: PyMOL writes
    a per-residue field with one ``alter`` per residue, or with a single
    dictionary lookup for either. Which one a quantity is measured at is a fact
    about the quantity, not an implementation detail.
    """

    keys: tuple[tuple[str, ...], ...]
    values: tuple[float, ...]
    granularity: Granularity

    def __post_init__(self) -> None:
        if len(self.keys) != len(self.values):
            raise ValueError(
                f"{len(self.keys)} keys but {len(self.values)} values — a "
                f"scalar field with a length mismatch would colour by an "
                f"offset, which renders plausibly and is wrong"
            )
        # Every key is a pair — (model, rank) for atoms, (chain, resi) for
        # residues. The contract used to be a 4-tuple (chain, resi, name, alt),
        # and a field still built that way matches nothing on the viewer side
        # while raising nothing: `alter` becomes a no-op and the structure is
        # coloured by whatever the B-factor column already held, over the new
        # quantity's domain, under a legend naming a quantity never drawn.
        # Documenting that was weaker than refusing it, and every other hazard
        # in this dataclass is refused.
        # `isinstance` before `len`, and sorted by repr: the malformed keys this
        # exists to catch are exactly the ones that break a naive check. A bare
        # `rank` has no len() at all, and the old 4-tuple contract's `resi` was
        # an int in some callers and a str in others, so sorting the raw keys
        # raised TypeError comparing the two. A guard that dies on its own
        # subject matter is worse than the silent wrongness it replaced.
        wrong = sorted(
            {repr(k): k for k in self.keys if not isinstance(k, tuple) or len(k) != 2}.values(),
            key=repr,
        )
        if wrong:
            raise ValueError(
                f"scalar field keys must have two parts — (model, rank) per "
                f"atom, (chain, resi) per residue — but got "
                f"{wrong[:3]}{' …' if len(wrong) > 3 else ''}. Build per-atom "
                f"keys with `atom.key`; a key of any other shape matches no "
                f"atom and colours the structure by the previous column instead."
            )
        if len(set(self.keys)) != len(self.keys):
            # A backend lowers this to a lookup table, so a repeated key
            # silently keeps one value and hands it to every atom that shares
            # the key. Nothing about the result looks unusual: the wrong atom
            # is drawn a legitimate colour from the right scale.
            seen: set[tuple[str, ...]] = set()
            # `key=repr` for the same reason as the arity guard above: these are
            # caller-supplied keys reached *because* something is wrong with
            # them, and a per-residue `resi` is an int in some callers and a str
            # in others, so ordering the raw values raises TypeError from inside
            # the error path. Fixing one of these two guards and not the other
            # is how this survived a round.
            clashes = sorted(
                {k for k in self.keys if k in seen or seen.add(k)},  # type: ignore[func-returns-value]
                key=repr,
            )
            raise ValueError(
                f"duplicate keys in a per-{self.granularity.value} scalar field: "
                f"{clashes[:5]}{' …' if len(clashes) > 5 else ''}. Each key must "
                f"identify one thing, or the value drawn for it is whichever "
                f"happened to be last — and that renders as an ordinary colour "
                f"on the right scale, so nothing looks wrong."
            )

    @classmethod
    def per_atom(cls, pairs) -> ScalarField:
        """From ``[(atom.key, value), …]`` — that is, ``[((model, rank), v), …]``.

        Use :attr:`wiggles_em.atoms.Atom.key` rather than assembling the tuple.
        A key of any other shape matches nothing on the viewer side: ``alter``
        becomes a no-op and the structure is coloured by whatever the B-factor
        column already held — over the new quantity's domain, under a legend
        naming a quantity that was never drawn. A wrong-arity key is now
        refused rather than merely documented, but a two-part key of the wrong
        *kind* still cannot be caught here.
        """
        keys, values = zip(*pairs, strict=True) if pairs else ((), ())
        return cls(tuple(keys), tuple(float(v) for v in values), Granularity.ATOM)

    @classmethod
    def per_residue(cls, pairs) -> ScalarField:
        """From ``[((chain, resi), value), …]``."""
        keys, values = zip(*pairs, strict=True) if pairs else ((), ())
        return cls(tuple(keys), tuple(float(v) for v in values), Granularity.RESIDUE)

    @property
    def span(self) -> tuple[float, float]:
        """Observed min and max. Not a domain — see :class:`ColorByScalar`."""
        return (min(self.values), max(self.values)) if self.values else (0.0, 0.0)

    def __len__(self) -> int:
        return len(self.values)


# --------------------------------------------------------------------------
# Ops
# --------------------------------------------------------------------------


@dataclass(frozen=True)
class SceneOp:
    """Base class. Every op is a frozen dataclass and carries no behaviour."""


@dataclass(frozen=True)
class ColorByScalar(SceneOp):
    """Colour a selection by a per-atom scalar over an **explicit** domain.

    ``domain`` is never optional and never inferred. Colouring occupancy over
    its observed range rather than a fixed ``(0.0, 1.0)`` turns a model that is
    0.95–1.0 everywhere into a full rainbow, implying variation that is not
    there. Whether a domain is fixed or data-derived is a property of the
    quantity, so the view states it.

    ``values`` is parallel to ``atoms`` — the atom list the view read, in the
    order it read it. Both viewers ramp per-atom scalars through the B-factor
    column (PyMOL's ``spectrum b``, Mol\\*'s ``uncertainty`` size and colour
    themes read ``B_iso_or_equiv``), so both backends smuggle it the same way.
    They differ in destructiveness, not mechanism, and that is a backend's
    business: PyMOL has one copy and must stash the originals, protean builds a
    display copy and re-sends it.

    ``palette`` is colour **stops**, low value first, not a palette name.
    PyMOL's ``spectrum`` takes ``"blue_white_red"``; naming that in a Scene
    made the ramp a thing a second viewer had to already know, and one with a
    single fixed ramp could not tell that ``red_white_blue`` asked for the
    reverse of it. Stops make the direction inspectable. Build them with
    :func:`ramp`, which still accepts the underscore spelling.
    """

    sel: Sel
    field: ScalarField
    domain: tuple[float, float]
    palette: tuple[Colour, ...] = RED_WHITE_BLUE

    def __post_init__(self) -> None:
        if isinstance(self.palette, str):
            raise ValueError(
                f"palette is colour stops, not a palette name: got "
                f"{self.palette!r}. Use ramp({self.palette!r}), which resolves "
                f"the underscore spelling. Passing the string would iterate it "
                f"character by character and resolve none of them."
            )
        if len(self.palette) < 2:
            raise ValueError(
                f"a ramp needs at least two stops to interpolate between; got {len(self.palette)}"
            )


@dataclass(frozen=True)
class SizeByScalar(SceneOp):
    """Vary thickness by a per-atom scalar — a putty.

    PyMOL lowers this to ``cartoon putty`` plus the putty scale settings. Mol\\*
    lowers it to a cartoon with the ``uncertainty`` size theme, which computes
    ``baseSize + B_iso_or_equiv * bfactorFactor`` over the same display copy.

    A backend that cannot vary thickness honestly must **refuse** rather than
    fall back to colour silently: a view that asked for thickness and got
    colour reads as though the quantity were unavailable.
    """

    sel: Sel
    field: ScalarField
    domain: tuple[float, float]
    scale_min: float = 0.4
    scale_max: float = 3.0


@dataclass(frozen=True)
class ColorFlat(SceneOp):
    """One colour across a selection."""

    sel: Sel
    colour: Colour


@dataclass(frozen=True)
class Show(SceneOp):
    """Show a representation."""

    sel: Sel
    rep: Rep


@dataclass(frozen=True)
class Hide(SceneOp):
    """Hide a representation."""

    sel: Sel
    rep: Rep


@dataclass(frozen=True)
class Label(SceneOp):
    """Label atoms in a selection.

    ``text`` is literal. Where a view wants a per-atom value in the label it
    passes ``fields``, naming atom properties the backend substitutes — PyMOL
    evaluates a Python expression against the atom, Mol\\* has to build the
    string itself, and a raw expression string would only work on one of them.
    """

    sel: Sel
    text: str
    fields: tuple[str, ...] = ()


@dataclass(frozen=True)
class Opacity(SceneOp):
    """Set opacity, 0 (invisible) to 1 (solid).

    Stated as opacity rather than transparency because the two are inverses and
    every viewer picks a different one; the conversion belongs in the backend.
    """

    sel: Sel
    value: float


@dataclass(frozen=True)
class Isosurface(SceneOp):
    """Contour a volume.

    ``level`` is stated with its :class:`Unit` and no default. ``carve_radius``
    limits the surface to within that many Å of ``carve_around``; a backend
    with no carve of its own answers it by cropping the volume before the data
    is ever sent, which is a preprocessing step rather than a drawing one.
    """

    name: str
    volume: str
    level: float
    unit: Unit
    style: Rep = Rep.MESH
    carve_around: Sel | None = None
    carve_radius: float | None = None
    #: The same contour expressed in the *other* unit, when the view could
    #: work it out and the backend cannot.
    #:
    #: A backend converts against the volume's header, which it reads from the
    #: ``load_map`` record. Ensemble frames have no such record — they are
    #: registered by ``load_ensemble`` — so a backend asked for the other unit
    #: had no way to get there and refused, telling the user to load the volume
    #: through ``load_map``, which for an ensemble frame is impossible. Nothing
    #: was drawn at all.
    #:
    #: ``latent_traverse_view`` holds one absolute level and converts it to each
    #: frame's own sigma using headers it already has. Carrying the absolute
    #: value alongside costs nothing and is the only copy that exists.
    equivalent: float | None = None


@dataclass(frozen=True)
class ColorSurfaceByMap(SceneOp):
    """Colour an isosurface by the values of a *second* volume.

    Local resolution: a density surface coloured by a per-voxel resolution
    field. ``breakpoints`` are in the second volume's own units (Å), ramped
    through ``palette``.

    A backend must refuse when the two volumes do not share a voxel grid.
    Sampling one at the other's coordinates does not render visibly broken — it
    renders smooth and plausible, which is worse.
    """

    surface: str
    volume: str
    breakpoints: tuple[float, ...]
    palette: tuple[Colour, ...]


@dataclass(frozen=True)
class Arrow:
    """One displacement arrow, in model coordinates.

    Radius is per-arrow because a long displacement drawn at the same width as
    a short one reads as equally confident; the views scale it with length.
    """

    start: tuple[float, float, float]
    end: tuple[float, float, float]
    colour: Colour
    radius: float = 0.15


@dataclass(frozen=True)
class Arrows(SceneOp):
    """Displacement arrows, as explicit geometry.

    Each :class:`Arrow` is start, end, colour and radius in model coordinates.
    PyMOL builds a CGO; Mol\\* builds a custom shape. Neither has a semantic
    notion of "an arrow between two atoms", so the view does the arithmetic and
    hands over geometry — which also keeps the arrowhead proportions identical
    between viewers instead of leaving each to invent them.

    ``name`` is the object the arrows land in, so a backend can replace the
    previous set rather than accumulating them.
    """

    segments: tuple[Arrow, ...]
    name: str = "wgf_arrows"


@dataclass(frozen=True)
class Frames(SceneOp):
    """A steppable sequence, one object visible at a time.

    A latent traversal: each frame is an isosurface already emitted by an
    :class:`Isosurface` op. ``build_timeline`` asks for viewer playback as
    well as the objects.

    Every frame holds the **same absolute level**. Contouring each frame at a
    fixed σ contours it against its own normalisation, which flattens away the
    density change the traversal exists to show.

    ``numbers`` is the frame each surface was made from, and it is **not**
    derivable from position: a frame whose header carries no usable RMS is
    skipped, so the surfaces run 1, 2, 4, 5 with a gap where 3 would be. The
    view knows those numbers and used to drop them, leaving the backend to
    number the timeline by position — which silently repointed ``frame 4`` at
    frame 5's density under a report promising the opposite.

    No default, deliberately. A default would let the scene and the backend
    hold different beliefs about the numbering, which is the bug itself.
    """

    names: tuple[str, ...]
    numbers: tuple[int, ...]
    build_timeline: bool = False

    def __post_init__(self) -> None:
        if len(self.names) != len(self.numbers):
            raise ValueError(
                f"{len(self.names)} frame names but {len(self.numbers)} numbers — "
                f"a timeline built from these would step to the wrong density"
            )
        # Integers first: everything below compares them — `sorted()` here and
        # `n < 1` further down — and a string sails through the length check
        # only to raise TypeError from inside the guard meant to reject it.
        # Third instance of that shape in this file; see the ScalarField guards.
        not_whole = sorted(
            (n for n in self.numbers if not isinstance(n, int) or isinstance(n, bool)),
            key=repr,
        )
        if not_whole:
            raise ValueError(
                f"frame numbers must be whole numbers, but got {not_whole} — a "
                f"timeline is numbered by integer frame, and anything else "
                f"cannot be placed on it"
            )
        # A backend lowers this to a number -> name mapping, so a repeated
        # number keeps one surface and drops the other, which is then never
        # enabled on any frame. Same silent-drop the ScalarField key check
        # exists for, and the same remedy: refuse rather than render.
        if len(set(self.numbers)) != len(self.numbers):
            seen: set[int] = set()
            clashes = sorted(
                {n for n in self.numbers if n in seen or seen.add(n)},  # type: ignore[func-returns-value]
                key=repr,
            )
            raise ValueError(
                f"duplicate frame numbers {clashes} — each surface must hold a "
                f"different frame, or one of them is silently never shown"
            )
        # The timeline runs 1..max, so anything below 1 falls outside it and
        # vanishes without trace.
        below = sorted(n for n in self.numbers if n < 1)
        if below:
            raise ValueError(
                f"frame numbers {below} are below 1, and a timeline is numbered "
                f"from 1 — those surfaces would never appear. Frame numbers must "
                f"be 1 or greater."
            )


@dataclass(frozen=True)
class Morph(SceneOp):
    """Interpolate between the states of a multi-state object.

    The view emitting this has already decided the interpolation is well posed
    — that atoms pair across states — which is the judgement in
    ``morph_states``. Whether the viewer can actually perform it is separate:
    ``cmd.morph`` is Incentive-only, so on open-source PyMOL this op cannot be
    honoured even though the request was sound. That is a licence fact about
    one viewer, so the backend reports it as a note rather than the view
    hedging in prose every viewer would have to read.
    """

    name: str
    obj: str
    steps: int = 30


@dataclass(frozen=True)
class Delete(SceneOp):
    """Remove objects this view created.

    Scoped to names the view made itself, never a wildcard over the session —
    MCPymol issue #15 was a cleanup that deleted unrelated objects and reported
    success.
    """

    names: tuple[str, ...]


@dataclass(frozen=True)
class Legend(SceneOp):
    """A claim about what is being shown. Draws nothing.

    In the scene rather than only in the report text so that invariants are
    checkable on the value: I1 asks whether a scene showing a volume carries a
    provenance legend, and a report string cannot answer that without a
    substring search.
    """

    text: str
    sense: Sense | None = None
    provenance: str | None = None


@dataclass(frozen=True)
class Scatter(SceneOp):
    """A 2-D scatter of latent coordinates.

    **Defined but never emitted, deliberately.** Invariant I2 is that no view
    draws an unlabelled latent plot, and the rendering users most want from
    latent space is the one the evidence says is most often wrong: motion is
    recoverable from a heterogeneity method, populations are not.

    It exists so the invariant test can name the thing it forbids. A test that
    asserts an absence needs the absent thing to have a name, or it passes
    because of a typo.
    """

    points: tuple[tuple[float, float], ...]
    labels: tuple[str, ...] = ()


# --------------------------------------------------------------------------
# Scene
# --------------------------------------------------------------------------


@dataclass(frozen=True)
class Scene:
    """An ordered list of ops. Order matters: later ops draw over earlier ones."""

    ops: tuple[SceneOp, ...] = ()

    def __init__(self, ops=()):
        object.__setattr__(self, "ops", tuple(ops))

    def __iter__(self):
        return iter(self.ops)

    def __len__(self) -> int:
        return len(self.ops)

    def of(self, *types: type) -> list[SceneOp]:
        """Every op of the given types, in order."""
        return [op for op in self.ops if isinstance(op, types)]

    def has(self, *types: type) -> bool:
        """Does this scene contain any op of the given types?"""
        return any(isinstance(op, types) for op in self.ops)

    @property
    def legends(self) -> list[Legend]:
        return [op for op in self.ops if isinstance(op, Legend)]

    @property
    def draws(self) -> bool:
        """Does this scene change the picture at all?

        A view that refuses returns a report and a scene that draws nothing.
        That is a success, not a failure, and it needs to be distinguishable
        from a view that drew something — hence a scene of pure legends is
        explicitly *not* drawing.
        """
        return any(not isinstance(op, Legend) for op in self.ops)


class Refused(Exception):
    """A backend cannot honour an op and will not approximate it.

    Raised rather than skipped. An op quietly dropped is the failure mode this
    whole package exists to prevent: the picture looks fine and means something
    other than what was asked for.
    """


__all__ = [
    "BLUE_WHITE_RED",
    "RED_WHITE_BLUE",
    "RED_YELLOW_GREEN",
    "Arrow",
    "Arrows",
    "ColorByScalar",
    "ColorFlat",
    "ColorSurfaceByMap",
    "Colour",
    "Delete",
    "Frames",
    "Granularity",
    "Hide",
    "Isosurface",
    "Label",
    "Legend",
    "Morph",
    "Opacity",
    "Refused",
    "Rep",
    "ScalarField",
    "Scatter",
    "Scene",
    "SceneOp",
    "Sel",
    "Sense",
    "Show",
    "SizeByScalar",
    "Unit",
    "ramp",
    "resolve_colour",
]
