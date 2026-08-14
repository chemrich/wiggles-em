# Splitting the code out again — `wiggles-em`

Planned 2026-08-12. Supersedes the arrangement in [`MOVING.md`](MOVING.md),
which is now a record of the first move rather than the current design.

`MOVING.md` put the code inside MCPymol because there was one viewer. There
are two now — MCPymol drives PyMOL, [protean](https://github.com/chemrich/protean)
drives Mol\* in a browser — so the code becomes its own package and both
consume it.

**The compendium does not move.** Research private, tools public was the right
split and still is; this only changes where "tools" lives. Ten open items in
[`PUBLICATION-GATE.md`](PUBLICATION-GATE.md) gate the prose, and a package two
public repos depend on has to be public and on PyPI. Hence three repos rather
than two.

---

## Why now, specifically

`v1.5.1` is MCPymol's latest tag and release, and **no tag contains `992edd7`**.
PR #54 landed after the release was cut, so no published MCPymol has ever
carried these tools and nobody has installed one that does. Reverting is free
today. The moment 1.6.0 ships it becomes an API removal with a deprecation
cycle attached.

That window is the whole reason to do this now rather than after tier 3 lands.

---

## The three repos

| Repo | Visibility | Holds |
|---|---|---|
| `chemrich/wiggles` (this one) | private | compendium, `SPEC.md`, publication gate, `check_gate.py`, `check_spec.py`, `livefire.py` |
| `chemrich/wiggles-em` | public, PyPI | the package: analysis, views, `Scene`, `FakeBackend` |
| `chemrich/MCPymol`, `chemrich/protean` | public, PyPI | a backend + `@mcp.tool` wrappers each |

`wiggles` is taken on PyPI (a signal-processing library, v2.5.8), so the
distribution and import name are both `wiggles_em`.

---

## The seam moves up a layer

`SPEC.md` §0 says Wiggles talks to PyMOL through one small protocol and nothing
else. That protocol abstracted the **transport**, and it did that job well —
protean's viewer bridge is the same shape (`{id, action, args}` over a
WebSocket against `{action, args, kwargs}` over TCP). What does not port is the
**vocabulary**: `call(port, "spectrum", "b", "red_white_blue", obj, minimum=0,
maximum=1)` names a `pymol.cmd` function, and Mol\* has never heard of it.

So §0's rule is replaced, not weakened:

> **Rule: a view computes and returns. It never calls a viewer.**
> Views are pure functions returning `(report, Scene)`. A backend lowers a
> `Scene` onto one viewer's API. Nothing above `backends/` imports a viewer,
> a socket, or a command string.

```python
def occupancy_view(atoms: list[Atom], obj: str) -> tuple[str, Scene]:
    ...
    return report, Scene([
        ColorByScalar(Sel.obj(obj), values=qs, domain=(0.0, 1.0),
                      palette="red_white_blue"),
        Show(Sel.obj(obj) & Sel.lt("q", FULL_OCCUPANCY), Rep.STICKS),
    ])
```

### Why declarative rather than a wider `Viewer` protocol

Three reasons, in order of how much they matter.

**1. It puts the σ trap where it belongs.** Every view currently holds an
absolute contour level and converts to σ against the right map's header,
because PyMOL normalises each map independently. Mol\* does not need that: it
has `Volume.IsoValue.absolute()` and `.relative()` and converts internally. If
a view emits `Isosurface(map, level=0.05, kind=ABSOLUTE)`, then converting is
a **property of the backend** — PyMOL's converts, protean's passes through —
instead of a rule every view has to remember. This is the single hardest-won
finding in the project and it stops being something a view can get wrong.

**2. Invariants become assertions on a value.** I1/I2/I3 are claims about what
gets drawn. Today they are checked by grepping a recorded command log
(`FakePort.ran("spectrum")`). Against a `Scene` they are checks on a data
structure: *this scene contains no `Scatter` op*, *every scene carries a
`Legend` naming its provenance*. A string match passes for the wrong reasons;
this does not.

**3. Views become testable with no viewer and no fake.** The 209 `FakePort`
tests already run without PyMOL, but they assert on I/O that was recorded.
Asserting on a returned value is smaller and says more.

The cost is real: all ten views change shape, and any view that needs to read
state *between* draws cannot. Checked all ten — every one reads up front and
draws at the end, so nothing is lost. If a later view needs interleaving, it
takes a reader argument and still returns a `Scene`.

### The op vocabulary

Extracted from every `call(port, …)` site across the ten tools. It is small.

| Op | Carries | PyMOL lowering | protean lowering |
|---|---|---|---|
| `ColorByScalar` | selection, per-atom values, **explicit domain**, palette | `alter b=` + `spectrum` | `color_by_scalar` |
| `ColorFlat` | selection, colour | `color` | `color` |
| `Show` / `Hide` | selection, representation | `show`/`hide` | `show`/`hide` |
| `Label` | selection, per-atom text | `label` | `label` |
| `Isosurface` | map, level, **kind: absolute or sigma**, carve radius, mesh/surface | `isomesh`/`isosurface` + conversion | `isosurface`, `IsoValue` direct |
| `ColorSurfaceByMap` | surface, map, breakpoints, palette | `ramp_new` + `set surface_color` | extend `color_by_volume` |
| `SizeByScalar` | selection, values, domain | `cartoon putty` + scale min/max | `cartoon` + `uncertainty` size theme |
| `Legend` | text, provenance | report text | report text |

Every op has a confirmed lowering on both backends.

### Selections

Scene ops cannot carry PyMOL selection strings — protean parses a subset with
real gaps. So `Sel` is a tiny algebra the backends lower:
`Sel.obj(name)`, `Sel.prop("name", "CA")`, `Sel.lt("q", 0.999)`,
`Sel.indices([...])`, and `&` / `|` / `~`.

Five node types covers every selection the ten tools build. `Sel.indices` is
the escape hatch — wiggles has already fetched the full atom table by the time
it draws, so it can always compute a set itself; backends that can express the
predicate compactly should, because `index 1+2+3+…` over 100k atoms is a
command line nobody wants.

---

## What lives where

```
wiggles_em/
  analysis/          pure; no viewer, no I/O beyond reading files
    mapinfo.py       MRC/CCP4 headers, voxel size (cella/m), isotropy tolerance
    provenance.py    Provenance, banners, I1
    qscore.py        wwPDB validation-report parsing
    density.py       sigma <-> absolute, carve defaults
    atoms.py         Atom, altloc_groups, group_by_residue
    ensembles.py     spread and RMSF maths
  scene.py           Scene, the ops, Sel, Rep, palettes
  views/             analysis + scene + report text; the invariants live here
  sources.py         StructureSource protocol: atoms(sel), n_states(), coords()
  backends/
    fake.py          records a Scene; the whole test suite
    pymol.py         lowers Scene onto pymol.cmd; owns the sigma conversion
                     and the B-factor stash hack
```

### Two things stop being shared

**`bfactors.py` moves into the PyMOL backend** — but not for the reason it
first appears. Mol\* has the *same* constraint PyMOL does: its `uncertainty`
size theme reads `B_iso_or_equiv` and that is the one per-atom numeric field it
will ramp over, so protean smuggles scalars through the B-factor column too.

What differs is the destructiveness. PyMOL has one copy of an object, so
`alter b=q` overwrites the crystallographic values and they have to be stashed
and restored. protean builds a **display copy** and re-sends it, leaving the
analysis copy intact — see `color_by_rmsf` in its `server.py`. So the op is
`ColorByScalar(values, domain, palette)` either way; stash-and-restore is the
PyMOL backend's lowering detail, and `preserve_bfactors` is a PyMOL-only
parameter that should not appear in protean's tool schema.

**`cmd.morph` being Incentive-only stops being a shared caveat.** It is a PyMOL
licensing limitation. protean has `load_trajectory` and `frame` and does not
inherit it, so `morph_states`' loudest disclaimer is backend-specific text.

### Where the reads come from

`PymolPort.query` is doing two unrelated jobs today: issuing draws, and reading
per-atom data that has no other route out of PyMOL. Only the first is a viewer
concern. In protean the atoms are already in Python — biotite parsed the file
server-side — so a viewer round-trip to read occupancy would be absurd.

Hence `StructureSource`: `atoms(selection)`, `n_states()`, `coords(state)`.
MCPymol's implementation is the existing `iterate_to_list` call. protean's
reads biotite. Views take a source and return a `Scene`; neither half touches a
socket.

---

## Update, 2026-08-12: the fork went stale in a day

Between the extraction and now, MCPymol's copy went through a code review that
landed **PR #57** ("the four data-destroying defects"), **PR #58** ("the
remaining review findings, and close the gap that hid them"), and the **merge
of #55**. `origin/main` is at `0d4c9a5`; the extraction was taken from `#55`'s
branch head at `ea645e4`, which predates all of it.

That is **711 insertions across 19 files** that `wiggles-em` does not have, and
an audit says it has **none of the nine fixes**. Worse, one of them it
reproduced independently: `render_selection` interpolates chain and residue
identifiers unquoted, which is exactly the blank-chain bug — and the
chain-grouping change made for the 1123-term disjunction sits directly on top of
it, so a negative residue number is now ambiguous in a second way.

**The lesson for this plan, not just for the code:** a fork of a moving,
actively-reviewed module diverges immediately, and re-deriving the same
mistakes is the expected outcome rather than bad luck. So the re-sync is not
"apply nine patches" — it is **port the upstream tests first**, and let a fix
that fails to arrive show up as a failing test.

### What has to be ported, and where each lands

The port is not mechanical, because the seam moved. Upstream's fixes sit where
the old architecture put them; three of them belong somewhere else now.

| Upstream fix | Where it lands in `wiggles-em` |
|---|---|
| `usable_rms` — MRC writes `rms=-1` for "statistics not computed", and it divides cleanly, so `to_sigma(0.05)` returned `-2.05` and a resolution ramp ran backwards | `analysis` — unchanged. **`to_absolute` needs the guard too**; it currently has none at all |
| Provenance: longest matching token wins, so `unsharpened` stops being read as `sharp` | `analysis` — unchanged |
| `loaded_map(obj, port)` evicts a record whose object has left the session | loader layer, which still holds a port — fits |
| First stash wins; a restore clears it | `analysis`. The backend already refuses to re-stash, but the guard belongs in `stash_bfactors` where every caller gets it |
| **`quote` / `residue_selection` / `residue_clause`** | **the PyMOL backend, not `atoms`.** This is PyMOL's selection grammar, and no other viewer shares it. `render_selection` is the single place that names a residue now, so it is the right home — but the quoting has to survive the `+`-list grouping |
| Contour follows the normalisation setting | **already present**, arrived at independently |
| `normalisation_state` moved to `maps` | **already present**, in the backend instead. Keep it there: it is a question only PyMOL can be asked |
| `latent_traverse_view` anchors on the first frame with a usable RMS, not frame 0 | view — unchanged |
| `ensemble_spread_view` superposes states before measuring | **needs a decision — see below** |

### The one that does not port cleanly

Upstream fixes `ensemble_spread_view` by calling `intra_fit`, because spread
measures whatever separates the states and a rigid-body drift is not
flexibility. That is right, and the finding holds for any viewer.

But `intra_fit` mutates the session, and `wiggles-em`'s view no longer talks to
a session — it receives coordinates. Three options, none free:

1. **The host fits before reading.** PyMOL calls `intra_fit`, protean uses
   biotite's `superimpose`. The view takes `superposed: bool` and says which
   happened. Honest, and it keeps the coordinate mutation where the session is.
2. **The package fits, purely.** It has the coordinates; Kabsch would make the
   view self-contained and mutate nothing. But Kabsch needs an SVD, and this
   package has **no runtime dependencies** on purpose — that is a load-bearing
   property for two consumers.
3. **The view refuses** unless told the states share a frame. Safe, useless.

Option 1, unless the no-dependency rule is worth trading. Worth deciding before
the port rather than during it.

### Two testing gaps upstream found, which apply here

`#58`'s most valuable finding was not a bug but the reason the class of bugs
survived: the live sweep supplied arguments from its own table, so it **never
called a tool the way its schema said it may**, and `spheroid` passed every run
with a default that always failed. `wiggles-em`'s suite has the same shape —
every test passes explicit arguments, so no documented default is ever
exercised. That test comes across.

The second: `WIRING_FAILURES` caught only Python binding errors, so "the viewer
accepted the signature and rejected the argument's *meaning*" passed silently.
`FakePort` has exactly that blind spot — it answers `"OK"` to any command it
does not recognise.

---

## Phases, and what constrains the order

**Revised.** A new phase 0 goes first, and the MCPymol revert moves behind it:
reverting a reviewed, fixed module so it can depend on a fork that is behind it
would be a downgrade, whatever the packaging argument says.

**0. Re-sync `wiggles-em` with `origin/main`.** Port the nine fixes and the
~380 lines of tests that came with them, placing each per the table above.
Tests first, so a fix that fails to arrive fails a test.

**The tier-3 code exists only on an unmerged branch.** No longer true — #55
merged as `0d4c9a5`. The constraint it created is gone, and `origin/main` is
now the better extraction source than anything this plan has been using.

1. **Stand up `wiggles-em`.** Port tiers 1–2 from MCPymol `origin/main` and
   tier 3 from the PR #55 branch. Straight lift first — same `port.py`, same
   views — with the suite passing. No refactor yet; a move and a redesign in
   one step means a failure could be either.
2. **Introduce `Scene` and convert the views.** Backend `fake.py` and
   `pymol.py`. The suite converts from asserting on recorded calls to asserting
   on returned scenes. `SPEC.md` §0 is rewritten here.
3. **Publish `wiggles-em` to PyPI.** ~~MCPymol cannot depend on it before
   this.~~ **DEFERRED INDEFINITELY — Charlie's call, 2026-08-13.** No date, no
   trigger condition. The repo is public
   (github.com/chemrich/wiggles-em) and that is where it stays for now.

4. **MCPymol: revert, then re-add thin.** **BLOCKED, and not merely delayed.**
   PyPI rejects uploads whose dependencies use direct URL references, so a
   `wiggles-em @ git+https://…` dependency cannot ship in a PyPI package.
   MCPymol is on PyPI. With phase 3 deferred indefinitely there is therefore no
   route to this phase at all, and the revert must not be attempted: it would
   delete a working, four-times-reviewed module and replace it with a
   dependency MCPymol cannot declare.

   **MCPymol keeps its own copy indefinitely.** That is the decision, not an
   oversight, and it has a running cost — see below.

5. **protean: the volume tools, then its backend and wrappers.**
   **This is now the next phase, and it is unblocked.** protean is private and
   not published to PyPI, so it *can* depend on
   `wiggles-em @ git+https://github.com/chemrich/wiggles-em@<sha>` today. The
   asymmetry is the whole point: the constraint was never "the package is not
   ready", it was "PyPI forbids URL dependencies", and only one of the two
   consumers is on PyPI.

6. **Re-point `check_spec.py`.** It reconciles `SPEC.md` against what *MCPymol
   registers* today. It should reconcile against what **`wiggles_em` exposes**,
   with per-consumer coverage as a separate, weaker check.

---

## The cost of keeping two copies, and what to do about it

With phase 4 off the table, `mcpymol.wiggles` and `wiggles_em` are two
maintained copies of the same code. This is the exact condition that produced
the "fork went stale in a day" entry above, and it is no longer a transient
state to be endured — it is the steady state.

**It has already cost us twice, in one day.** Two defects were found in
`wiggles-em`, fixed there, and then found to exist in MCPymol and fixed *again*
as separate PRs:

| defect | wiggles-em | MCPymol |
|---|---|---|
| quoting does not escape a negative `resi` | `2ee28fc` | #59, `0031b93` |
| provenance lost its category priority | G3 + H4 | #60, `8357e02` |

Neither was caught by MCPymol's own review rounds (#57, #58). Both were found
only because the extracted copy was reviewed hard, and both had shipped to PyPI.

So the extraction is earning its keep as a *review surface* even while it cannot
be consumed as a dependency. That is worth stating plainly, because "the split
did not achieve its packaging goal" and "the split was not worth doing" are
different claims and only the first is true.

**The discipline this requires, until phase 4 becomes possible:**

- A fix to a view, a selection, a unit conversion or a provenance rule lands in
  **both** repos or in neither. Not "port it later" — later is how the fork went
  stale the first time.
- The port is not a cherry-pick. The two have diverged structurally: wiggles-em
  has the `Scene` seam and MCPymol does not, so `quote` lives in
  `backends/pymol.py` in one and `atoms.py` in the other. Port the *test* first,
  then make it pass in that repo's own idiom.
- Divergence is checkable and nobody is checking it. The audit that belongs in
  CI: every public name in `mcpymol.wiggles` has a counterpart somewhere in
  `wiggles_em`. Run 2026-08-13 — 69 public names upstream, 102 in the package,
  two without counterparts (`residue_selection`, `residue_clause`, both
  deliberately replaced by `Sel.residues` + `render_selection`).


## What protean needs building

Worth building on their own merits, not as a favour to wiggles:
`color_by_volume` currently parses **OpenDX only** (APBS output), while Mol\*
4.18 ships providers for `ccp4`, `dsn6`, `dx`, `cube`, `dscif` and `segcif`.
Anyone showing a cryo-EM or crystallographic map hits that gap with or without
wiggles.

| Tool | Notes |
|---|---|
| `load_volume(source, format=auto)` | The `rawData` → parse → state-tree path is already proven by `color_by_volume` (`dispatch.ts:1548`). Mostly format plumbing and gzip. |
| `isosurface(volume, level, kind=absolute\|sigma, carve_around=None, style=mesh\|surface)` | `kind` maps onto `Volume.IsoValue.absolute`/`.relative`. Exposing both explicitly is the single best defence against contouring noise. |
| `color_surface_by_volume(surface, volume, domain, palette)` | Local resolution. `color_by_volume` retargeted from a structure to a volume representation. |
| `volume_info(volume)` | Dims, voxel size, min/max/σ — **read back off the viewer**, per protean's own rule that replies report state rather than echo arguments. |

### The altloc problem

`selections.py:237` lists `alt` as unsupported: *"coordinates are parsed keeping
one conformer per atom site, so no altloc field survives to select on. Loading
every conformer is possible but would make buried areas and potentials be
computed over atoms that overlap each other."*

The reason is correct and the conclusion is too broad. It argues for the
**analysis** path deduplicating conformers; it does not argue for discarding
the field at parse time. As it stands `altloc_view` cannot be hosted by protean
at all, and `occupancy_view` shows a model whose partial-occupancy atoms have
already been silently thinned.

This is the compendium's own thesis pointed at protean — the entry is
`multiconformer`, and "a viewer threw the alternates away and the picture
looked fine" is exactly the failure the project catalogues. Resolution: keep
altlocs in the parsed table, deduplicate at the point of use in buried-area and
electrostatics, and let `alt` become selectable. That is protean's call and it
is not small.

---

## Spikes — two resolved, one open

- ~~**`SizeByScalar` on Mol\*.**~~ Resolved: `mol-theme/size/uncertainty.js`
  computes `baseSize + B_iso_or_equiv * bfactorFactor`, so cartoon plus the
  `uncertainty` size theme over a display copy is a putty. Still worth
  confirming it applies to `cartoon` and not only to ball-and-stick.
- ~~**Carve radius on Mol\*.**~~ Resolved, by not being a viewer feature.
  Mol\*'s `selection-box` lives only in the volume-streaming behaviour, which
  wants a remote density server. protean has the map and biotite server-side,
  so carve is a numpy crop with a corrected origin before the bytes are sent —
  which also makes the transfer smaller. It belongs on `load_volume`, not on
  `Isosurface`, because it changes what data exists rather than how it is drawn.
  **This makes `carve` an op parameter with no protean lowering at all**, which
  is the first case of a Scene op that one backend answers by preprocessing.
- **Backend parity test.** `interpret()` had a test asserting both adapters
  translate identically, because divergence would make behaviour depend on how
  the package was installed. The same argument now applies to backends, and it
  is a harder test: the same `Scene` must mean the same picture. At minimum,
  every backend must reject the ops it cannot honour rather than skipping them
  silently.

---

## Decisions carried forward unchanged

Nothing here reopens `occupancy_view` vs `composition_view` (separate forever,
neither infers the other), provenance defaulting to UNKNOWN and never being
inferred, structured actions over `do()`, or motion-recoverable /
populations-not. The seam changes; the argument does not.

One is *strengthened*: "structured actions, never `do()`" was about `cmd.do`
reporting success for literal nonsense. A `Scene` has no `do()` to reach for.
