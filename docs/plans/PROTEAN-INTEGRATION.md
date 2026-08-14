# protean ← wiggles-em: implementation plan

Phase 5 of [`SPLIT.md`](SPLIT.md), which became the *next* phase when phase 3
(PyPI) was deferred indefinitely and phase 4 (MCPymol revert) went with it.

**Why protean and not MCPymol:** the constraint was never "the package is not
ready", it was "PyPI rejects direct URL dependencies". MCPymol is on PyPI and
therefore cannot take a git dependency at all. protean is not published, so it
can. Only one of the two consumers is blocked, and it is not this one.

**What this is really for.** Every view returns a `Scene` and never calls a
viewer; `backends/` lowers a Scene onto a specific one. That is the claim the
whole split rests on, and **only PyMOL has ever exercised it**. `backends/fake.py`
is a strictness check, not a second viewer. Until a Mol\*-shaped viewer lowers a
Scene, "viewer-neutral" is a design intention rather than a demonstrated fact.

---

## Decisions

### D1. Dependency form — settled by CI, not by preference

protean's CI runs `uv sync` after `actions/checkout@v4`, which checks out
protean and nothing else. A path dependency on `../wiggles-em` therefore
**fails CI the moment a PR is opened**. It works on a laptop where both repos
sit side by side, and nowhere else.

| option | consequence |
|---|---|
| Path dependency in a PR | CI fails immediately. Non-starter. |
| **Pinned git dependency, committed; path override locally** | CI green, works for any cloner, instant local edits retained. Cost: bump the pin when wiggles-em changes. |
| Path dependency, never PR'd, long-lived branch | Avoids the question. protean moves fast — a branch went 0 to 10 commits behind within an hour — so the branch rots. This is how the fork went stale the first time. |
| Teach CI to check out wiggles-em too | Works, but it is CI complexity propping up a dependency form that is wrong for a shared repo regardless. |

**Decision: commit the git pin, override locally.**

```toml
[project]
dependencies = ["wiggles-em", ...]

[tool.uv.sources]
wiggles-em = { git = "https://github.com/chemrich/wiggles-em", rev = "<sha>" }
```

and locally, after `uv sync`, shadow it with the working copy:

```
uv pip install -e ../wiggles-em
```

**Verified 2026-08-14, and the answer is no.** The editable install does *not*
survive `uv sync` — it is reverted to the pinned revision, and the fallback
named here does not exist either. Measured on uv 0.12.3:

```
$ uv pip install -e ../wiggles-em && uv sync
 - wiggles-em==0.1.0 (from file:///Users/charlie/code/wiggles-em)
 + wiggles-em==0.1.0 (from git+https://github.com/chemrich/wiggles-em@9381dbe7)

$ cat uv.toml
[sources]
wiggles-em = { path = "../wiggles-em", editable = true }
$ uv sync
error: The `sources` field is not allowed in a `uv.toml` file.
```

**The working local loop is two lines, not one:**

```
uv pip install -e ../wiggles-em
export UV_NO_SYNC=1
```

Both are needed. `uv run` syncs before running, so without `UV_NO_SYNC` every
test run silently puts the pin back — as a single line of output that is easy
to read past. This is written into the comment above `[tool.uv.sources]` in
protean's `pyproject.toml`, which is where someone hits it.

The branch `add-wiggles-em-dependency` (`d8400bb`) carried the **path** form and
went 15 commits behind within a day. Superseded by `wiggles-em-git-pin`
(**PR #67**), rebased onto `033aa85`; the old branch can be deleted on merge.

### D2. One PR or two

The dependency alone is **inert** — nothing imports `wiggles_em` yet. There is
nothing to review behaviourally, so the only meaningful verification is *CI
going green with the git pin*, which is exactly what wants proving before
anything is built on top. Merging it alone also lets backend work start from
`main` rather than stack on an unmerged branch, which matters at protean's pace.

**Decision: dependency alone first, backend as its own PR.**

### D3. Which view first

`occupancy_view` needs only `load_structure` + `color`. That is the minimum
that exercises the mechanism `scene.py` predicted and nobody has tested: a
per-atom scalar riding the B-factor column of a re-sent display copy, read by
Mol\*'s `uncertainty` theme.

`altloc_view` is more visually interesting but adds `hide`/`show`/`label`
without testing anything new about the seam. Volume views were blocked on the
`cryoem-volumes` branch, which is a protean gap rather than a wiggles-em one;
it **merged as protean PR #69**, and `Isosurface` now lowers (PR 73).

**Decision: `occupancy_view` end to end, then `altloc_view`, then
`ensemble_spread_view`.**

### D4. Ops protean cannot honour

`Frames`, `Morph` and `Arrows` have no protean equivalent. The design says
refuse rather than approximate, and `Refused` exists for exactly this.

**Decision: refuse explicitly, with a message naming what is missing.**

**The price, stated plainly:** `Frames` is what `latent_traverse_view` needs,
and a latent traversal is arguably the most valuable cryo-EM view protean could
offer. Refusing means **protean cannot do latent traversals** until a movie
timeline exists in Mol\*. That is a protean roadmap item, not a wiggles-em one.

### D5. protean going public, and PyPI later

Public plus a path dependency is broken for every cloner — the same failure as
CI, for the same reason. Public plus a git dependency on wiggles-em is fine,
since wiggles-em is public.

**Dormant, not dead:** if protean is ever published to PyPI it inherits
MCPymol's exact problem, and a git-pinned wiggles-em becomes unshippable. The
PyPI question returns the moment protean wants to be pip-installable. Nothing to
decide now; worth not being surprised by.

---

## Op → action mapping

**Corrected 2026-08-14 against a written backend** (PR #68). The first version
of this table was read off protean's `_call(...)` sites — which name the
*actions* correctly and say nothing about whether the arguments those actions
take can carry what an op means. Three rows were wrong in exactly that way, and
all three were "available".

| Scene op | protean action | status |
|---|---|---|
| `ColorFlat` | `select` + `color` | available |
| `ColorByScalar` | `load_structure` + `show` (uncertainty theme) | available |
| `SizeByScalar` | — | **refuse** — see below |
| `ColorSurfaceByMap` | `color_by_volume` | **refuse** — needs a volume colour theme |
| `Opacity` | `select` + `opacity` | available |
| `Show` | `select` + `show` | available |
| `Hide` | `select` + `hide` | available **only for `Rep.EVERYTHING`** |
| `Label` | — | **refuse** — see below |
| `Delete` | `remove` | available |
| `Legend` | — (report text, no viewer call) | n/a |
| `Isosurface` | volume action + `unit` | available (protean PR 73); **refuse** with a carve |
| `Frames` | — | refuse |
| `Morph` | — | refuse |
| `Arrows` | — | refuse |
| `Scatter` | — | refuse (forbidden by I2 anyway) |

**`SizeByScalar` — the bridge exposes no size *theme*.** `show` takes a scalar
`size`, which becomes Mol\*'s `sizeFactor`: one number for the whole
representation, not a per-atom ramp. Mol\* itself *has* an `uncertainty` size
theme — the size counterpart of the colour theme this integration is built on —
and the viewer simply does not select it. The op's own contract forbids falling
back to colour, so it refuses.

**`Label` — the bridge's `label` takes no text.** It draws structural labels at
a level of `chain`, `residue` or `element`. A `Label` op carries literal text
plus atom fields to interpolate (`"B %.2f"` with `fields=("q",)`), and there is
nowhere to put it.

**`Hide` hides a component whole.** A specific representation cannot be hidden
while its siblings stay drawn, so anything but `Rep.EVERYTHING` refuses rather
than over-hiding — which would remove more of the picture than was asked for.

Ops emitted by the first three target views, and what that now means:

    occupancy_view         ColorByScalar, Show, Legend        renders in full
    altloc_view            Hide, Show, ColorFlat, Label       renders with
                                                              label=False
    ensemble_spread_view   ColorByScalar, SizeByScalar        BLOCKED on a
                                                              size theme

So **step 3 is unaffected** and step 5 splits: `altloc_view` is a matter of
passing `label=False`, while `ensemble_spread_view` needs a protean change
first — teaching `show` to select a size theme. That is a small, well-scoped
addition to `dispatch.ts`, and it is a protean roadmap item rather than a
finding about the seam.

**What this cost, and the lesson.** Reading `_call` sites answered "does an
action with this name exist", and the question that mattered was "can its
arguments carry what the op means". Two different questions, one of which looks
like the other. A vocabulary audit is not an interface audit.

---

## Steps

1. ~~**Convert the branch to a git pin** and open the dependency PR. CI must be
   green. Confirm the local editable override survives `uv sync`.~~ **Done
   2026-08-14 — PR #67, merged as `a679690`.** Rebased onto `033aa85`, pinned to wiggles-em
   `9381dbe7`. Verified in a fresh clone with no sibling checkout, which is the
   shape CI runs in: `uv sync` installs from git, and ruff, ruff format, mypy
   strict and pytest (587 passed, 292 skipped) are green. All three CI jobs
   passed on the PR, including the real-browser differential job (14m2s) — the
   one a path dependency could never have reached. The editable override does
   **not** survive `uv sync`; see D1 for what does.

   **Step 2 starts from `a679690`, not from a branch.** That was the point of
   D2 — merging the inert dependency alone means the backend work begins on
   `main`, which matters at protean's pace.
2. ~~**`backends/molstar.py` in protean**, mirroring
   `wiggles_em/backends/pymol.py` in shape: a class taking whatever protean
   uses to reach the viewer, a `render(scene)`, one private method per op, and
   `Refused` for the four it cannot honour.~~ **Done 2026-08-14 — PR #68.**
   Seven refusals, not four; see the corrected table above. 31 offline tests,
   each guard checked by breaking it. `render` is `async`, which is the one
   shape difference from the PyMOL backend and follows from the bridge being
   async.

   **Two findings worth carrying:**

   - **The views need no port.** `occupancy_view` takes `list[Atom]` and
     returns a `Scene`; nothing implements `PymolPort`. `atoms_for()` reads a
     biotite `AtomArray` into those atoms directly. The plan had budgeted for a
     transport adapter that turns out not to exist — the extraction from
     MCPymol left the atom-based views genuinely viewer-free.
   - **`Colour` named one viewer's palette. ~~Open.~~ CLOSED 2026-08-14**, in
     wiggles-em PR #3 (`df2252e`). `Colour = str | tuple`, and every string a
     view emitted was a *PyMOL* colour name — `grey70`, `skyblue` — which a
     second viewer could only honour by reimplementing PyMOL's table, which is
     what `_COLOUR_NAMES` in the backend is. That was a gap in the seam, not in
     the viewer: a viewer-neutral value should not name a viewer.

     `scene.resolve_colour()` is now the single table and views resolve at the
     op site, so **a Scene carries only RGB**. Names remain valid *arguments* —
     a view's signature is not the seam. `ColorByScalar.palette` went the same
     way: it was a PyMOL *spectrum* name (`blue_white_red`) and is now colour
     stops, low value first, built by `ramp()`.

     Two corrections to what this section used to say. `latent_traverse_view`
     did **not** already emit RGB; it emitted its `color="skyblue"` default.
     And the values themselves were wrong — `skyblue` is PyMOL's
     `(0.2, 0.5, 0.8)`, not `(0.34, 0.63, 0.83)`, and `lightblue` is
     `(0.75, 0.75, 1.0)`, not a pale cyan. Both wrong values came from the two
     tables being consolidated, which **agreed with each other**; that read as
     corroboration and was one error copied twice.
     `test_the_palette_matches_pymols_own` now asks a real PyMOL instead.

     **protean must act on this at the pin bump**, and it fails open until it
     does: `op.palette != "red_white_blue"` compares against a tuple now, so
     every scene gains a "was not applied" note and the two ramp directions
     stop being distinguishable. Its note also claims "the ordering and the
     domain are honoured", which a single-direction `uncertainty` theme cannot
     deliver for both ramps — that is the real defect, and stops make the
     direction inspectable so it can refuse instead of drawing backwards.
     Separately, protean's own `_COLOUR_NAMES` carries the same two wrong
     values, which is live today regardless of the pin.

   Also: **wiggles-em ships no `py.typed`**, so it is fully annotated and
   entirely invisible to mypy. protean carries an override for now; the fix is
   one file upstream plus a pin bump.
3. **`occupancy_view` end to end.** The interesting part is `ColorByScalar`:
   build a display copy, write the scalar into the B-factor column, re-send via
   `load_structure`, colour with the `uncertainty` theme. protean does **not**
   need the stash/restore dance PyMOL needs — it has a display copy, not the
   user's only copy.
4. **Hands-on check in the browser.** Not a test: does the picture look like
   occupancy, is the ramp the right way round, does the legend match what is
   drawn. This is the step the whole plan exists to reach.
5. **`altloc_view`, then `ensemble_spread_view`.**
6. **Refusals, with tests.** A scene containing `Frames` must refuse with a
   message naming the missing capability, not render three-quarters of itself.
7. **Re-point `check_spec.py`** (phase 6) once two consumers exist, with
   per-consumer coverage as a separate, weaker check — a tool protean cannot
   host is a known gap, not spec drift.

## What protean asks of wiggles-em

From the protean side's handover, 2026-08-14. These are changes *here*, and the
first one blocks something rather than merely tidying.

1. **`to_sigma` / `to_absolute` should take measured statistics, not a
   `MapHeader`** (`density.py:87`, `:103`).

   protean found that Mol\*'s `grid.stats` — the four numbers a viewer reports
   for a volume — are, for CCP4/MRC, *stored header fields* passed straight
   through unexamined. A fixture written with deliberately false header
   statistics (−999 / 999 / 42 / 7) failed on its first run with `min came back
   as the header's false value -999.0` **and the dimensions correct**, so the
   volume had genuinely parsed and every number describing it was the file's
   claim rather than its contents. protean now walks the voxels and reports the
   header's own four numbers separately under `stated`.

   These converters take a `MapHeader`, so they convert against exactly the
   numbers shown to be unreliable, and protean holds a trustworthy sigma it
   cannot hand them. Taking a small stats value instead would let any host feed
   viewer-measured statistics.

   How far it reaches: **Mol\*'s own default isosurface is `relative: 2`**,
   computed as `relativeValue * grid.stats.sigma + grid.stats.mean`. Any viewer
   contouring a map with a stale header puts the surface in the wrong place and
   it looks entirely normal. This is the σ trap this package already documents,
   one layer further out than the version we wrote down.

2. **`Isosurface.equivalent` is unused by this backend.** Its docstring says a
   backend "converts against the volume's header, which it reads from the
   `load_map` record", and that ensemble frames have no such record — which is
   what `equivalent` routes around. protean has no such constraint: it reaches
   computed statistics for any loaded volume, frame or not. Not a request to
   remove it, since MCPymol may still need it, but **"a backend converts against
   the header" is no longer true of every backend** and the contract should say
   so.

3. **The pin is stale, and the colour change makes it urgent.** protean pins
   `9381dbe`; `main` is well past it. Nothing is broken — protean builds against
   the pin and its CI is green — but see the `Colour` note in step 2: the palette
   change fails *open* on protean's side once the pin moves, so the bump and the
   `molstar.py` guard fix belong in one change.

---

## Cautions carried from the wiggles-em rounds

- protean verifies rendering by **reading pixels**, not return values. A backend
  test that asserts the action was *sent* proves the same thing MCPymol's mocked
  suite proved, which is not much. Assert on the canvas where it matters.
- The seam's claim is that a Scene needs no PyMOL knowledge. If the Mol\*
  backend ends up needing something only PyMOL provides, that is a finding about
  the seam and should be recorded as one rather than worked around in the
  backend.
- protean has **no docs-enforcement test**; its README tool count is
  hand-maintained. Adding tools without noticing is easy here in a way it is not
  in MCPymol.
