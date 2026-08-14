# Changelog

All notable changes to this project are documented here.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added

- The package, extracted from MCPymol's `mcpymol.wiggles`. Fourteen views
  across three tiers: occupancy and ensembles, maps and local resolution, and
  heterogeneity-method output.
- No runtime dependencies. Two hosts depend on this package, so a pin here
  would become a constraint in both of them.
- `wiggles-map-info`, a CLI over `map_info` — `mapinfo.py` was already written
  with a `main()`, and reading a real deposition's header is the quickest check
  that a map is what somebody says it is.

- `wiggles_em.scene` — the viewer-neutral seam. Twelve ops, a five-shape
  selection algebra, and `ScalarField`. Views return `(report, Scene)` and
  never call a viewer.
- `wiggles_em.backends` — `PymolBackend` and `FakeBackend`, the only code in
  the package that knows a viewer exists. A backend refuses an op it cannot
  honour rather than skipping it: a dropped op leaves a picture that looks fine
  and means something else.

- CI, in `.github/workflows/ci.yml`. Ruff and mypy once; the suite on Python
  3.10, 3.11, 3.12 and 3.13 — the four versions the classifiers advertise and
  that nothing had ever tested. Three guards, each mutation-tested against the
  failure it exists to catch: the matrix job asserts it got the interpreter it
  asked for, the collected count must equal 543 exactly so a file that
  accidentally gains a `network` or `live` marker cannot vanish behind a green
  tick, and `uv run pytest` must collect the same count as
  `uv run python -m pytest`, which is the G9 regression.

  Only the marker case is *silent*; a broken `pythonpath` exits non-zero and
  says so. The guards are written to survive that difference — they check
  pytest's own exit status before counting, because a count taken through a
  pipe reports `grep`'s status and would let a collection error read as a
  plausible-looking number.

- `ruff format --check` in CI, and the 16 of 44 files that had drifted are
  formatted. The README never listed `ruff format`, so it had never run here
  and nothing was enforcing it. Behaviour-neutral: the suite reports the same
  524 passed, 19 skipped, 11 deselected before and after.

### Changed

- **The package no longer hosts its own tools.** `tools.py` held MCPymol's
  `@mcp.tool` wrappers and did not come across; the views are the public
  surface and a host wraps them. `__init__` re-exports each view from the
  module it lives in rather than from the tool layer.

- **Views take data, not a port.** `occupancy_view(atoms, obj)` rather than
  `occupancy_view(port, obj)`. Reading atoms is the host's job — in MCPymol
  that is an `iterate_to_list` round trip, in protean the atoms are already in
  Python because biotite parsed the file.
- **Loaders are not views.** `load_map` and `load_ensemble` issue a load and
  then query back to confirm the object arrived, which a one-shot `Scene`
  cannot express, so they keep taking a port.
- **Unit conversion moved into the backend.** Views state levels in the unit
  the data is in and tag it; `PymolBackend` converts to σ against that map's
  own header. Mol\* takes absolute levels natively and will not convert at all.
- **The B-factor stash moved into `PymolBackend`**, along with the note about
  it. Both viewers route per-atom scalars through the B-factor column, but
  PyMOL has one copy of an object and must stash the originals while protean
  re-sends a display copy — so the "call `restore_bfactors`" advice is PyMOL's,
  not occupancy's.
- `local_resolution_view` lost `ramp_name`. A ramp is a PyMOL object; no other
  viewer has one to name.

### Fixed

Four findings from a code review of the seam conversion, all of them behaviour
that lived inside a view and was dropped or split when the view stopped calling
the viewer. Each is silent — the render looks ordinary and means something else.

- **Atom identity is `(model, index)`.** It was chain + residue + name +
  altloc, which collides on PDB insertion codes (`bfactors._key` documented
  this and it was carried into the new scalar path) and across a selection
  spanning two loaded models. A collision handed one atom another's occupancy
  and drew it a legitimate colour on the right scale. `ScalarField` now also
  refuses duplicate keys outright, so no future key choice can collapse
  quietly.
- **`normalised` is read once and passed to both the view and the backend.**
  The view took it as an argument for its report while `PymolBackend` queried
  the session independently, so the colour key could list σ breakpoints while
  the surface was ramped in Ångström. It is now a required argument on
  `local_resolution_view` — a default is what let the two drift.
- **A contour level converts to the unit the viewer actually wants.** With
  `normalize_ccp4_maps` off PyMOL reads an isosurface level as an absolute map
  value, and the backend was converting to σ regardless — contouring far above
  dmax and yielding an empty surface while the report stated the level it
  believed it had asked for.
- **Latent surfaces keep their own frame numbers.** Numbering ran over the
  frames that survived the rms=0 filter, so a single skipped frame made
  `_03` hold frame 4's density and the user read each density against the wrong
  latent coordinate. The report now names which frames were skipped.

The remaining six from the same review. Three are contracts the seam left
unstated, and the fix for each is to remove a default rather than document one.

- **`Sel.residues` lowers to one term per chain**, using PyMOL's
  `+`-separated residue list, rather than one parenthesised `or` term per
  residue. 1123 scored residues on 9C0K was a ~30 KB selection PyMOL evaluated
  term by term, sent three times, past the 10 s port timeout. The list is never
  collapsed into ranges — turning 5 and 7 into 5-7 would silently add 6.
- **The movie is entered.** `mdo` commands run when a frame is entered, and the
  `frame 1` that did so was dropped when the wiring moved into the backend, so
  every isosurface stayed enabled and a traversal rendered as one superimposed
  blob.
- **Frame wiring names every sibling explicitly** instead of disabling a
  reconstructed `prefix_*` glob, which also switched off an unrelated
  `v_model` in the session on each frame step.
- **`deformation_view` requires `start_state` and `end_state`.** The view no
  longer reads coordinates, so nothing can check the numbers against the
  arrays, and an `end_state` defaulting to the last state described a 1→20
  transition over a 1→2 displacement.
- **`composition_view` takes a `count_atoms` callable** and hands it the exact
  `Sel` it is about to colour. A pre-counted dict keyed on the table's own text
  let a host answer for a differently-scoped selection — `chain B` in a session
  with a second structure loaded — and the guard passed on atoms that were not
  in the object at all.
- **`Sel.first` replaces `rank 0` for labels.** An atom index is numbered over
  the whole object, so ANDing `rank 0` with a per-part selection matched
  nothing for every part but one, and no percentages were ever drawn.

### Fixed — re-synced with MCPymol `origin/main`

The extraction was taken from PR #55's branch head. MCPymol's copy then went
through a review that landed **#57** ("the four data-destroying defects"),
**#58** ("the remaining review findings, and close the gap that hid them") and
the merge of #55 — 711 insertions this package had none of. Ported, with the
upstream tests, each placed where the seam now puts it rather than where it
came from.

- **Selection identifiers are quoted.** A blank chain left `chain  and resi 2`,
  where PyMOL takes `and` as the chain name and the selection stops being
  scoped to the object — so `qscore_view` on a validation file with no chain
  attribute ran `alter <every atom loaded>`, and `restore_bfactors` only
  restores the object it is given. A negative residue number is a *range*:
  `resi -3` means 1-3. Both were found upstream against PyMOL 3.1.0, and this
  package had reproduced both independently. Quoting lives in the PyMOL
  backend, not in `atoms`: it is PyMOL's grammar and no other viewer shares it.
  The `+`-list grouping survives for plainly numeric identifiers and anything
  else is quoted and OR-ed inside its chain's term, since a quoted value in a
  `+` list is grammar nobody here has checked.
- **A negative RMS is not a sigma scale.** MRC writes `rms=-1` for "statistics
  not computed", and it divides cleanly, so `to_sigma(0.05)` returned `-2.05`
  and a resolution ramp ran backwards. `to_absolute` had no guard at all.
- **The first B-factor stash wins, and a restore clears it.** Every caller
  reads `b` *after* an earlier view may have overwritten it, so a second stash
  saved that view's output as the user's data and `restore_bfactors` then wrote
  it back and reported success.
- **Provenance takes the longest matching token**, so `unsharpened` stops being
  read as `sharp`.
- **`loaded_map(obj, port)` evicts a record whose object has left the session**,
  and `density_view` names the file its header came from — the only thing that
  makes a same-name substitution visible.
- **`latent_traverse_view` anchors on the first frame with a usable RMS**, not
  frame 0.
- **`ensemble_spread_view` requires `superposed`.** Spread measures whatever
  separates the states, including a rigid-body offset that is not flexibility.
  Fitting stays with the host, which has the session and a superposition
  routine already — but the flag is a claim, so `internal_distance_change`
  checks it against the data. A rigid motion preserves every interatomic
  distance exactly, so a large positional spread with no internal rearrangement
  under it is a rigid offset whatever the caller said.

### Changed — source-breaking

Both hosts construct these types directly, so each of these needs a change on
their side. They are grouped here rather than buried in the fix list because
the changelog is the only place a consumer looks.

- **`Atom.index` is now `Atom.rank`**, and `ATOM_EXPR` asks PyMOL for `rank`.
  `index` is renumbered when atoms are removed — checked live on PyMOL 3.1.0,
  where `remove hydro` moved it and left `rank` alone — so a value stashed
  against it silently repoints at a different atom the moment anyone tidies a
  structure. `Atom.key` is `(model, rank)`; build per-atom `ScalarField` keys
  with `atom.key` rather than assembling the tuple.
- **`PymolBackend` and `draw` require a keyword-only `normalised`.** It had a
  default of `None` and `draw` never passed it, so a host that had correctly
  read `normalize_ccp4_maps off` still got a backend assuming the opposite.
  Pass what `normalisation_state` returned — including `None`, which now has to
  be chosen rather than fallen into.
- **`Frames` requires `numbers`**, the frame each surface was made from. It is
  not derivable from position: a frame whose header has no usable RMS is
  skipped, so surfaces run 1, 2, 4, 5 and numbering by position steps `frame 4`
  to frame 5's density.
- **`radial_spread` is gone**, replaced by `internal_distance_change`. It
  measured distance to the state's own centroid, which is blind to tangential
  motion: a counter-rotating twist — a real conformational change — scored
  exactly zero and was reported as rigid-body motion whatever the threshold.
  `RIGID_RATIO` moved 10.0 → 2.0 with it; the two are calibrated together and
  neither carries over to the other's quantity.
- **`rms_meaning` is public** (was `_rms_meaning`), because three modules print
  it and a second copy would be a second place for "rms=0 means flat" to be
  wrong.

### Notes

Between 2026-08-09 and this split the code lived inside MCPymol as
`mcpymol.wiggles` (PR #54, tiers 1–2; PR #55, tier 3). It moved out when a
second viewer — [protean](https://github.com/chemrich/protean), driving Mol\*
in a browser — made "a submodule of one viewer" the wrong shape. No released
MCPymol ever carried it: `v1.5.1` was cut before PR #54 landed, so nothing was
removed from anybody's install.
