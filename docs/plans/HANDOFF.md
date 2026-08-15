# Handoff — 2026-08-14

Rewritten at a context boundary, replacing the 2026-08-13 version, which is now
wrong in three places: `wiggles-em` has a GitHub repo, phase 4 is not the next
phase, and the ten review findings are long since fixed.

**On rejoining, in one line:** the package is public and fixed, MCPymol is fixed
and keeps its own copy forever, and the live work is getting protean to consume
`wiggles-em` so the seam's central claim can finally be tested.

---

## Work in `~/code/wiggles-em`. There is one directory now.

Consolidated 2026-08-14, and **finished on 2026-08-16: there is one repository
now.** `~/code/wiggles-em` holds everything — the package, the plan, and the
unpublished material. **There are no symlinks left**, and `~/code/wiggles` is a
**frozen archive**: it still contains every file at the commit it was frozen at,
and it is not to be edited.

    ~/code/wiggles-em/
      src/, tests/          tracked, public
      docs/plans/           tracked, public — the plan and method documents,
                            including this one
      tools/                tracked, public — check_divergence.py
      private/              IGNORED WHOLESALE, no allow-list
        compendium/         the 17 entries, README, _TEMPLATE
        SPEC.md  PUBLICATION-GATE.md  REPOSITORIES.md  MOVING.md
        FLATIRON-REVIEW.md  FRAME-QUALITY.md
        tools/              check_gate.py, check_spec.py, livefire.py

The plan moved on 2026-08-14 because keeping it next door made it **stale**:
editing it meant reaching into another repo, so it did not get edited, and it
went on describing a closed gap as open. The compendium followed on 2026-08-16
for the same reason, once `private/` existed to receive it.

**The gate is unchanged.** Its ten open items are still the reason the
compendium is unpublished — its load-bearing findings rest on preprints, and
publishing them under this repo's name would turn each into a citable public
assertion. Living inside a public repository changes nothing, because nothing in
`private/` is tracked.

`check_gate.py` was updated to say so correctly. It used to ask `gh repo view`
whether the *repository* was public, which was right while the compendium was
the whole content of a private repo. That proxy inverts here — visibility says
PUBLIC while every gated file is unpublished — so it now asks git whether the
material is **tracked**. Untracked cannot have been pushed, whatever the repo's
visibility.

Both scripts resolve their root as `Path(__file__).parent.parent`, so the
mirrored layout under `private/` means they run unmodified:

    python3 private/tools/check_gate.py
    uv run python private/tools/check_spec.py

**Note what publishing the plan disclosed**, since it was a deliberate choice
and not an oversight: these documents describe protean in detail — its
op-support table, file paths and known defects — and protean is private, going
public, not yet flipped.

**`docs/plans` is deny-by-default in `.gitignore`, with the tracked documents
allow-listed back in.** Not the reverse. This directory is a drop point between
conversations, and a file left here by one of them carried its own assurance
that "wiggles-em gitignores this, so nothing here can reach the public repo by
accident" — true when written, false hours later, and staged by the next
`git add -A`. An allow-list means a new file stays private until somebody
deliberately adds a line.

**The symlinks are gone, and with them the rule that used to sit here.** Until
2026-08-16 the argument was *symlinks, not copies* — a copy would leave those
documents untracked in a public repo while the version-controlled originals sat
next door, which is the two-maintained-copies trap this project has already paid
for twice. What made that argument moot was freezing the other side: `wiggles`
is an archive, so there is exactly one live copy and nothing to diverge from.

**That safety is conditional, and the condition is the freeze.** Editing
`~/code/wiggles` re-creates the fork this arrangement was designed to avoid. Do
not.

The practical gain is that everything is now an ordinary file. Claude Code's
file tools refuse to write through a symlink, so editing `SPEC.md` used to mean
naming `~/code/wiggles/SPEC.md` explicitly — which is precisely why those
documents went stale. Nothing needs that any more.

**No trailing slashes** in `.gitignore` remains the standing rule even with no
symlinks left: a trailing slash matches a directory and not a symlink to one,
which is how protean came to commit `viewer/node_modules` as a mode-120000 link
past its own rule.

---

## Four repos

| Repo | State |
|---|---|
| `chemrich/wiggles` (this one, **private**) | Compendium, `SPEC.md`, the plan docs — all now tracked and pushed. The git home for everything symlinked above. Still gated by `PUBLICATION-GATE.md` |
| `chemrich/wiggles-em` (**PUBLIC**) | The working directory. `main` past `9381dbe` — a CI PR merged, plus the ignore rules. 524 tests, MIT. **Not on PyPI and not going there** |
| `chemrich/MCPymol` (public, PyPI) | `main` at `8357e02`. Two defects found here and fixed upstream this session (#59, #60). **Keeps its own copy of wiggles indefinitely** |
| `chemrich/protean` (**private**, going public — PR #66 merged) | Moves fastest of all. Consumes `wiggles-em` as a git pin at `9381dbe` (#67), and carries `backends/molstar.py` (#68) |

**The pin is now behind `wiggles-em`'s `main`.** protean resolves `9381dbe`
while main has moved past it. That is safe — a sha is a sha — but a wiggles-em
change does not reach protean until the pin is bumped.

---

## What to read, in order

1. **[`PROTEAN-INTEGRATION.md`](PROTEAN-INTEGRATION.md)** — the live plan. Five
   decisions with their implications, the op→action mapping, and seven steps.
2. **[`GUIDELINES.md`](GUIDELINES.md)** — how to work on this without repeating
   nine rounds of it. Every rule has a specific incident behind it.
3. **[`SPLIT.md`](SPLIT.md)** — the design and the amended phase plan.
4. `FIX-PLAN.md` — historical now, except its final section, which records that
   the exit criteria were never met.

---

## Decisions taken this session

- **PyPI is deferred indefinitely.** No date, no trigger.
- **Phase 4 (revert MCPymol's copy) is dead, not delayed.** PyPI rejects direct
  URL dependencies and MCPymol is on PyPI, so with no release there is no route
  to it. **Do not attempt `git revert 992edd7`** — it would delete a working,
  reviewed module and replace it with a dependency MCPymol cannot declare.
- **protean is the next phase and is unblocked**, because it is unpublished and
  can take a git dependency.
- **Preprints may admit a method**, not merely be cited — Charlie, 2026-08-13.
  Written into `compendium/README.md`'s Conventions, which is what governs
  entries. Preprint-admitted methods take `emerging` in the Status column.
- **Bug-hunting stopped after nine rounds**, on judgement. See the end of
  `FIX-PLAN.md`; criterion 6 was never met.

---

## The seam's central claim is still untested

Every view returns a `Scene` and never calls a viewer; `backends/` lowers it. That
is what the whole split rests on, and **only PyMOL has ever lowered one**.
`backends/fake.py` is a strictness check, not a second viewer.

The spike run this session says it should work. protean's viewer vocabulary,
read off its `_call` sites, covers every op the atom-based views emit:

    occupancy_view         ColorByScalar, Show, Legend          all available
    altloc_view            Hide, Show, ColorFlat, Label, Legend all available
    ensemble_spread_view   ColorByScalar, SizeByScalar, Legend  all available

Per-atom scalars ride the B-factor column of a re-sent display copy, read by
Mol\*'s `uncertainty` theme — which is exactly what `scene.py` predicted:
"they differ in destructiveness, not mechanism". PyMOL stashes originals
because it has one copy; protean re-sends a display copy and needs no stash.

**Missing:** `Isosurface` needs a volume action (the `cryoem-volumes` branch).
`Frames`, `Morph` and `Arrows` have no equivalent and must be **refused**. The
price of that: `Frames` is what `latent_traverse_view` needs, so protean cannot
do latent traversals until Mol\* grows a movie timeline.

---

## Traps that cost time here — do not rediscover

- **protean's CI checks out only protean**, so a `../wiggles-em` path dependency
  fails the moment a PR opens. The branch as pushed carries the path form and
  **must be converted to a git pin before PR**.
- **`/code-review <number>` resolves the number against the wrong repository.**
  It did so three times — once to wiggles-em, twice to protean. Pass a path or
  work inside the target repo.
- **The bash tool resets cwd between calls.** `cd` inside a compound command
  affects only that command; a `git log` after a `cd` reports the other repo.
- **Never work in `~/code/MCPymol` or `~/code/protean` directly.** Charlie has
  live sessions in both. Clone to scratch — a worktree writes to their `.git`.
- **`uv run pytest` and `python -m pytest` behave differently** without
  `pythonpath = ["src", "."]`. Both are green now; keep them both green.

---

## Open work, roughly in order

1. **protean integration** — `PROTEAN-INTEGRATION.md`. **Steps 1 and 2 are
   done**: the git pin merged as `a679690` (PR #67), and `backends/molstar.py`
   is PR #68. **A Scene has now been lowered by something that is not PyMOL**,
   which is what the whole split rested on and what had never been tested.

   Next is **step 3**, `occupancy_view` end to end, then looking at it in a
   browser — the step the plan exists to reach. Nothing calls the backend yet.

   **The op→action table in the plan was wrong in three rows** and is corrected
   there. `SizeByScalar` and `Label` cannot be honoured at all, and `Hide` only
   for `Rep.EVERYTHING`. Step 3 is unaffected; `ensemble_spread_view` is
   blocked on a protean-side size theme.
2. ~~**Re-verify the overdue preprints**~~ **DONE 2026-08-15** (wiggles
   `0a6b7fc`). Three of the five had been published: EMReady2 (Nat Commun
   17:6493), LocScale-2.0 (Nat Commun) and `cryospire` (NeurIPS 38). CryoLVM
   is at an unrecorded v2; `limits` is still v1 and still preprint-only, which
   is the case the rule exists for. CryoBench's `[unverified]` NeurIPS
   placement was confirmed at the same time.

   **Verified against Crossref, not the preprint servers** — an arXiv page
   shows a journal reference only if an author adds one, and arXiv showed none
   for `cryospire`, which has been in NeurIPS proceedings since 2025.

   **Left open:** `emerging` means *preprint-admitted*, and `cryospire` no
   longer is. Whether peer review alone promotes it is a decision about what
   the column means.
3. ~~**`check_spec.py`**~~ **DONE 2026-08-15** — reconciles against
   `wiggles_em.TOOLS`, with per-consumer coverage reported and never failed on.
   It found four drifts on its first run: tier 3 had shipped while `SPEC.md`
   still called `load_ensemble`, `latent_traverse_view`, `deformation_view` and
   `composition_view` unbuilt.
4. ~~**The divergence audit belongs in CI.**~~ **DONE 2026-08-15** —
   `tools/check_divergence.py`, run weekly by `.github/workflows/divergence.yml`
   and on demand. Four layers: modules, public names, signatures, and
   **behaviour**. Only the last would have caught either real incident, since
   both were behavioural changes to functions that existed with matching
   signatures on both sides — so it runs identical inputs through both copies
   and compares. Verified by reintroducing the historical provenance bug: it
   reports `postprocess_emready.mrc` as `nn_enhanced` upstream and `sharpened`
   here, which is exactly the incident.

   **Not on push or pull_request**, deliberately: the job installs another repo
   at its moving `main`, so wiring it into PR CI would train everyone to ignore
   a red tick for reasons unrelated to the commit. And it tracks `@main` rather
   than a pin, because an audit pinned to a frozen upstream reports agreement
   *while* the two drift — PyPI cannot serve it either, since the released
   `mcpymol` predates the subpackage.

   **Two divergences are recorded as OUTSTANDING**: `to_sigma`/`to_absolute`
   take a `MapStats` here and a `MapHeader` upstream, from this morning's
   change. They print under "Known divergences" every run rather than being
   exempted, and the success line says so instead of claiming agreement.
5. ~~**`9381dbe` is unreviewed**~~ **REVIEWED 2026-08-15** — the tenth round,
   against current `main` rather than the commit in isolation, since a defect
   fixed since is not worth reporting. The round-9 fixes hold: the putty
   duplicate still collapses to one note, and both `test_localres` remedy
   assertions still fail when the remedy is withheld or stops naming the
   literal command.

   **One real finding.** `9381dbe` reported K1 as mutation-checked, and on
   current `main` **both halves survived** — each rescues the other's test, so
   neither was guarded. Removing `_noted.clear()` as redundant leaves the suite
   green and makes a reused backend emit *zero* notes on a second render of the
   same view: B-factors overwritten, nothing said. Both halves now have their
   own guard, each mutation-verified.

   Also: `_note_once`'s docstring claimed "once per object" while keying on
   `(object, text)`, and `render_op` — a second public door that does not reset
   the ledger — had no docstring saying it is not an entry point.
6. Ten items in `PUBLICATION-GATE.md` still gate the compendium prose.

---

## What the nine rounds were actually about

Worth carrying, because it is the most transferable thing here and it is
written up properly in `GUIDELINES.md`:

- **Mutation testing proves necessity, not sufficiency.** Round 3 deliberately
  broke every fix, watched every test fail, and still shipped nine defects.
- **Corrections and designs need different defences.** "Transcribe the repro
  first" has nothing to transcribe when the fix is invented, and inventions
  produced the worst findings.
- **Stale prose adjacent to changed code** appeared in three consecutive rounds,
  each stale claim inviting a specific wrong action.
- **Fix the shape, not the symptom.** Three of round 7's four findings were a
  two-site defect fixed at one site.
- **Follow a remedy through and check the user is better off**, not merely that
  it runs. Two remedies pointed at each other; one produced a byte-identical
  refusal.
- **A test that passes on the negation of its claim is worse than no test.**
  This is what let the original bug ship in two repositories, and round 9 found
  it again in the tests written to verify round 8.
