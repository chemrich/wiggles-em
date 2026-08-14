# Fix plan — wiggles-em, after the 2026-08-13 review

Ten confirmed defects in [`REVIEW-2026-08-13.md`](REVIEW-2026-08-13.md). Read
that first; this file is how to fix them without producing a third round.

**The premise of this plan:** fix-then-trust-the-suite has now failed twice, in
the same way both times. Round 1 refactored and got ten findings. Round 2 fixed
those ten and got ten more, **most of them defects in the fixes**. So the plan
adds mechanisms aimed at the failure *classes*, and treats "the suite is green"
as meaning nothing on its own.

---

## Decisions needed before writing code

**1. `radial_spread` — remove it, or replace the invariant?**

It cannot be tuned into correctness. For a counter-rotating twist (the
F1-ATPase class of genuine conformational change) it is *exactly zero*, so the
detector fires on a correct measurement whatever `RIGID_RATIO` is. Options:

- **Remove it.** The report states whether the host fitted, and stops there. The
  original finding — rigid drift misread as flexibility — goes back to being
  the host's problem, which is where upstream left it with `intra_fit`.
- **Replace the invariant.** What is actually rigid-invariant is the *internal
  distance matrix*. If every pairwise distance is unchanged across states, the
  motion is entirely rigid; that is exact rather than heuristic, and it catches
  counter-rotation. Cost is O(N²), so it needs sampling for large N — and a
  sampled check is a heuristic again, though a much better founded one.
- **Keep the flag, drop the check.** Weakest; leaves a claim nobody verifies.

Recommendation: **replace the invariant**, with a bounded random sample of atom
pairs and a stated sample size in the report. It is the only option that keeps
the property worth having — a wrong `superposed` claim gets caught — without
resting on a quantity that is blind to a whole class of real motion.

**2. `index` → `rank` for atom identity.** The review is right and I was wrong:
`index` is renumbered by `remove solvent`, `rank` is the original input order
and is not. This changes `ATOM_EXPR` a second time and touches every test row
again. Confirm before doing it, because it is the second churn of the same
surface in two days.

**3. The strategic one: is the fork still the right shape?**

MCPymol's `mcpymol.wiggles` has been through two review rounds and is working.
This fork has been through two and is not. An honest alternative:

- **Abandon this fork.** Re-extract from MCPymol `origin/main` *with no
  refactor at all*, publish that, and then introduce the Scene seam in small
  reviewed steps against a package that is already known-good.

The case for continuing instead: the Scene work is sound in design — none of
the twenty findings across both rounds were "the seam is wrong", they were
lowering bugs, partial applications and report divergences. Re-extracting
throws away work that the reviews did not fault.

Recommendation: **continue, but fix in the order below and re-review.** Revisit
if round three finds another ten.

---

## Why the last two rounds failed

Five classes, all visible in the current ten. The mechanisms below map onto these.

| Class | Example from this round |
|---|---|
| **Partial application** — a change made at one call site, not all | `usable_rms` added in `density`, left as truthiness in `localres` and the backend |
| **Two fixes contradicting** | latent names surfaces by original frame number; `_frames` renumbers over survivors |
| **Report/render divergence** | report names frame 1 as the anchor; the code anchors on the first *usable* frame |
| **Wrong primitive chosen** | `index` (renumbered) instead of `rank` (stable); `radial_spread` instead of an internal-distance invariant |
| **Docs left behind the code** | `Granularity.ATOM` still documents `(chain, resi, name, alt)` after the lookup moved to `(model, index)` |

Every one of them passed 359 tests, ruff, mypy, and four deliberate-break checks.

---

## Mechanisms to add — do these before the fixes

**M1. Tests come from the review's repro text, not from my understanding.**
Every finding in `REVIEW-2026-08-13.md` carries a concrete "Reproduced:"
scenario. Transcribe it into a test *before* reading the code around it. This
breaks the loop where the author writes a test for the same mental model that
produced the bug — which is how ten fixes passed a suite the same session.

**M2. A completeness grep per changed predicate, recorded in the commit.**
When a predicate or key changes, grep for every occurrence of the *old* form
and either fix or justify each. Concretely, for this round:

```
rg 'header\.rms' --type py          # must all be usable_rms() now
rg 'chain, resi, name, alt'         # docs and code, after index/rank
rg 'normalisation_state|normalised' # view, backend, draw(), tests
rg 'enumerate\(.*names'             # any numbering over a filtered list
```

The grep output goes in the commit message. "I checked" is not checkable; a
recorded grep is.

**M3. A report-vs-render consistency suite.** Four of the ten are the report
asserting something the scene or backend does not do. That is testable as a
class rather than case by case: for each view, assert that factual claims in
the report are borne out by the emitted ops — a named frame number exists as a
surface, a stated unit matches the emitted `Unit`, a named object appears in a
selection. New `tests/test_report_consistency.py`.

**M4. Fix interacting pairs in one commit, not separately.** The latent
numbering and `_frames` numbering are one bug in two files; so are the three
`normalised` sites. Fixing them apart is what let them disagree.

**M5. Settle the three live-PyMOL unknowns with `tools/livefire.py`.** These
cannot be resolved by reasoning and are currently assumptions:

- Is `resi 1+2+3` valid while `resi "1"+"2"` is not? (The `+`-list split rests
  on this.)
- Does `remove solvent` renumber `index` but not `rank`? (Decision 2 rests on
  this.)
- Does `alter` with a `stored` dict lookup reach every atom in a selection
  built from a quoted blank chain?

**Ask before running** — the live suite clears the session, and Charlie's
PyMOL usually has work in it.

**M6. Close the two gaps upstream's #58 identified**, both still open here:

- *Defaults are never exercised.* Every test passes explicit arguments, so no
  documented default is ever called. Add a test that calls each view the way
  its signature says it may.
- *`FakePort` answers `"OK"` to anything it does not recognise*, so "the viewer
  accepted the signature and rejected the argument's meaning" passes silently.
  Make unknown commands raise, with an allowlist for the ones genuinely
  fire-and-forget.

---

## The fixes, grouped and ordered

Order is by blast radius, and pairs that interact are one unit.

**Group A — data destruction. Do first.**

1. **#1 `preserve_bfactors=False` leaves no stash.** The next view then stashes
   the clobbered column and `restore_bfactors` writes it back as "the
   originals". Likely fix: record that an object's B-factors are *known
   destroyed*, so a later stash refuses rather than saving garbage.
2. **#2 `index` → `rank`** (pending decision 2). Second churn of `ATOM_EXPR`;
   do it in one commit with the docs (#10) since they are the same change.

**Group B — partial application. One commit each, each with its M2 grep.**

3. **#3 `localres` truthiness → `usable_rms`**, plus the same guard surviving
   in the backend.
4. **#10 `Granularity.ATOM` / `ScalarField.per_atom` docs** — fold into #2.

**Group C — contradictions, fixed as pairs.**

5. **#5 + latent numbering.** `_frames` must number by the surface's own frame
   number, not by position. The report already promises this.
6. **#6 `normalised` across view, backend and `draw()`** — the helper cannot
   pass it, which is the whole bug. Either `draw()` reads the session itself,
   or it takes the value, or it goes away.

**Group D — semantics regressed by a fix.**

7. **#7 provenance priority.** Longest-token-wins dropped
   most-cautionary-first. Needs *both*: category priority (NN_ENHANCED before
   SHARPENED) **and** longest-token within a category.
8. **#8 latent report names the wrong anchor** — report and docstring both.

**Group E — hard failure.**

9. **#4 latent frames have no `load_map` record**, so an unnormalised session
   cannot render an ensemble at all. The absolute level is already computed in
   the view and converted away; carry it instead of recomputing.

**Group F — unsound design (decision 1).**

10. **#9 `radial_spread`.**

---

## Re-review, and what to tell it

Same settings that worked: **high effort, and the instruction to distrust
"covered by a test" claims** — that is what turned verification into mutation
testing, cost 83 minutes instead of 15, and produced five findings that say
"Reproduced:". Do not economise on it.

Add three things this round:

- **"These commits are fixes for a prior review. The last two rounds of fixes
  each introduced new defects, mostly by applying a change at some call sites
  and not all. Check completeness of every changed predicate."**
- **"For each fix, check that its test would actually fail against the
  pre-fix code."** M1 makes this true by construction; the review should verify
  it rather than take it.
- **"Check that report text and emitted ops agree"** — the class M3 covers, in
  case M3's own coverage is thin.

---

## Exit criteria

Not "the suite is green". That has been true at the end of both failed rounds.

- Every one of the ten has a test transcribed from the review's repro, and that
  test fails against the current code before the fix lands.
- Every changed predicate has a recorded grep showing no old form survives.
- `tests/test_report_consistency.py` exists and covers all fourteen views.
- The three live-PyMOL unknowns are settled by observation, or explicitly
  marked as assumptions in the code that rests on them.
- `FakePort` raises on unknown commands; defaults are exercised.
- A re-review at high effort returns **no confirmed correctness findings**.

Only then is the package worth publishing — and publishing is what unblocks
MCPymol depending on it, so this is the real critical path, not Phase 7.


---

## Outcome, recorded 2026-08-14 — criterion 6 was NEVER met

Nine review rounds ran. **Not one came back clean.**

    round   1    2    3    4    5    6    7    8    9
    found  10   10    9    8    6    4    3    5   (fixed, unreviewed)

The work stopped on **judgement, not on the bar being cleared.** Anyone reading
"nine review rounds" as "this was verified" is reading it wrong, and that
misreading is exactly the failure this project kept finding in its own
docstrings — `quote()` said "verified in both directions" and was not, #58's
changelog said "impossible to reintroduce" and it was not.

**Five of the six criteria were met**: transcribed tests that failed pre-fix,
recorded greps, the report-consistency suite, the three live-PyMOL unknowns
settled by observation, and `FakePort` raising on unknown commands with the
defaults exercised.

**What justified stopping**, stated so it can be argued with:

- Severity fell and stayed down. Rounds 2–3 ended in data destruction and
  silently wrong scientific claims. **Nothing since round 4 has produced a
  wrong number or a wrong picture.** Rounds 7–9 found prose, error paths, and
  loose test assertions.
- The character of findings changed. Round 9 was one code defect reachable only
  by a host reusing a backend, and four defects in the *tests* verifying the
  previous round. Review had begun finding the verification looser than the
  code.
- The package has no consumer, so message defects cost nothing today.
- The next validation is a different kind: hands-on use through protean, which
  tests what nine rounds never touched — whether a `Scene` lowers onto a
  non-PyMOL viewer at all.

**What that leaves open.** `9381dbe` (the K-series) is unreviewed, and every
round so far found something in the previous round's fixes — nine for nine. The
honest expectation is that a tenth round would find one or two more things of
the same character. That is a reason to stop, not a reason to believe there is
nothing there.
