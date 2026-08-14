# Working guidelines — written after three rounds of fixes for fixes, kept current

Findings per round: **10, 10, 9, 8, 6, 4, 3.** Rounds 4 onward ran under these rules.

The count declines slowly and the composition matters more than the total.
Round 4 introduced roughly one new defect across seven commits; round 5
introduced six across nine. **Round 5 was worse per commit, and it was the round
where I moved fastest** — nine commits in one sitting, several of them touching
message text I did not re-read.

Two things did improve, and both are worth watching:

- **Designs began surviving review.** Round 4's central change was a new
  statistic and a re-derived threshold; an independent reviewer exercised it
  across six motion types and found it sound. Round 3's designs were where all
  of its worst defects came from.
- **Severity fell.** Rounds 2 and 3 ended in data destruction and silently wrong
  scientific claims. Round 5's worst were a refusal message that made three
  false statements about the user's session, and a guard that crashed on the
  input it existed to reject. Both are visible to the user rather than silent.

One caveat on reading these numbers: each review's scope differed. Round 6
looked only at round 5's commits, so "all six were mine" is a scoping artifact,
not a trend. New defects per commit is the honest metric; total findings also
counts what earlier reviews left behind.

The problem was never carelessness — round 3 ran a deliberate-break mutation
check on every single fix and still produced nine. Something about the method
regenerates defects, and these rules are aimed at that rather than at effort.

Every rule below is derived from a specific thing that went wrong. Where a rule
has no incident behind it, it has been left out.

---

## Five philosophies

**1. A fix is a claim about the system, not about the symptom.**

Round 3 verified every fix twice: the symptom disappears, and the test fails
without the fix. Both checks are local. Neither asks what the change made the
author newly responsible for. `internal_distance_change` did exactly what it was
designed to do on the case the review named — and was compared against a
threshold calibrated for a different statistic, so it silently stopped firing on
ensembles with many states.

**2. Necessity is not sufficiency, and mutation testing only proves necessity.**

Deliberately breaking a fix and watching the test fail proves the fix is
load-bearing *for the reported symptom*. It says nothing about scope. Round 3
ran that check on all ten fixes and felt well defended by it. Treating a
necessity proof as a completeness proof is the single most expensive error in
this project's history.

**3. Corrections and designs need different defences, and only corrections
currently have one.**

M1 — transcribe the review's repro before reading the code — is a real defence
against sharing the bug's mental model. It works. But when a fix is a *design*
rather than a correction, there is no independent repro to transcribe and M1
provides no coverage at all.

Sort round 3's findings by severity: the HIGH was an invented statistic, the
worst MEDIUM was invented state, the next was a consequence of an invented
refusal. **All three of the worst were designs.** Anything invented needs an
adversary, a property test, and a reviewer who is not its author.

**4. The dangerous artifacts are the ones nothing forces you to look at.**

Report prose, docstrings, changelog entries, threshold constants. None of them
can fail a test, so they drift silently and then get quoted as evidence.
`quote()` claimed "verified against PyMOL 3.1.0, in both directions" for a fix
that had never been verified in either. `Granularity.ATOM` documented a key the
code had stopped using. `RIGID_RATIO` explained itself in terms of a function
that no longer existed. Two reports hard-coded `rms=0` for a map whose rms was
`-1`.

**5. Anything that can hold a belief twice will eventually hold two different
ones.**

Defaults, duplicated derivations, prose restating a contract. `normalised`
defaulted to `None` and `draw()` never passed it, so the view and the backend
disagreed. `Frames` carried names but not numbers, so the backend re-derived the
numbering and got it wrong. The fix in both cases was to delete the second
holder, not to synchronise the two.

---

## Choosing what goes in one change

"Smaller PRs" is the right instinct but not actionable on its own. Here is the
rule, derived by sorting round 3's nine findings by what kind of change produced
them.

**A change ships alone if it does any of these:**

1. **Introduces persistent state.** (`_DESTROYED` — new state has a lifecycle,
   and the lifecycle was only walked along the paths already being edited.)
2. **Adds or removes a required parameter or field on a shared type.**
   (`Frames.numbers`, required `normalised` — the blast radius is every
   construction site, which drowns anything else in the same diff.)
3. **Changes what a number or key *means*.** (`index` → `rank`, range → RMS —
   every consumer has a calibration you are silently inheriting.)
4. **Converts a refusal into a survival, or the reverse.** (`localres`
   surviving `rms=-1` moved the failure downstream into a layer nobody looked
   at.)
5. **Rewrites a block rather than editing it.** (`_frames` — the rewrite
   dropped a guard that existed in the code it replaced.)

**Everything else — genuine corrections — may be batched**, but not across more
than one module.

Applied to round 3, this turns seven commits into roughly nine, nearly all of
them alone. That is the concrete answer to "how do we actually make them
smaller": not a line budget, but a categorical rule about which changes are too
entangling to share a diff.

**Land each one before starting the next.** Round 3 finished #3, started #5, and
then could not cleanly separate them — a self-inflicted loss of bisectability
that was noticed at commit time and too late to undo.

---

## Before writing the fix: the obligations ledger

Write down, in the commit message draft, **what this change makes you newly
responsible for** — not what it fixes. Every round-3 defect is an unlisted entry
in this ledger:

| If the change… | …you now owe |
|---|---|
| adds state | who creates it, who clears it, and what *external* event invalidates it |
| changes a quantity | every consumer, **including thresholds calibrated on the old one** |
| adds a required argument | every construction site, and a decision at each |
| loosens a refusal | the downstream that must now cope with the survivor |
| rewrites a block | the invariants in the replaced code, named explicitly |
| edits report text | every sibling saying the same thing elsewhere |

**If the ledger has entries in more than one row, it is more than one change.**

---

## While writing it

**When fixing a logic error, grep for the shape — not the symptom.** The
reported instance is usually one of several. Round 7 found four defects and
three were the same meta-error: a two-site defect fixed at one site.

- H2 fixed `sorted()` over raw caller keys in the arity guard and left the
  identical bug in the duplicate-key guard **three lines below**, in the same
  method. A later audit found a *third* instance in `Frames`.
- H1 made a refusal's remedy takeable on the defaulted path; three of the four
  situations that reach that refusal still offered only the unusable one.
- G6 made a held stash outrank the destroyed mark in `stash_bfactors`; the
  backend's `_stash` had the same ordering and still reported the opposite.

The ledger already says "edits report text → find every sibling saying the same
thing". This is the same instruction for logic: the mirror branch, the sibling
call site, the guard three lines down. Ask *what shape is this bug* before
asking whether it is fixed — and start the fix from that audit, not from the
line the review named.

**Grep twice: for the predicate and for the string.** Round 3's M2 grep was
`header.rms`, which found the code and missed both hard-coded `"rms=0"` strings
in report prose. Predicates and their prose drift apart independently.

**Property tests wherever a threshold or a statistic lives.** Sweep the
parameter the calibration depends on. The HIGH finding is invisible on a clean
fixture and instantly fatal under a sweep of state count × noise. A noiseless
test of a statistical quantity tests nothing statistical.

**The deletion test: would this test still pass if I deleted the thing I just
built?** If yes, it is incomplete. Round 3's `radial_spread` replacement was
first covered only by "a twist is not flagged" — which passes just as well if
the detector is removed entirely. The positive case got added because the author
happened to think of it, not because anything required it.

**Execute every remedy a message names.** A message that tells the user what to
do is a promise the code makes, and promises are testable. Two of them in one
session were not kept: `destroyed_note` said "reload the object" when nothing
on any load path cleared the mark, and `localres`' refusal said "give the level
with `units='absolute'`" when that path raised before the refusal was even
reached. **The second was written immediately after fixing the first**, which is
what makes this a rule rather than an anecdote — knowing about the failure mode
did not prevent repeating it two commits later.

Name the argument rather than describing it, so the advice is copy-pasteable and
a test can run exactly what the user is told to run.

And **follow the remedy through — check the user ends up better off**, not
merely that the call succeeds. Round 8 found two remedies that satisfied
"execute it in a test" and were still wrong: one led to a second refusal telling
the user to undo the first (the two pointed at each other), and one produced a
byte-for-byte identical refusal, because the condition it addressed was not the
one that had failed. A remedy that costs a reload and changes nothing is worse
than no remedy.

And **execute it from the state the user is actually in, not one you build.**
The test written for this rule passed only because it invented a contour value
of its own. On the path that mattered — no level given — the user has no value
to supply, because computing it is the thing that just failed. Supplying what a
remedy asks for has already assumed the user can.

**When a question is about reality, get an observation.** No amount of reasoning
settles whether PyMOL renumbers `index`, or whether quoting escapes a negative
`resi`. M5 answered three such questions live in under an hour, **corrected the
review** on one of them (`remove solvent` renumbers nothing), and found an
eleventh defect that had shipped to PyPI under a false "verified" claim. This is
the highest-yield hour in the whole project; do it earlier and more often.

**Retuning a constant is a stop signal.** `RIGID_RATIO` went 3.0 → 10.0 after a
false positive. That retuning was treating a symptom of a quantity that was
blind to a whole class of motion; no threshold could have worked. When you find
yourself adjusting a number because of a bad result, the null hypothesis is that
the thing being measured is wrong.

**A behaviour change includes its docstring in the diff.** Not as context — as
part of the change. When you edit a line you read the line, and the paragraph
around it keeps its old promises: after G3 the comment still said `_TOKENS`
order decided classification; after G6 the docstring still promised "Returns 0";
after G7 the attribute was still documented "by surface name" while both writes
keyed it by volume. Four of round 6's six findings were prose adjacent to code
I had just changed, and this is the third consecutive round that class has
appeared.

The stale prose is worse than absent prose in every one of those cases, because
each invited a specific wrong action: reorder `_TOKENS`, probe destruction with
`== 0`, look up a key never written. If the enclosing docstring is not in the
hunk, say why not.

**Never reformat what you did not change.** `ruff format` touched 11 files, 9
unrelated to the change; unrelated churn hides the real diff from the reviewer
you are depending on.

---

## Before landing

**When your test fails against your code, you have no information about which is
wrong.** Both encode the same model. Round 3 wrote **four** test fixtures that
were wrong before the code was: an assertion using `"original"` as a proxy for a
claim, a fixture reading the wrong port channel, a fixture modelling `rank` as
renumbering, and a "counter-rotating" pair that rotated both atoms the same way
(a rigid rotation, correctly scoring zero). Resolve these against an external
referent — a live session, the review's own repro text, physical reasoning —
never against your own reading of the code.

**A check that cries wolf loses its licence.** Two of M3's three consistency
checks fired on correct code in their first draft: one demanded every number in
a report be an emitted level (reports legitimately state one level both ways),
the other read "Skipped: frame 3" as an instruction to type `frame 3`. Narrow a
check until it fires only on a real contradiction, and record the false
positives so nobody re-derives them.

**Run the suite over degenerate input, not only well-formed input.** Every
divergence in the round-2 review appeared exactly when something was *skipped* —
a frame dropped, a header without statistics. M3's class check caught nothing
until a gapped ensemble was added, and then caught finding #5 independently.

**"The suite is green" is not evidence.** It was true at the end of all three
failed rounds. 359 tests, ruff and mypy clean, and four deliberate-break checks
were all true of a package with ten confirmed defects in it.

---

## Who reviews

**Do not grade your own round when the round contains designs.** The author's
verification uses the author's model, which is the thing under suspicion. The
review that found the round-2 defects worked because it was told to distrust
"covered by a test" claims — that instruction turned every verification into a
mutation test, cost 83 minutes instead of 15, and produced five findings that
say "Reproduced:".

Point the next review hardest at whatever the last round *designed*, and tell it
which changes those were.

---

## What worked and should be kept

Not everything here is a failure mode. These earned their place:

- **Removing a default rather than documenting it.** Every time a default was
  deleted, call sites were forced to decide, and several of them decided
  differently from the default — which is the bug surfacing. It works.
- **Collapsing two paths into one.** The `_PLAIN_RESI` split existed because a
  grammar claim had never been checked; checking it removed a whole branch and
  with it a permanent drift risk. Prefer changes that *reduce* the number of
  seams over changes that fix behaviour at one.
- **Live observation over documentation.** See M5.
- **Testing a failure class rather than its instances.** M3 caught #5 on its
  own, from a different direction than the test written for #5.
- **Recording what a mechanism cost and what it caught**, so the next person can
  judge whether to keep paying for it.
