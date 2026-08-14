# M5 — the three live-PyMOL unknowns, settled by observation (2026-08-13)

Run against Charlie's live PyMOL **3.1.0** on :9876, session empty beforehand,
scratch objects under `wfm5*` prefixes, all deleted afterwards.

FIX-PLAN.md listed three unknowns that "cannot be resolved by reasoning". All
three are now settled. **Two came back against what the plan assumed, and one
of them is an eleventh defect that the review did not find** — live in
`wiggles-em` *and* in MCPymol's `origin/main`, which is public on PyPI.

---

## Q1. Is `resi 1+2+3` valid while `resi "1"+"2"` is not?

**No — the premise is false. Both are valid.**

| Expression | Atoms matched |
|---|---|
| `resi 1+2` | 4 |
| `resi "1"+"2"` | 4 |
| `resi 1+"2"` (mixed) | 4 |
| `(resi "1" or resi "2")` | 4 |

`_PLAIN_RESI` in `backends/pymol.py:85` splits plain digits onto a `+` list and
everything else onto quoted `or` clauses, and its comment says a quoted value
in a `+` list "relies on grammar this package has not checked". It has now been
checked: the grammar accepts it. **The split is unnecessary complexity**, not a
bug — but see Q1b, because the branch it feeds is where the real defect lives.

## Q1b. NEW DEFECT — quoting does *not* escape a negative residue number

`quote()`'s docstring (`backends/pymol.py:108`) says a negative `resi` is read
as a range, "verified against PyMOL 3.1.0 upstream, **in both directions**" —
i.e. that quoting fixes it. **It does not.**

Structure: chain A residues **-3, 1, 2**, two atoms each, six total. So
"literal residue -3" is 2 atoms and "the range up to 3" is 6 — distinguishable.
The earlier probe used a structure with *no* residue -3, where both readings
return the same count, which is why this was never caught.

Run through the package's own `render_selection`:

| `Sel` | Lowered to | Atoms | Expected |
|---|---|---|---|
| `Sel.residues([("A","-3")])` | `(chain "A" and resi "-3")` | **6** | 2 |
| `Sel.residues([("A","-3"),("A","1")])` | `(chain "A" and (resi 1 or resi "-3"))` | **6** | 4 |
| `Sel.prop("resi","-3")` | `resi "-3"` | **6** | 2 |
| `Sel.residues([("A","1")])` | `(chain "A" and resi 1)` | 2 | 2 |

Three of four wrong. This is the **same data-destroying class as findings #1
and #2**: one residue's value written across three, `alter` reporting success.
Negative residue numbers are an expression-tag remnant in most NMR and EM
entries, so this is a common case, not an exotic one.

**The escape that does work is a backslash**, and it composes with a `+` list:

| Expression | Atoms |
|---|---|
| `chain "A" and resi \-3` | **2** ✅ |
| `chain "A" and resi \-3+1` | **4** ✅ |
| `chain "A" and resi "-3"` | 6 ❌ |

Because `\-3+1` works, the fix is *simpler* than the current code, not more
complex: escape the leading `-` and the `_PLAIN_RESI`/`odd` split can collapse
entirely.

**Upstream carries it too.** `git show origin/main:src/mcpymol/wiggles/atoms.py`
has the same `return f'"{value}"'` under the same claim ("Both were verified
against PyMOL 3.1.0, in both directions. Quoting fixes both"). MCPymol is public
on PyPI, so this needs reporting there regardless of what happens to the fork.

## Q2. Does `remove solvent` renumber `index` but not `rank`?

**The conclusion holds — but the review's specific repro does not.**

| Removal | `index` renumbered | `rank` renumbered |
|---|---|---|
| `remove hydro` | **yes** (6→4, 7→5, 8→6) | no |
| `remove resi 2` (a middle residue) | **yes** (5→3, 6→4) | no |
| `remove solvent` | **no** | no |

PyMOL sorts solvent to the *end* of index order, so removing it renumbers
nothing after it. Finding #2's stated repro — "the user then runs
`remove solvent`, which is routine before a figure" — **would not reproduce**.
`remove hydro`, which the finding also names, does.

This matters for M1: transcribing that repro literally would produce a test that
passes against the pre-fix code, and it would have been recorded as a
verification. That is exactly the failure class the plan exists to prevent.

`rank` was also confirmed **unique** per object, and it stays non-contiguous
after removal (0,2,4,5,6,7) — which is precisely what makes it stable.

**Decision 2 is therefore confirmed: `index` → `rank`.** Use `remove hydro` as
the repro.

## Q3. Does `alter` with a `stored` dict reach every atom of a quoted blank chain?

**Yes.** Against a structure with a real blank chain, a chain A and a chain B:

- `chain "" and resi 5` → 2 atoms, correctly scoped.
- `alter` with `b=stored.d.get('|'.join((model, str(rank))), b)` reached
  **2/2** blank-chain atoms and leaked onto **0** others.

So the quoting fix is sound for the blank-chain case — it is only the *negative
number* case above where quoting is the wrong escape.

---

## What this changes in FIX-PLAN.md

1. **A new Group A item**: the negative-`resi` escape. Same severity as #1/#2,
   and it also invalidates a docstring claim in two repos.
2. **Decision 2 confirmed**, with a corrected repro (`remove hydro`).
3. **Q1's premise is dead**: the `_PLAIN_RESI` split can collapse when the
   escape lands, since `\-3+1` parses.
4. **Report to MCPymol** — public, on PyPI, same bug, same false claim.

Scripts: `scratchpad/m5_probe.py`, `m5_q2.py`, `m5_q1q3.py`, `m5_negresi.py`.
