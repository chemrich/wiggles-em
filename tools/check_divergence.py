#!/usr/bin/env python3
"""Check that MCPymol's copy of this code and this package have not diverged.

**Two maintained copies is the permanent steady state.** PyPI rejects direct
URL dependencies and MCPymol is on PyPI, so it cannot depend on this package
while this one stays off PyPI — which is the standing decision. MCPymol
therefore keeps `mcpymol.wiggles` indefinitely, and nothing has ever checked
the two agree.

That is not hypothetical. Two fixes have already been found here, fixed here,
and then found and fixed *again* upstream: the negative-`resi` selection escape
and the provenance category priority. Both were live in a shipped release in
between.

## What each layer can and cannot catch

Ordered by how much a wrong answer costs, which is the reverse of how easy the
layer is to write.

1. **Modules** — a module on one side with no counterpart. Cheapest, weakest.
2. **Public names** — a function added on one side and not the other. This is
   the audit that existed as a one-off, and it is genuinely useful: it catches
   a *new* tool that never got ported.
3. **Signatures** — a parameter added, removed or renamed on one side. Catches
   a class of drift names cannot see.
4. **Behaviour** — identical inputs through both, outputs compared. **This is
   the only layer that would have caught either real incident**, because both
   were behavioural changes to functions that existed, with matching
   signatures, on both sides. Layers 1-3 would have reported agreement
   throughout.

A probe that cannot run is reported as **UNRUNNABLE and fails the check**,
never silently skipped: "no failures" and "no checks" must not look alike.

## Why this tracks MCPymol's `main` rather than a pin

PyPI's `mcpymol` does not contain `mcpymol.wiggles` at all — the released
version predates the merge — so installing from PyPI cannot run this check.
A pinned git ref would be worse than useless: the audit exists to detect the
two copies drifting apart, and pinning one side means it reports agreement
*while* they drift. It has to follow `@main`.

    uv run python tools/check_divergence.py

Exit status is 0 when the copies agree, 1 when they do not, and 1 when the
comparison could not be made.
"""

from __future__ import annotations

import inspect
from types import ModuleType
from typing import Any

#: Modules expected on one side only, with the reason. Both entries are
#: structural: upstream owns the MCP tool boundary, this package owns the
#: viewer-neutral seam, and neither wants the other's.
MODULE_EXEMPT = {
    "tools": "upstream only — MCPymol's MCP tool boundary; this package registers nothing",
    "scene": "here only — the viewer-neutral seam; upstream calls PyMOL directly",
    "backends": "here only — Scene lowering, which upstream has no need of",
}

#: Public names upstream with no counterpart here, with the reason.
NAME_EXEMPT = {
    "residue_selection": "replaced here by Sel.residues + render_selection",
    "residue_clause": "replaced here by Sel.residues + render_selection",
    "quote": "here in backends/pymol.py — quoting is PyMOL selection syntax",
    "normalisation_state": "here in backends/pymol.py — reads a PyMOL setting",
}


#: Divergences that are known and accepted, with the reason. **Not the same as
#: agreement** — these are printed under their own heading every run so they
#: stay visible, and each one is either a deliberate consequence of the split
#: or an outstanding port. An entry whose reason begins "OUTSTANDING" is work
#: someone still owes.
ACCEPTED = {
    "atoms.Atom": (
        "model/rank are the Scene's atom key — `Atom.key` is `(model, str(rank))`, "
        "and upstream has no Scene to key against"
    ),
    "density.to_sigma": (
        "OUTSTANDING — `MapStats` landed here 2026-08-15 so a host can convert "
        "against measured statistics rather than a header's claim. Not yet ported "
        "upstream; MCPymol still converts against the header alone"
    ),
    "density.to_absolute": ("OUTSTANDING — see density.to_sigma"),
}


def _public(module: ModuleType) -> dict[str, Any]:
    """Public names *defined by* this module, not re-exported into it.

    Filtering on ``__module__`` matters: without it every ``from x import y``
    counts as a name the module owns, the two sides' import styles differ, and
    the audit reports divergence that is not there.
    """
    out = {}
    for name, obj in vars(module).items():
        if name.startswith("_"):
            continue
        if getattr(obj, "__module__", "").rsplit(".", 1)[-1] != module.__name__.rsplit(".", 1)[-1]:
            continue
        out[name] = obj
    return out


def _submodules(package: ModuleType) -> set[str]:
    import pkgutil

    return {m.name for m in pkgutil.iter_modules(package.__path__)}


def compare(upstream: ModuleType, here: ModuleType) -> tuple[list[str], list[str], list[str]]:
    """``(failures, notes, accepted)`` for the structural layers, 1 to 3."""
    import importlib

    failures: list[str] = []
    notes: list[str] = []
    accepted: list[str] = []
    architectural = 0

    up_mods, here_mods = _submodules(upstream), _submodules(here)
    for name in sorted(up_mods - here_mods):
        if name in MODULE_EXEMPT:
            notes.append(f"module `{name}` upstream only — {MODULE_EXEMPT[name]}")
        else:
            failures.append(
                f"module `mcpymol.wiggles.{name}` has no counterpart here. Port it, "
                f"or record it in MODULE_EXEMPT with a reason."
            )
    for name in sorted(here_mods - up_mods):
        if name in MODULE_EXEMPT:
            notes.append(f"module `{name}` here only — {MODULE_EXEMPT[name]}")
        else:
            notes.append(f"module `{name}` here only, unrecorded — upstream may be behind")

    for name in sorted(up_mods & here_mods):
        up = importlib.import_module(f"{upstream.__name__}.{name}")
        mine = importlib.import_module(f"{here.__name__}.{name}")
        up_names, my_names = _public(up), _public(mine)

        for missing in sorted(set(up_names) - set(my_names)):
            if missing in NAME_EXEMPT:
                notes.append(f"`{name}.{missing}` upstream only — {NAME_EXEMPT[missing]}")
            else:
                failures.append(
                    f"`mcpymol.wiggles.{name}.{missing}` has no counterpart in "
                    f"`wiggles_em.{name}`. A fix or feature landed upstream and not "
                    f"here — port it, or record it in NAME_EXEMPT."
                )

        for shared in sorted(set(up_names) & set(my_names)):
            a, b = up_names[shared], my_names[shared]
            if not (callable(a) and callable(b)) or inspect.isclass(a) != inspect.isclass(b):
                continue
            try:
                sig_a = list(inspect.signature(a).parameters)
                sig_b = list(inspect.signature(b).parameters)
            except (TypeError, ValueError):
                continue
            # The split itself is not drift. Upstream views take a `port` and
            # call PyMOL; here they take values and return a Scene. Comparing
            # those signatures compares two things designed to differ, and it
            # buried the one real finding under eleven false ones the first
            # time this ran.
            if "port" in sig_a and "port" not in sig_b:
                architectural += 1
                continue
            if sig_a != sig_b:
                key = f"{name}.{shared}"
                if key in ACCEPTED:
                    accepted.append(f"`{key}` — {ACCEPTED[key]}")
                    continue
                failures.append(
                    f"`{key}` takes {sig_a} upstream and {sig_b} here. "
                    f"A signature that drifts is a behaviour that drifted with it."
                )
    if architectural:
        notes.append(
            f"{architectural} signature(s) differ only by upstream's `port` argument — "
            f"the split itself, not drift. Counted rather than hidden."
        )
    return failures, notes, accepted


def probe_behaviour(upstream: ModuleType, here: ModuleType) -> tuple[list[str], list[str]]:
    """Layer 4 — identical inputs both sides, outputs compared.

    The only layer that would have caught either real incident. Add a probe
    whenever a rule is fixed in one copy, because that is exactly the moment
    the two are most likely to disagree.
    """
    failures: list[str] = []
    ran: list[str] = []

    # Provenance from filename tokens. This one *has* diverged: longest-token
    # matching dropped the category priority, and `postprocess_emready.mrc` —
    # the ordinary name for EMReady over a RELION postprocess map — was
    # reported SHARPENED instead of NN_ENHANCED. Live on PyPI until 2026-08-13.
    names = [
        "postprocess_emready.mrc",
        "run_postprocess_cryolvm.mrc",
        "postprocess_locscale.mrc",
        "unsharpened_map.mrc",
        "emd_30913_sharpened.mrc",
        "deepemhancer_output.mrc",
    ]
    try:
        import mcpymol.wiggles.provenance as up_p

        import wiggles_em.provenance as my_p

        header = my_p.MapHeader(**_MINIMAL_HEADER)
        up_header = up_p.MapHeader(**_MINIMAL_HEADER)
    except Exception as exc:
        failures.append(
            f"UNRUNNABLE: the provenance probe could not be set up ({exc!r}). A "
            f"probe that cannot run must not read as agreement."
        )
    else:
        for filename in names:
            try:
                a = up_p.gather_evidence(up_header, filename).suggested
                b = my_p.gather_evidence(header, filename).suggested
            except Exception as exc:
                failures.append(f"UNRUNNABLE: gather_evidence({filename!r}) raised {exc!r}")
                continue
            ran.append(filename)
            if getattr(a, "value", a) != getattr(b, "value", b):
                failures.append(
                    f"provenance of {filename!r}: upstream says "
                    f"{getattr(a, 'value', a)!r}, here says {getattr(b, 'value', b)!r}. "
                    f"This is the rule that has already diverged once."
                )
    return failures, ran


#: The smallest header both sides accept. Statistics are deliberately usable so
#: a probe fails on the rule under test rather than on an rms=-1 refusal.
_MINIMAL_HEADER = {
    "path": "probe.mrc",
    "byte_order": "<",
    "mode": 2,
    "magic": "MAP ",
    "nversion": 20140,
    "nx": 4,
    "ny": 4,
    "nz": 4,
    "mx": 4,
    "my": 4,
    "mz": 4,
    "nxstart": 0,
    "nystart": 0,
    "nzstart": 0,
    "cella": (4.0, 4.0, 4.0),
    "cellb": (90.0, 90.0, 90.0),
    "mapc": 1,
    "mapr": 2,
    "maps": 3,
    "dmin": 0.0,
    "dmax": 1.0,
    "dmean": 0.5,
    "rms": 0.25,
    "ispg": 1,
    "nsymbt": 0,
    "exttyp": "",
    "origin": (0.0, 0.0, 0.0),
    "nlabl": 0,
    "labels": (),
}


def main() -> int:
    try:
        import mcpymol.wiggles as upstream
    except ImportError as exc:
        print("NOT CHECKED: the two copies could not be compared.")
        print(f"    mcpymol.wiggles is not importable ({exc}).")
        print("    PyPI's mcpymol does not carry it; install from source:")
        print("      uv pip install 'mcpymol @ git+https://github.com/chemrich/MCPymol@main'")
        print()
        print("Exiting non-zero: an unrun audit must not read as agreement.")
        return 1

    import wiggles_em as here

    failures, notes, accepted = compare(upstream, here)
    behaviour_failures, probed = probe_behaviour(upstream, here)
    failures += behaviour_failures

    print("Divergence audit — mcpymol.wiggles vs wiggles_em")
    print(f"  modules upstream : {len(_submodules(upstream))}")
    print(f"  modules here     : {len(_submodules(here))}")
    print(f"  behaviour probes : {len(probed)} run")
    print()

    if notes:
        print(f"Recorded differences ({len(notes)})")
        for note in notes:
            print(f"  - {note}")
        print()

    if accepted:
        print(f"Known divergences ({len(accepted)}) — recorded, not resolved")
        for item in accepted:
            print(f"  ~ {item}")
        print()

    if failures:
        print(f"DIVERGENCE ({len(failures)})")
        for failure in failures:
            print(f"  ! {failure}")
        return 1

    owed = [a for a in accepted if "OUTSTANDING" in a]
    if owed:
        print(f"OK — nothing NEW has diverged. {len(owed)} known divergence(s) are still")
        print("owed, listed above; this is not the same as the copies agreeing.")
    else:
        print("OK. The two copies agree on modules, names, signatures and probed behaviour.")
    print("Only the *probed* behaviour is compared — see this module's docstring for")
    print("what each layer can and cannot catch.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
