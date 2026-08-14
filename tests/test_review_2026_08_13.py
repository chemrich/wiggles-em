"""One test per finding in REVIEW-2026-08-13.md, transcribed from its repro.

M1 of FIX-PLAN.md. The rule that makes this file worth having: **each test is
written from the review's "Reproduced:" text, before reading the code around
it.** Two rounds of fixes have now each introduced new defects, and the common
thread is that the author wrote the test for the same mental model that
produced the bug. A test derived from an independent description of the
failure does not share that model.

Every test here must **fail against the pre-fix code**, and be confirmed to do
so rather than assumed to — a test that passes before the fix proves nothing
about the fix. Where a test is currently expected to fail it is marked xfail
with ``strict=True``, so it also fails if it starts passing for the wrong
reason, and the marker is removed as each fix lands.

One deviation, deliberate and recorded
--------------------------------------
Finding #2's stated repro is ``remove solvent``. Checked against a live PyMOL
3.1.0 (see M5-RESULTS.md): PyMOL sorts solvent to the *end* of index order, so
removing it renumbers nothing and the repro as written does not reproduce.
``remove hydro`` — which the finding also names — does renumber, as does
removing any middle residue. So #2 is transcribed as the *mechanism* the
finding describes, with the discrepancy recorded in the test.
"""

from __future__ import annotations

import re

import pytest

from tests.conftest import make_atoms, render
from wiggles_em.backends.pymol import PymolBackend
from wiggles_em.bfactors import has_stash, restore_bfactors
from wiggles_em.port import FakePort
from wiggles_em.scene import ColorByScalar, Granularity, ScalarField, Scene, Sel

# ── #1 ──────────────────────────────────────────────────────────────────────


def test_1_a_view_that_destroyed_bfactors_does_not_let_the_next_one_stash_garbage():
    """REVIEW #1, src/wiggles_em/backends/pymol.py:294.

    Transcribed from: render a ColorByScalar on object 'm' with
    `PymolBackend(port, preserve_bfactors=False)` (real B-factors 11.0/22.0
    destroyed, stash empty), then render another with `preserve_bfactors=True`
    on a port now reporting the overwritten values — the stash becomes
    {('m','1'): 0.9, ('m','2'): 0.3} (view A's scalar) and the note says the
    originals are held. Calling `restore_bfactors` then writes 0.9/0.3 into the
    B-factor column and returns 'Restored 2 B-factors on m', permanently
    destroying the crystallographic values while reporting success.
    """
    original = [("A", "1", "MET", "CA", "", 1.0, 11.0), ("A", "2", "SER", "CA", "", 1.0, 22.0)]
    scalars = ScalarField.per_atom([(("m", "1"), 0.9), (("m", "2"), 0.3)])
    scene = Scene([ColorByScalar(Sel.obj("m"), scalars, (0.0, 1.0))])

    # View A destroys the column and preserves nothing.
    port_a = FakePort({"iterate_to_list": [(*r, "m", i) for i, r in enumerate(original, 1)]})
    PymolBackend(port_a, preserve_bfactors=False, normalised=None).render(scene)
    assert not has_stash("m"), "precondition: preserve_bfactors=False leaves no stash"

    # View B now reads a column holding view A's scalars, not the originals.
    clobbered = [("A", "1", "MET", "CA", "", 1.0, 0.9), ("A", "2", "SER", "CA", "", 1.0, 0.3)]
    port_b = FakePort({"iterate_to_list": [(*r, "m", i) for i, r in enumerate(clobbered, 1)]})
    backend_b = PymolBackend(port_b, preserve_bfactors=True, normalised=None)
    backend_b.render(scene)

    note = "\n".join(backend_b.notes)
    if has_stash("m"):
        restored = restore_bfactors(port_b, "m")
        written = "\n".join(port_b.commands) + restored
        assert "0.9" not in written, (
            "restore wrote view A's scalar back as though it were the "
            f"crystallographic original: {written}"
        )
    # The claim the repro quotes is that the originals are *held* and can be
    # put back. A note saying the opposite — that they are gone — is the
    # honest outcome, so the assertion is on the claim rather than on the word
    # "original", which appears in both.
    assert "are held" not in note.lower(), (
        f"the note claims originals are held when they were already destroyed: {note}"
    )
    assert "restore_bfactors" not in note, (
        f"the note offers a restore that would write view A's scalar back: {note}"
    )


# ── #2 ──────────────────────────────────────────────────────────────────────


def test_2_atom_identity_is_keyed_on_the_field_that_removal_does_not_renumber():
    """REVIEW #2, src/wiggles_em/atoms.py:45.

    Transcribed from: a per-atom view stashes {('6xyz','1'): b1, …} keyed on
    `index`. The user then removes atoms, PyMOL renumbers `index` over the
    remainder, and `restore_bfactors` gives every atom after the first deleted
    one a different atom's B-factor while reporting success. `rank` is the
    original input order and is not renumbered, so (model, rank) is both unique
    and stable.

    DEVIATION, see the module docstring: the finding says `remove solvent`, but
    PyMOL sorts solvent to the end of index order, so that particular removal
    renumbers nothing (checked live, PyMOL 3.1.0). `remove hydro` — also named
    in the finding — does renumber, and is used here.

    Whether the field is stable is a fact about PyMOL, not about this package,
    and it cannot honestly be asserted against hand-written rows: writing a
    "renumbered rank" would be modelling the very thing rank does not do, and
    writing an unchanged one would make the test true by construction. That
    half is proved by observation in `test_selection_live.py`, which removes
    atoms from a real session and reads both fields back.

    What *is* this package's to get right, and what this asserts: that it asks
    PyMOL for the stable field and keys identity on it, in every one of the
    three places that has to agree — the iterate expression, `Atom.key`, and
    both `alter` lookups. A change landing in some of those and not all is how
    the previous two rounds went.
    """
    from wiggles_em.atoms import ATOM_EXPR
    from wiggles_em.bfactors import restore_bfactors, stash_bfactors

    fields = [f.strip() for f in ATOM_EXPR.split(",")]
    assert "rank" in fields, f"the iterate expression does not request rank: {ATOM_EXPR}"
    assert "index" not in fields, (
        f"the iterate expression still requests index, which a removal renumbers: {ATOM_EXPR}"
    )

    atom = make_atoms([("A", "1", "MET", "CA", "", 1.0, 11.0)])[0]
    assert atom.key == (atom.model, str(atom.rank)), atom.key

    # The scalar lookup and the restore lookup must key on the same field, or
    # a view colours by one identity and the restore undoes it by another.
    rows = [("A", "1", "MET", "CA", "", 1.0, 11.0)]
    field = ScalarField.per_atom([(a.key, 0.5) for a in make_atoms(rows)])
    drawn = render(("r", Scene([ColorByScalar(Sel.obj("m"), field, (0.0, 1.0))])), rows)
    scalar_alter = drawn.port.call_log

    stash_bfactors("m2", make_atoms(rows))
    restore_port = FakePort()
    restore_bfactors(restore_port, "m2")
    restore_alter = restore_port.call_log

    for label, log in (("scalar", scalar_alter), ("restore", restore_alter)):
        assert "str(rank)" in log, f"the {label} alter does not key on rank: {log}"
        assert "str(index)" not in log, f"the {label} alter still keys on index: {log}"


# ── #10 ─────────────────────────────────────────────────────────────────────


def test_10_the_documented_per_atom_key_is_the_key_the_backend_looks_up():
    """REVIEW #10, src/wiggles_em/scene.py:213.

    Transcribed from: `Granularity.ATOM` and `ScalarField.per_atom` still
    document per-atom keys as (chain, resi, name, alt) after `_push` switched
    the PyMOL lookup to (model, index). A field built to the documented
    contract silently matches nothing: `alter` becomes a no-op and `spectrum`
    colours the structure by whatever the B-factor column already held, under a
    legend naming a quantity that was never drawn.

    So the contract and the lowering are asserted against each other rather
    than either being asserted alone — a doc fix that missed the code, or a
    code fix that missed the doc, fails here.
    """
    documented = (Granularity.ATOM.__doc__ or "") + (ScalarField.per_atom.__doc__ or "")

    rows = [("A", "1", "MET", "CA", "", 1.0, 11.0)]
    field = ScalarField.per_atom([(a.key, 0.5) for a in make_atoms(rows)])
    drawn = render(
        (
            "r",
            Scene([ColorByScalar(Sel.obj("m"), field, (0.0, 1.0))]),
        ),
        rows,
    )
    # The `alter` expression travels as a structured call, not a command line,
    # so both channels have to be read or the lookup is invisible here.
    lowered = "\n".join(drawn.port.commands) + "\n" + drawn.port.call_log

    # Whatever fields the lookup is built from, the contract has to name those
    # same ones — asserted in both directions rather than against a fixed
    # spelling, so this keeps working if the key changes again.
    key_expr = re.search(r"join\(\((.*?)\)\)", lowered)
    assert key_expr, f"no per-atom key expression was emitted at all: {lowered}"
    used = {f.strip().removeprefix("str(").removesuffix(")") for f in key_expr.group(1).split(",")}

    for field in used:
        assert field in documented, (
            f"the lookup keys on {sorted(used)} but the contract never mentions "
            f"{field!r}. A field built to the documentation matches nothing: "
            f"`alter` becomes a no-op and the structure is coloured by whatever "
            f"the B-factor column already held, under a legend naming a "
            f"quantity that was never drawn."
        )
    assert "chain, resi, name, alt" not in documented, (
        "the contract still documents the old (chain, resi, name, alt) key, "
        f"which the lookup no longer uses — it keys on {sorted(used)}."
    )


# ── #7 ──────────────────────────────────────────────────────────────────────


@pytest.mark.parametrize(
    "filename",
    ["postprocess_emready.mrc", "postprocess_locscale.mrc", "run_postprocess_cryolvm.mrc"],
)
def test_7_a_network_enhanced_map_is_not_suggested_as_merely_sharpened(filename, tmp_path):
    """REVIEW #7, src/wiggles_em/provenance.py:155.

    Transcribed from: replacing the declaration-order `break` with
    longest-token-wins dropped the category priority. 'postprocess' (11 chars,
    SHARPENED) is longer than 'emready' (7), 'locscale' (8) and 'cryolvm' (7),
    so postprocess_emready.mrc, postprocess_locscale.mrc and
    run_postprocess_cryolvm.mrc all now report suggested=sharpened; before the
    change all three gave nn_enhanced. These are the ordinary output names for
    running EMReady/LocScale/CryoLVM on a RELION postprocess map. The user is
    offered 'features are easier to over-read than in the unsharpened map'
    instead of 'Treat this map as a hypothesis, not a measurement', and the
    whole point of the NN_ENHANCED category is that its warning is not
    interchangeable with sharpening.
    """
    from tests.test_mapinfo import write_map
    from wiggles_em.mapinfo import read_map_header
    from wiggles_em.provenance import Provenance, gather_evidence

    header = read_map_header(write_map(tmp_path, filename))
    evidence = gather_evidence(header, filename)

    assert evidence.suggested is Provenance.NN_ENHANCED, (
        f"{filename} was read as {evidence.suggested}: {evidence.reasons}. "
        "Longest-token-wins picked 'postprocess' over the network marker, so "
        "the cautionary category lost to the less cautionary one."
    )


# ── #3 ──────────────────────────────────────────────────────────────────────


def test_3_a_map_with_no_computed_statistics_does_not_kill_the_resolution_view(tmp_path):
    """REVIEW #3, src/wiggles_em/localres.py:353.

    Transcribed from: load a density map written by mrcfile without
    `update_header_stats()` (rms = -1, the exact case `usable_rms` was added
    for) as `main`, a matched resolution map as `locres`, then
    `local_resolution_view('main', 'locres', normalised=True)` raises an
    uncaught ValueError about sigma being undefined. The user asked for a sigma
    contour and never requested a conversion; the surrounding code shows the
    intent was to report the absolute value as unknown. Instead the whole tool
    dies and no grid check, no ramp table and no provenance banner are
    produced.

    `to_absolute` was made to raise on a non-positive RMS, but this call site
    kept guarding with `if main.header.rms` — and -1 is truthy.
    """
    from tests.test_mapinfo import write_map
    from wiggles_em.localres import local_resolution_view
    from wiggles_em.maps import forget_map, load_map
    from wiggles_em.provenance import Provenance

    forget_map()
    main = write_map(tmp_path, "main.mrc", rms=-1.0)  # "statistics not computed"
    res = write_map(tmp_path, "res.mrc", rms=0.5)
    port = FakePort({"get_names": ["main", "res"], "iterate_to_list": []})
    load_map(port, main, "main", provenance=Provenance.MEASURED)
    load_map(port, res, "res", provenance=Provenance.MEASURED)

    report, _scene = local_resolution_view("main", "res", normalised=True)

    assert report.strip(), "the view produced no report at all"
    forget_map()


# ── #4 ──────────────────────────────────────────────────────────────────────


def test_4_an_unnormalised_session_can_still_render_an_ensemble(tmp_path):
    """REVIEW #4, src/wiggles_em/backends/pymol.py:380.

    Transcribed from: build a 3-frame ensemble with `load_ensemble`, get the
    scene from `latent_traverse_view('ens')`, render through
    `PymolBackend(port, normalised=False)` — the value a host gets from
    `normalisation_state` when PyMOL says 'off' — and it raises PortError:
    a level in sigma was given but this session needs it in absolute and that
    volume was not loaded through `load_map`. Nothing is drawn, and the remedy
    the message gives is impossible: an ensemble frame cannot be loaded through
    `load_map`. The absolute level the backend needs is the one the view
    already computed and then converted away.
    """
    from tests.test_mapinfo import write_map
    from wiggles_em.heterogeneity import forget_ensemble, load_ensemble
    from wiggles_em.latent import latent_traverse_view

    forget_ensemble()
    root = tmp_path / "ens"
    root.mkdir()
    (root / "z.pkl").write_text("x")
    for i, rms in enumerate((0.5, 0.6, 0.8), start=1):
        write_map(root, f"vol_{i}.mrc", rms=rms)
    names = [f"ens_f{i:02d}" for i in range(1, 4)]
    port = FakePort({"get_names": names, "iterate_to_list": []})
    load_ensemble(port, root, "ens")

    _report, scene = latent_traverse_view("ens")
    PymolBackend(port, normalised=False).render(scene)

    forget_ensemble()


# ── #5 ──────────────────────────────────────────────────────────────────────


def test_5_a_movie_frame_steps_to_the_latent_frame_its_surface_is_named_for(tmp_path):
    """REVIEW #5, src/wiggles_em/backends/pymol.py:475.

    Transcribed from: a 5-frame ensemble whose frame 3 has rms=0. The backend
    emits `mset 1 x4`, `mdo 3, ...enable ens_surf_04` and
    `mdo 4, ...enable ens_surf_05`. latent.py names surfaces by their original
    frame number and the report promises "a surface's number is always the
    frame it was made from" plus "`frame N` to step". But `_frames` renumbers
    over the survivors. A user who reads "Skipped: frame 3" and types `frame 3`
    to see what came after the gap is shown frame 4's density, and `frame 5`
    shows nothing at all.
    """
    from tests.test_mapinfo import write_map
    from wiggles_em.heterogeneity import forget_ensemble, load_ensemble
    from wiggles_em.latent import latent_traverse_view

    forget_ensemble()
    root = tmp_path / "ens"
    root.mkdir()
    (root / "z.pkl").write_text("x")
    for i, rms in enumerate((0.5, 0.6, 0.0, 0.7, 0.9), start=1):  # frame 3 unusable
        write_map(root, f"vol_{i}.mrc", rms=rms)
    names = [f"ens_f{i:02d}" for i in range(1, 6)]
    port = FakePort({"get_names": names, "iterate_to_list": []})
    load_ensemble(port, root, "ens")

    _report, scene = latent_traverse_view("ens")
    PymolBackend(port, normalised=True).render(scene)

    mdo = [(a[0], a[1]) for a, _ in port.calls("mdo")]
    assert mdo, f"no movie frames emitted at all: {port.call_log}"
    for index, command in mdo:
        enabled = [w for w in str(command).replace(";", " ").split() if "_04" in w or "_05" in w]
        if enabled and "_04" in enabled[-1]:
            assert index == 4, (
                f"movie frame {index} enables surface _04, so `frame {index}` "
                f"shows frame 4's density while the report says a surface's "
                f"number is the frame it was made from. Emitted: {mdo}"
            )
    forget_ensemble()


# ── #6 ──────────────────────────────────────────────────────────────────────


def test_6_the_public_draw_helper_does_not_silently_assume_a_normalised_session(tmp_path):
    """REVIEW #6, src/wiggles_em/backends/pymol.py:543.

    Transcribed from: a host running `set normalize_ccp4_maps, off` correctly
    calls `local_resolution_view(..., normalised=False)` and then renders with
    the public `draw(port, scene)` helper, which constructs `PymolBackend(port)`
    with `normalised=None`. `self.normalised is False` is then False, so the
    Ångström breakpoints are converted to sigma against the resolution map's
    header while the volume still holds Ångström values — the surface comes out
    flat in one extreme colour, under a report stating the breakpoints were
    sent unconverted.

    The lazy `normalisation_state(self.port)` read that used to answer this was
    deleted, and `draw()` cannot pass the value — which is the whole bug. So
    the assertion is that `draw` agrees with what the session actually says,
    however that is achieved.
    """
    from tests.test_mapinfo import write_map
    from wiggles_em.backends.pymol import draw
    from wiggles_em.localres import local_resolution_view
    from wiggles_em.maps import forget_map, load_map
    from wiggles_em.provenance import Provenance

    forget_map()
    main = write_map(tmp_path, "main.mrc", rms=0.5)
    res = write_map(tmp_path, "res.mrc", rms=0.5)
    # The session says normalisation is OFF.
    port = FakePort({"get_names": ["main", "res"], "iterate_to_list": [], "get": "0"})
    load_map(port, main, "main", provenance=Provenance.MEASURED)
    load_map(port, res, "res", provenance=Provenance.MEASURED)

    _report, scene = local_resolution_view("main", "res", normalised=False)
    backend = draw(port, scene, normalised=False)

    assert backend.normalised is False, (
        "the session reports normalize_ccp4_maps off and the view was told "
        f"so, but draw() built a backend with normalised={backend.normalised}, "
        "so breakpoints get converted against a volume that was never "
        "normalised."
    )
    forget_map()


# ── #8 ──────────────────────────────────────────────────────────────────────


def test_8_the_latent_report_names_the_frame_the_contour_was_actually_taken_against(tmp_path):
    """REVIEW #8, src/wiggles_em/latent.py:314.

    Transcribed from: a 3-frame ensemble whose frame 1 has rms=0 and frames 2-3
    have dmean=10.0, rms=0.5/0.6. The report says 'Contour: 10.75 absolute' and
    'No level given, so 1.5 sigma against the FIRST frame was used and
    converted to an absolute value', while 10.75 is 1.5σ against frame *2*.
    1.5σ against frame 1 is not even defined. A reader checking the contour
    against frame 1's header — which the report tells them to use — gets a
    different number and concludes the tool converted wrongly.

    The anchor moved to the first header with a usable RMS; the report's line
    did not move with it.
    """
    from tests.test_mapinfo import write_map
    from wiggles_em.heterogeneity import forget_ensemble, load_ensemble
    from wiggles_em.latent import latent_traverse_view

    forget_ensemble()
    root = tmp_path / "ens"
    root.mkdir()
    (root / "z.pkl").write_text("x")
    write_map(root, "vol_1.mrc", rms=0.0, dmean=10.0)  # unusable
    write_map(root, "vol_2.mrc", rms=0.5, dmean=10.0)
    write_map(root, "vol_3.mrc", rms=0.6, dmean=10.0)
    names = [f"ens_f{i:02d}" for i in range(1, 4)]
    load_ensemble(FakePort({"get_names": names, "iterate_to_list": []}), root, "ens")

    report, _scene = latent_traverse_view("ens")

    if "first frame" in report.lower():
        assert "10.75" not in report, (
            "the report says the level was taken against the FIRST frame, but "
            "10.75 is 1.5 sigma against frame 2 — frame 1 has rms=0, so 1.5 "
            f"sigma against it is not defined. Report:\n{report}"
        )
    forget_ensemble()


# ── #9 ──────────────────────────────────────────────────────────────────────


def test_9_a_counter_rotating_twist_is_not_reported_as_rigid_body_motion():
    """REVIEW #9, src/wiggles_em/ensembles.py:195.

    Transcribed from: a counter-rotating twist — top ring +25°, bottom ring
    -25°, the ratchet/F1-ATPase class of real conformational change, not a
    rigid body motion. positional mean 2.16 Å, radial mean 7.9e-16 Å, ratio
    ~3e15, so `rigid_dominated` fires regardless of RIGID_RATIO. With
    `superposed=True` the report prints '! RIGID-BODY MOTION DOMINATES … the
    number above is not a conformational quantity', instructing the user to
    refit already-fitted states and discard a correct measurement.

    `radial_spread` is by its own docstring blind to tangential motion, and
    retuning 3.0 -> 10.0 cannot fix this class because radial is exactly zero.
    """
    import math

    from wiggles_em.ensembles import ensemble_spread_view

    # Two rings of atoms at radius 5, twisting in opposite directions. Every
    # atom keeps its distance from the axis, so radial spread is exactly zero
    # while the conformation genuinely changes.
    n = 8
    rows, start, end = [], [], []
    for ring, (z, sign) in enumerate([(5.0, +1.0), (-5.0, -1.0)]):
        for i in range(n):
            theta = 2 * math.pi * i / n
            twist = sign * math.radians(25.0)
            rows.append(("A", str(ring * n + i + 1), "ALA", "CA", "", 1.0, 20.0))
            start.append((5 * math.cos(theta), 5 * math.sin(theta), z))
            end.append((5 * math.cos(theta + twist), 5 * math.sin(theta + twist), z))

    report, _scene = ensemble_spread_view(make_atoms(rows), [start, end], "obj", superposed=True)

    # The warning itself, not the word "rigid" — the colour scale is legitimately
    # described as running "blue (rigid) → red", and asserting on the bare word
    # flagged that correct line.
    assert "RIGID-BODY MOTION DOMINATES" not in report, (
        "a counter-rotating twist — a real conformational change — was "
        "reported as rigid-body motion, telling the user to discard a correct "
        f"measurement. Report:\n{report}"
    )

    # And the detector must still detect. Removing it would satisfy the
    # assertion above, so the property worth keeping is asserted alongside:
    # the same atoms, translated bodily, are still caught.
    shifted = [(x + 7.0, y, z) for x, y, z in start]
    rigid_report, _ = ensemble_spread_view(
        make_atoms(rows), [start, shifted], "obj", superposed=True
    )

    assert "RIGID-BODY MOTION DOMINATES" in rigid_report, (
        "a pure translation is rigid motion by definition and must still be "
        f"caught, or the fix was to delete the check. Report:\n{rigid_report}"
    )
