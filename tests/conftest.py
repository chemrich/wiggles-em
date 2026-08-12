"""Shared test fixtures and the two ways a view gets exercised."""

from __future__ import annotations

import pytest

from wiggles_em.atoms import Atom
from wiggles_em.backends.fake import FakeBackend
from wiggles_em.backends.pymol import PymolBackend
from wiggles_em.bfactors import clear_stash
from wiggles_em.port import FakePort
from wiggles_em.provenance import forget


@pytest.fixture(autouse=True)
def _isolate_module_state():
    """Both the B-factor stash and the provenance registry are module-level.

    Without this, a view in one test leaves state that changes what another
    test sees — the sort of coupling that makes a suite pass in one order and
    fail in another. Provenance especially: a leaked declaration would make
    "UNKNOWN by default" pass for the wrong reason.
    """
    clear_stash()
    forget()
    yield
    clear_stash()
    forget()


def atom_rows(rows) -> list[tuple]:
    """Pad ``(chain, resi, resn, name, alt, q, b)`` rows out to ATOM_EXPR.

    Atom identity is ``(model, index)`` — chain/residue/name/altloc collide on
    insertion codes and across a selection spanning two models. Most tests do
    not care which model an atom is in, so they write seven fields and get a
    distinct index each; a test *about* identity writes all nine itself.
    """
    return [row if len(row) == 9 else (*row, "m", i) for i, row in enumerate(rows, start=1)]


def make_atoms(rows) -> list[Atom]:
    """Build atoms from seven- or nine-field rows."""
    return [
        Atom(chain=c, resi=i, resn=n, name=a, alt=alt, q=q, b=b, model=model, index=index)
        for c, i, n, a, alt, q, b, model, index in atom_rows(rows)
    ]


def iterate_response(rows) -> dict:
    """The port response a backend gets when it reads atoms back."""
    return {"iterate_to_list": atom_rows(rows)}


class Drawn:
    """The result of a view, plus what each backend made of it.

    Views return ``(report, Scene)`` and draw nothing themselves, so a test has
    two distinct things it can assert on and they answer different questions:

    - ``scene`` / ``ops`` — what the view *decided*. Invariants belong here.
      A claim about what is drawn is a claim about this value, and it cannot
      pass because of a substring that happened to match.
    - ``port`` — what PyMOL would actually receive. Regression tests belong
      here: they are what stops the seam from changing the picture.
    """

    def __init__(self, report: str, scene, port: FakePort, pymol: PymolBackend, fake: FakeBackend):
        self.report = report
        self.scene = scene
        self.port = port
        self.pymol = pymol
        self.fake = fake

    @property
    def ops(self):
        return list(self.scene)

    @property
    def notes(self) -> str:
        return "\n".join(self.pymol.notes)

    @property
    def full_report(self) -> str:
        """Report plus backend caveats — what a host would actually print."""
        return self.report + ("\n" + self.notes if self.notes else "")

    def ran(self, needle: str) -> bool:
        return self.port.ran(needle) or any(
            needle in str(a) or any(needle in str(x) for x in args)
            for a, args, _ in self.port.queries
        )


def render(result, rows=(), *, port=None, preserve_bfactors: bool = True, normalised=None) -> Drawn:
    """Render an already-computed ``(report, scene)`` through both backends.

    Rendering through *both* is deliberate. ``FakeBackend`` is strict about ops
    no viewer may honour, so a view that emitted a forbidden op fails here even
    if the PyMOL lowering happened to swallow it.
    """
    report, scene = result

    fake = FakeBackend()
    fake.render(scene)

    if port is None:
        port = FakePort(iterate_response(rows))
    pymol = PymolBackend(port, preserve_bfactors=preserve_bfactors, normalised=normalised)
    pymol.render(scene)

    return Drawn(report, scene, port, pymol, fake)


def draw(view, rows, *args, preserve_bfactors: bool = True, port=None, normalised=None, **kwargs) -> Drawn:
    """Run an atoms-first view over ``rows`` and render what it returned."""
    return render(
        view(make_atoms(rows), *args, **kwargs),
        rows,
        port=port,
        preserve_bfactors=preserve_bfactors,
        normalised=normalised,
    )
