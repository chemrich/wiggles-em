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


def make_atoms(rows) -> list[Atom]:
    """Build atoms from ``(chain, resi, resn, name, alt, q, b)`` rows."""
    return [Atom(chain=c, resi=i, resn=n, name=a, alt=alt, q=q, b=b) for c, i, n, a, alt, q, b in rows]


def iterate_response(rows) -> dict:
    """The port response a backend gets when it reads atoms back."""
    return {"iterate_to_list": list(rows)}


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


def draw(view, rows, *args, preserve_bfactors: bool = True, **kwargs) -> Drawn:
    """Run a view over ``rows`` and render its scene through both backends.

    Rendering through *both* is deliberate. ``FakeBackend`` is strict about ops
    no viewer may honour, so a view that emitted a forbidden op fails here even
    if the PyMOL lowering happened to swallow it.
    """
    atoms = make_atoms(rows)
    report, scene = view(atoms, *args, **kwargs)

    fake = FakeBackend()
    fake.render(scene)

    port = FakePort(iterate_response(rows))
    pymol = PymolBackend(port, preserve_bfactors=preserve_bfactors)
    pymol.render(scene)

    return Drawn(report, scene, port, pymol, fake)
