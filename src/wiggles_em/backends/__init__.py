"""Backends: the only code in this package that knows a viewer exists.

A backend takes a :class:`~wiggles_em.scene.Scene` and draws it. Everything
above here is viewer-neutral, so a backend is where the viewer-specific
knowledge is allowed to live — and where it is *required* to live.

Two rules every backend follows, both learned the hard way:

**Refuse rather than approximate.** An op a backend cannot honour raises
:class:`~wiggles_em.scene.Refused`. It is never skipped. A dropped op leaves a
picture that looks fine and means something other than what was asked for,
which is the failure this whole package exists to prevent.

**Normalisation is yours.** Levels and domains arrive in the units the data is
in. PyMOL contours in σ against a per-map normalisation, so
:class:`~wiggles_em.backends.pymol.PymolBackend` converts. A backend whose
viewer takes absolute levels passes them through. No view converts anything.
"""

from wiggles_em.backends.fake import FakeBackend
from wiggles_em.backends.pymol import PymolBackend

__all__ = ["FakeBackend", "PymolBackend"]
