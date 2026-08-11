"""A backend that records instead of drawing.

Mostly a convenience: a :class:`~wiggles_em.scene.Scene` is already a value, so
a test can assert on what a view *returned* without rendering it at all, which
is the point of the seam. What this adds is the ability to ask "would a real
backend have been asked to do something it cannot do" without a viewer, and to
answer questions about a scene after it has been flattened.

It honours the refusal rule: :attr:`strict` makes it raise on ops that no
viewer could satisfy, so a test cannot pass because the fake was more
permissive than everything real.
"""

from __future__ import annotations

from wiggles_em.scene import Legend, Refused, ScalarField, Scatter, Scene, SceneOp


class FakeBackend:
    """Records every op it is given.

    Ops land in :attr:`ops` in order. With ``strict=True`` (the default), a
    :class:`~wiggles_em.scene.Scatter` raises — it is defined so I2 can name
    what it forbids, and a backend that quietly accepted one would let a view
    emit the thing the invariant exists to prevent.
    """

    def __init__(self, strict: bool = True) -> None:
        self.ops: list[SceneOp] = []
        self.strict = strict

    def render(self, scene: Scene) -> None:
        for op in scene:
            self.render_op(op)

    def render_op(self, op: SceneOp) -> None:
        if self.strict and isinstance(op, Scatter):
            raise Refused("Scatter: no view may emit a latent scatter (I2)")
        self.ops.append(op)

    # -- test conveniences -------------------------------------------------

    def of(self, *types: type) -> list[SceneOp]:
        """Every recorded op of the given types, in order."""
        return [op for op in self.ops if isinstance(op, types)]

    def has(self, *types: type) -> bool:
        """Was any op of these types recorded?"""
        return any(isinstance(op, types) for op in self.ops)

    @property
    def legends(self) -> list[Legend]:
        return [op for op in self.ops if isinstance(op, Legend)]

    @property
    def fields(self) -> list[ScalarField]:
        """Every scalar field drawn, for checking domains and lengths."""
        return [op.field for op in self.ops if hasattr(op, "field")]

    @property
    def transcript(self) -> str:
        """One op per line — for assertion messages."""
        return "\n".join(repr(op) for op in self.ops)


__all__ = ["FakeBackend"]
