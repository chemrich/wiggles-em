"""How much of each state there is — and, just as importantly, how that was worked out.

State occupancy is recoverable from cryo-EM. That is the premise of this
module, and it is a correction to what this package used to say.

**The thing that is unreliable is one particular way of computing it.** Counting
which particle images fall into which 3D class, or histogramming them in a
latent space, commits the *base-rate fallacy*: deciding which conformation an
image belongs to requires already knowing the population distribution, and doing
it without that knowledge biases the answer. Under the noise levels of real
cryo-EM the bias is severe — Evans *et al.* show latent-space histogramming
return a nearly flat distribution for a genuinely trimodal one, missing the
middle mode entirely.

**Methods that solve the inverse problem instead get it right.** Deconvolving
the latent space (RECOVAR) recovers all three modes; ensemble reweighting
against the image stack recovers the correct distribution outright. Both use the
statistics of the whole dataset rather than committing to a per-image
assignment.

So a population weight is only as good as its provenance, and two weights that
render identically as a bar height can differ in whether they mean anything at
all. That is the same problem :mod:`wiggles_em.provenance` exists for, and it
gets the same answer: **the source travels with the number, it is declared and
never inferred, and the default is UNKNOWN.**

One caution that makes labelling non-optional rather than nice to have. Evans
*et al.* note that at low noise, counting *is* reliable — but that "determining
the exact noise regime of real data is not possible". So the number can never
tell you whether it was safe to compute that way. Only the label can.

**Refs.**

- Evans L, Dingeldein L, Covino R, Gilles MA, Thiede E, Cossio P. *Counting
  particles could give wrong probabilities in cryo-electron microscopy.*
  bioRxiv 2025.03.27.644168. **Preprint.**
- Tang WS *et al.* *Ensemble reweighting using cryo-EM particle images.*
  *J Phys Chem B* — the reweighting route.
- Tang WS, Soules J, Rangan A, Cossio P. *CryoLike.* *Acta Cryst D*,
  doi:10.1107/s2059798325009350 — the likelihood engine, peer-reviewed.
- Gilles MA, Singer A. *PNAS* 122(9):e2419140122 (2025) — RECOVAR, the
  deconvolution route. See the compendium entry `recovar`.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from enum import Enum

#: Gas constant in kJ/(mol·K). Free energies here are kJ/mol, matching the
#: convention used elsewhere in this project's tooling.
GAS_CONSTANT_KJ = 8.314462618e-3

#: Weights below this are treated as "not resolvably populated" rather than as
#: a probability. A state at 1e-9 is not a measurement of a very rare state; it
#: is a number that fell out of a normalisation, and taking its logarithm
#: produces a confident-looking free energy from nothing.
NEGLIGIBLE_WEIGHT = 1e-6


class WeightSource(Enum):
    """How a set of population weights was computed.

    The categories are read off Evans *et al.* Figure 1, which benchmarks
    exactly these estimators against a known answer. They are not a taxonomy
    invented here.
    """

    THREE_D_CLASSIFICATION = "3d_classification"
    """Particle counts per 3D class. Biased under noise even with known poses."""

    HARD_ASSIGNMENT = "hard_assignment"
    """Each image assigned to its single most likely conformation, then counted."""

    SOFT_ASSIGNMENT = "soft_assignment"
    """Each image's likelihood spread across conformations, then summed. Still
    per-image, so still subject to the same bias."""

    LATENT_HISTOGRAM = "latent_histogram"
    """Counts of embedded images per latent-space bin. The commonest way to read
    a population off a heterogeneity method, and the one Evans *et al.* show
    returning a flat distribution for a trimodal truth."""

    DECONVOLVED = "deconvolved"
    """The latent density deconvolved against the noise model — RECOVAR's
    route. Recovers modes that histogramming loses."""

    ENSEMBLE_REWEIGHTED = "ensemble_reweighted"
    """Weights fitted to the whole image stack by image-to-structure
    likelihood, e.g. via CryoLike. Recovers the distribution directly."""

    UNKNOWN = "unknown"
    """Not declared and not inferable. The honest default, and what any weight
    arriving without a stated method gets."""

    @property
    def is_per_image_assignment(self) -> bool:
        """Does this work by deciding what each image is, then tallying?

        The property that makes an estimator vulnerable, rather than the name of
        any one tool. A future method that assigns per image belongs here too.
        """
        return self in (
            WeightSource.THREE_D_CLASSIFICATION,
            WeightSource.HARD_ASSIGNMENT,
            WeightSource.SOFT_ASSIGNMENT,
            WeightSource.LATENT_HISTOGRAM,
        )

    @property
    def solves_inverse_problem(self) -> bool:
        """Does this fit the distribution to the whole dataset at once?

        The property that makes an estimator trustworthy. ``UNKNOWN`` is false
        here *and* false for :attr:`is_per_image_assignment` — not knowing is
        its own category, not a middle ground.
        """
        return self in (WeightSource.DECONVOLVED, WeightSource.ENSEMBLE_REWEIGHTED)

    @property
    def caveat(self) -> str:
        """The one-line warning that belongs in any readout of these weights."""
        return {
            WeightSource.THREE_D_CLASSIFICATION: (
                "Counted per 3D class. Biased under noise even when poses are "
                "known, and real data's noise regime cannot be determined, so "
                "the size of the bias here is unknown."
            ),
            WeightSource.HARD_ASSIGNMENT: (
                "Each image assigned to one conformation, then counted. Assigning "
                "an image requires already knowing the distribution being "
                "measured, so this is biased by construction under noise."
            ),
            WeightSource.SOFT_ASSIGNMENT: (
                "Per-image likelihoods summed across conformations. Softening the "
                "assignment does not remove the bias — it is still decided one "
                "image at a time."
            ),
            WeightSource.LATENT_HISTOGRAM: (
                "Counted per latent-space bin. This is the estimator shown to "
                "return a nearly flat distribution for a genuinely trimodal one, "
                "losing the sparsely populated middle state entirely."
            ),
            WeightSource.DECONVOLVED: (
                "Latent density deconvolved against the noise model, using the "
                "whole dataset. Recovers modes that histogramming loses."
            ),
            WeightSource.ENSEMBLE_REWEIGHTED: (
                "Fitted to the full image stack by image-to-structure likelihood. "
                "The most defensible route to a population currently available."
            ),
            WeightSource.UNKNOWN: (
                "NOT DECLARED. How these weights were computed decides whether "
                "they mean anything, and it was not stated. Treat them as "
                "unusable for any quantitative claim until the method is known."
            ),
        }[self]


@dataclass(frozen=True)
class Populations:
    """Normalised state occupancies, with the method that produced them.

    Weights are stored **normalised**, because everything downstream — a bar
    height, a free energy — is a statement about probability, and a set of
    weights that does not sum to one silently changes what those mean.
    :attr:`raw_total` keeps what they summed to before, since a total far from
    one is usually a sign the caller passed counts, or a subset.

    Construct through :meth:`declare`. There is no inference anywhere in this
    class: nothing looks at a filename, a method name or the shape of the
    distribution to guess how it was made.
    """

    probabilities: tuple[float, ...]
    source: WeightSource
    raw_total: float
    uncertainty: tuple[float, ...] | None = None
    """One standard deviation on each probability, on the same normalised
    scale. Optional — most methods do not report it — and its absence is stated
    rather than filled in with a default."""

    temperature_k: float | None = None
    """The temperature the ensemble was equilibrated at **before** vitrification.

    Deliberately not defaulted, and deliberately not the temperature of the
    grid. Populations in a cryo-EM dataset are trapped from the equilibrium that
    existed in solution before plunge-freezing, so a Boltzmann reading of them
    refers to *that* temperature. Filling in 298 K, or worse the cryogenic
    temperature, would put a plausible number on an assumption the user never
    made.
    """

    @classmethod
    def declare(
        cls,
        weights: list[float] | tuple[float, ...],
        source: WeightSource = WeightSource.UNKNOWN,
        *,
        uncertainty: list[float] | tuple[float, ...] | None = None,
        temperature_k: float | None = None,
    ) -> Populations:
        """Validate and normalise a set of weights.

        Raises:
            ValueError: empty, negative, non-finite, summing to zero, or an
                uncertainty whose length does not match.
        """
        values = tuple(float(w) for w in weights)
        if not values:
            raise ValueError("no weights given; a population of nothing is not a population")
        if any(not math.isfinite(w) for w in values):
            raise ValueError("weights contain a NaN or an infinity")
        if any(w < 0 for w in values):
            negatives = [i for i, w in enumerate(values, start=1) if w < 0]
            raise ValueError(
                f"weights {negatives} are negative. A negative occupancy is not a "
                f"small one — it usually means a fit was run without a "
                f"non-negativity constraint, and normalising it would hide that."
            )
        total = math.fsum(values)
        if total <= 0:
            raise ValueError("weights sum to zero, so they cannot be normalised")

        if uncertainty is not None:
            errors = tuple(float(u) for u in uncertainty)
            if len(errors) != len(values):
                raise ValueError(
                    f"{len(values)} weights but {len(errors)} uncertainties — these "
                    f"must be parallel, or the error bars land on the wrong states"
                )
            if any(u < 0 or not math.isfinite(u) for u in errors):
                raise ValueError("uncertainties must be finite and non-negative")
            scaled: tuple[float, ...] | None = tuple(u / total for u in errors)
        else:
            scaled = None

        if temperature_k is not None and (not math.isfinite(temperature_k) or temperature_k <= 0):
            raise ValueError("temperature must be a positive number of kelvin")

        return cls(
            probabilities=tuple(w / total for w in values),
            source=source,
            raw_total=total,
            uncertainty=scaled,
            temperature_k=temperature_k,
        )

    @property
    def n_states(self) -> int:
        return len(self.probabilities)

    @property
    def is_quantitative(self) -> bool:
        """May these weights carry a numeric claim?

        True only for estimators that fit the whole dataset. Everything else —
        including ``UNKNOWN`` — is drawable but must not be read as a
        measurement.
        """
        return self.source.solves_inverse_problem

    def banner(self) -> str:
        """The block every population readout carries."""
        lines = [
            f"  Population source: {self.source.value.upper()}",
            f"  {self.source.caveat}",
        ]
        if not self.is_quantitative:
            lines.append(
                "  NOT QUANTITATIVE: read these as an ordering, not as numbers. "
                "Recomputing them by deconvolution or ensemble reweighting is what "
                "would make them numbers."
            )
        if self.uncertainty is None:
            lines.append(
                "  No uncertainty was supplied, so none is shown. That is not the "
                "same as the weights being precise."
            )
        if abs(self.raw_total - 1.0) > 0.01:
            lines.append(
                f"  Normalised from a total of {self.raw_total:g}; the values shown "
                f"are fractions of that total."
            )
        return "\n".join(lines)

    def relative_free_energy(
        self, reference: int = 0
    ) -> list[tuple[float, float | None, float | None]]:
        """ΔG of each state relative to ``reference``, in kJ/mol.

        ``ΔG_i = -RT ln(p_i / p_ref)``. Only *differences* are recoverable from
        occupancies, which is why there is a reference state and no absolute
        value anywhere.

        Returns one ``(delta_g, minus, plus)`` per state, where ``minus`` and
        ``plus`` are the distances to the low and high ends of the error bar.
        Both are ``None`` when no uncertainty was supplied.

        **The error bars are asymmetric on purpose, and this is the whole reason
        temperature and uncertainty are modelled at all.** A logarithm turns a
        symmetric error on a probability into a lopsided one on a free energy,
        and the effect explodes as the probability approaches its own error bar.
        A state at p = 0.02 ± 0.01 has a well-defined lower bound and an
        *unbounded* upper one — it might be arbitrarily unfavourable. Those are
        exactly the sparsely populated states this whole approach exists to
        find, so drawing them with symmetric error bars would misrepresent
        precisely the states that matter most. ``plus`` is ``None`` in that
        case, meaning unbounded, rather than a large finite number.

        Raises:
            ValueError: no temperature was declared, the reference index is out
                of range, or the reference state is not resolvably populated.
        """
        if self.temperature_k is None:
            raise ValueError(
                "no temperature was declared, so these occupancies cannot be turned "
                "into free energies. It is the temperature the ensemble was "
                "equilibrated at before vitrification, not the temperature of the "
                "grid, and this package will not guess it: pass temperature_k= to "
                "Populations.declare()."
            )
        if not 0 <= reference < self.n_states:
            raise ValueError(f"reference state {reference} is outside 0..{self.n_states - 1}")
        p_ref = self.probabilities[reference]
        if p_ref < NEGLIGIBLE_WEIGHT:
            raise ValueError(
                f"state {reference} has probability {p_ref:g}, which is too small to "
                f"reference against — every other free energy would be measured from "
                f"a state that was barely observed. Choose a populated reference."
            )

        rt = GAS_CONSTANT_KJ * self.temperature_k
        out: list[tuple[float, float | None, float | None]] = []
        for i, p in enumerate(self.probabilities):
            if p < NEGLIGIBLE_WEIGHT:
                out.append((math.inf, None, None))
                continue
            delta_g = -rt * math.log(p / p_ref)
            if self.uncertainty is None:
                out.append((delta_g, None, None))
                continue
            sigma = self.uncertainty[i]
            # Higher probability -> lower (more favourable) free energy, so the
            # upper end of the probability range gives the *minus* side.
            low = -rt * math.log((p + sigma) / p_ref)
            if p - sigma < NEGLIGIBLE_WEIGHT:
                out.append((delta_g, delta_g - low, None))
            else:
                high = -rt * math.log((p - sigma) / p_ref)
                out.append((delta_g, delta_g - low, high - delta_g))
        return out
