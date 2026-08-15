"""Tests for state occupancies and the free energies derived from them.

The claim under test is not "these numbers are right" — that is the method's
job — but "a number never travels without the method that produced it, and the
two families of method are never confused". Evans *et al.* show latent
histogramming losing a genuinely present middle mode while deconvolution and
ensemble reweighting recover it, so the distinction between those families is
the load-bearing part of this module.
"""

from __future__ import annotations

import math
import re

import pytest

from wiggles_em.populations import (
    GAS_CONSTANT_KJ,
    NEGLIGIBLE_WEIGHT,
    Populations,
    WeightSource,
)

# ── the two families ─────────────────────────────────────────────────────────


def test_every_source_lands_in_exactly_one_family_or_neither():
    """`UNKNOWN` is its own category, not a middle ground.

    If a source were both per-image and inverse-solving, the readout would tell
    a user two contradictory things about the same number.
    """
    for source in WeightSource:
        assert not (source.is_per_image_assignment and source.solves_inverse_problem)
    assert WeightSource.UNKNOWN.is_per_image_assignment is False
    assert WeightSource.UNKNOWN.solves_inverse_problem is False


def test_the_families_match_the_paper_they_came_from():
    counting = {
        WeightSource.THREE_D_CLASSIFICATION,
        WeightSource.HARD_ASSIGNMENT,
        WeightSource.SOFT_ASSIGNMENT,
        WeightSource.LATENT_HISTOGRAM,
    }
    assert {s for s in WeightSource if s.is_per_image_assignment} == counting
    assert {s for s in WeightSource if s.solves_inverse_problem} == {
        WeightSource.DECONVOLVED,
        WeightSource.ENSEMBLE_REWEIGHTED,
    }


def test_every_source_has_its_own_caveat():
    """Boilerplate caveats would make the label decorative."""
    caveats = {s: s.caveat for s in WeightSource}
    assert len(set(caveats.values())) == len(caveats)
    for source, caveat in caveats.items():
        assert len(caveat) > 60, f"{source} has a stub caveat"


def test_only_inverse_solving_weights_are_quantitative():
    for source in WeightSource:
        pops = Populations.declare([1.0, 1.0], source)
        assert pops.is_quantitative is source.solves_inverse_problem


def test_a_non_quantitative_readout_says_so_and_says_what_would_fix_it():
    pops = Populations.declare([0.8, 0.2], WeightSource.LATENT_HISTOGRAM)
    banner = pops.banner()

    assert "NOT QUANTITATIVE" in banner
    assert "deconvolution" in banner and "reweighting" in banner
    assert "flat distribution" in banner, "the specific known failure, not a vague warning"


def test_undeclared_weights_are_unknown_rather_than_assumed_good():
    pops = Populations.declare([0.5, 0.5])

    assert pops.source is WeightSource.UNKNOWN
    assert pops.is_quantitative is False
    assert "NOT DECLARED" in pops.banner()


# ── normalisation ────────────────────────────────────────────────────────────


def test_counts_become_probabilities_and_the_original_total_is_kept():
    pops = Populations.declare([300.0, 100.0], WeightSource.HARD_ASSIGNMENT)

    assert pops.probabilities == pytest.approx((0.75, 0.25))
    assert pops.raw_total == pytest.approx(400.0)
    assert "Normalised from a total of 400" in pops.banner()


def test_already_normalised_weights_do_not_advertise_a_rescale():
    assert "Normalised from" not in Populations.declare([0.6, 0.4]).banner()


def test_uncertainty_is_rescaled_with_the_weights():
    """Otherwise an error bar quoted on counts lands on a probability axis."""
    pops = Populations.declare([300.0, 100.0], uncertainty=[30.0, 10.0])

    assert pops.uncertainty == pytest.approx((0.075, 0.025))


@pytest.mark.parametrize(
    "kwargs,message",
    [
        ({"weights": []}, "not a population"),
        ({"weights": [1.0, -0.2]}, "negative"),
        ({"weights": [0.0, 0.0]}, "sum to zero"),
        ({"weights": [1.0, math.nan]}, "NaN"),
        ({"weights": [1.0, 1.0], "uncertainty": [0.1]}, "parallel"),
        ({"weights": [1.0, 1.0], "uncertainty": [0.1, -0.1]}, "non-negative"),
        ({"weights": [1.0, 1.0], "temperature_k": 0.0}, "positive number of kelvin"),
    ],
)
def test_malformed_input_is_refused(kwargs, message):
    with pytest.raises(ValueError, match=message):
        Populations.declare(**kwargs)


def test_a_negative_weight_is_refused_rather_than_clamped():
    """A negative occupancy means an unconstrained fit, not a small population.

    Clamping to zero would produce a plausible distribution from a broken one.
    """
    with pytest.raises(ValueError, match="non-negativity constraint"):
        Populations.declare([0.9, -0.1, 0.2])


# ── free energy ──────────────────────────────────────────────────────────────


def test_free_energy_needs_a_declared_temperature():
    pops = Populations.declare([0.75, 0.25], WeightSource.ENSEMBLE_REWEIGHTED)
    with pytest.raises(ValueError, match="not the temperature of the grid"):
        pops.relative_free_energy()


def test_free_energy_matches_the_boltzmann_relation():
    pops = Populations.declare([0.75, 0.25], WeightSource.ENSEMBLE_REWEIGHTED, temperature_k=298.15)
    dg = pops.relative_free_energy()

    rt = GAS_CONSTANT_KJ * 298.15
    assert dg[0][0] == pytest.approx(0.0), "the reference state is its own zero"
    assert dg[1][0] == pytest.approx(-rt * math.log(0.25 / 0.75))
    assert dg[1][0] > 0, "the rarer state is the less favourable one"


def test_the_reference_state_can_be_chosen():
    pops = Populations.declare([0.75, 0.25], temperature_k=298.15)
    forward = pops.relative_free_energy(reference=0)[1][0]
    backward = pops.relative_free_energy(reference=1)[0][0]

    assert forward == pytest.approx(-backward)


def test_error_bars_are_asymmetric_and_the_upper_one_can_be_unbounded():
    """The point of modelling uncertainty here at all.

    A logarithm turns a symmetric error on a probability into a lopsided one on
    a free energy, and when the error reaches the probability the upper end
    stops existing. Those are the sparsely populated states this approach is
    for, so a symmetric bar would mislead exactly where it matters most.
    """
    # A rare state whose error bar is as wide as the state itself — the case
    # this whole approach is aimed at, and where "2.3 +/- 0.4 kJ/mol" would be
    # a fabrication.
    pops = Populations.declare(
        [0.98, 0.02],
        WeightSource.ENSEMBLE_REWEIGHTED,
        uncertainty=[0.02, 0.02],
        temperature_k=298.15,
    )
    _, minus, plus = pops.relative_free_energy()[1]

    assert minus is not None
    assert plus is None, "p - sigma reaches zero, so the upper bound is unbounded"

    # Halve that error and the upper bound exists again — the transition is a
    # property of the numbers, not a special case for one hard-coded input.
    narrower = Populations.declare([0.98, 0.02], uncertainty=[0.01, 0.01], temperature_k=298.15)
    assert narrower.relative_free_energy()[1][2] is not None

    # And with room to spare it is finite, still asymmetric, and the larger side
    # is the one running towards lower probability.
    comfortable = Populations.declare([0.6, 0.4], uncertainty=[0.05, 0.05], temperature_k=298.15)
    _, m, p = comfortable.relative_free_energy()[1]
    assert m is not None and p is not None
    assert p > m, "the low-probability side of the bar is the long one"


def test_no_uncertainty_means_no_error_bars_rather_than_zero_width_ones():
    pops = Populations.declare([0.6, 0.4], temperature_k=298.15)
    dg, minus, plus = pops.relative_free_energy()[1]

    assert dg == pytest.approx(-GAS_CONSTANT_KJ * 298.15 * math.log(0.4 / 0.6))
    assert minus is None and plus is None
    assert "No uncertainty was supplied" in pops.banner()


def test_a_negligible_state_gets_infinite_free_energy_not_a_confident_number():
    """Taking the log of a normalisation artefact would invent a measurement."""
    pops = Populations.declare([1.0, NEGLIGIBLE_WEIGHT / 10], temperature_k=298.15)

    assert pops.relative_free_energy()[1][0] == math.inf


def test_a_negligible_reference_is_refused():
    pops = Populations.declare([1.0, NEGLIGIBLE_WEIGHT / 10], temperature_k=298.15)
    with pytest.raises(ValueError, match="too small to reference against"):
        pops.relative_free_energy(reference=1)


def test_an_out_of_range_reference_is_refused():
    pops = Populations.declare([0.6, 0.4], temperature_k=298.15)
    with pytest.raises(ValueError, match=re.escape("outside 0..1")):
        pops.relative_free_energy(reference=2)
