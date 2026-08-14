"""A sigma scale is two numbers and a claim about where they came from.

``to_sigma`` and ``to_absolute`` used to take a :class:`MapHeader`, so they
converted against the *file's own* mean and RMS and there was no way to hand
them anything else. protean found why that matters: Mol\\*'s ``grid.stats`` —
the four numbers a viewer reports for a volume — are, for CCP4/MRC, **stored
header fields passed straight through unexamined**. A fixture written with
deliberately false header statistics failed on its first run reporting the
header's fake minimum *with the dimensions correct*, so the volume had genuinely
parsed and every number describing it was the file's claim rather than its
contents.

A host that has walked the voxels holds trustworthy statistics and, before
this, could not give them to these functions. That is the blockage this module
tests the fix for.

The stakes are not subtle. Mol\\*'s own default isosurface is ``relative: 2``,
computed as ``relativeValue * sigma + mean`` against exactly those statistics,
so a map with a stale header puts the surface in the wrong place in any viewer
and it looks entirely normal.
"""

from __future__ import annotations

import pytest
from test_mapinfo import write_map

from wiggles_em.density import MapStats, StatsSource, to_absolute, to_sigma
from wiggles_em.mapinfo import read_map_header


@pytest.fixture
def stale(tmp_path):
    """A header whose statistics do not describe the data.

    The shape of a real failure: a map that has been cropped or rescaled keeps
    whatever header nobody updated. Here the header claims mean 0, rms 1 while
    the data actually has mean 10, rms 4.
    """
    header = read_map_header(write_map(tmp_path, "stale.mrc", dmean=0.0, rms=1.0))
    measured = MapStats.measured(mean=10.0, rms=4.0)
    return header, measured


def test_stated_and_measured_disagree_when_the_header_is_stale(stale):
    """The whole point. Same absolute level, two answers, and the difference is
    the header being wrong rather than anything about the level."""
    header, measured = stale

    assert to_sigma(MapStats.stated(header), 18.0) == pytest.approx(18.0)
    assert to_sigma(measured, 18.0) == pytest.approx(2.0)


def test_measured_statistics_can_be_supplied_at_all(stale):
    """Before this they could not be: the signature took a MapHeader, so the
    file's claim was the only thing convertible."""
    _header, measured = stale
    assert to_absolute(measured, 2.0) == pytest.approx(18.0)


def test_a_scale_says_where_its_numbers_came_from(stale):
    """`Unit` exists so a level cannot be a bare number; this exists so a
    *scale* cannot be two bare numbers. A viewer reporting header fields as
    though it had measured them is the failure that motivated it."""
    header, measured = stale
    assert MapStats.stated(header).source is StatsSource.STATED
    assert measured.source is StatsSource.MEASURED


def test_round_trip_holds_for_measured_statistics(stale):
    _header, measured = stale
    assert to_sigma(measured, to_absolute(measured, 1.5)) == pytest.approx(1.5)


def test_a_non_positive_rms_is_still_refused_however_it_arrived():
    """`usable` is the same rule wherever the numbers came from. MRC writes
    rms=-1 for "statistics not computed" and it divides perfectly cleanly."""
    for rms in (0.0, -1.0):
        scale = MapStats.measured(mean=0.0, rms=rms)
        assert not scale.usable
        with pytest.raises(ValueError, match="sigma"):
            to_sigma(scale, 1.0)
        with pytest.raises(ValueError, match="sigma"):
            to_absolute(scale, 1.0)


def test_measured_statistics_rescue_a_map_a_header_would_refuse(tmp_path):
    """A map with rms=-1 cannot be converted from its header at all — but a
    host that walked the voxels has a usable scale for the very same map. The
    refusal was a fact about the header, not about the data."""
    header = read_map_header(write_map(tmp_path, "nostats.mrc", rms=-1.0))

    assert not MapStats.stated(header).usable
    with pytest.raises(ValueError):
        to_sigma(MapStats.stated(header), 1.0)

    measured = MapStats.measured(mean=2.0, rms=0.5)
    assert to_sigma(measured, 3.0) == pytest.approx(2.0)
