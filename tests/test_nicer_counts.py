"""Tests for ``scripts/nicer_counts.py``: the adopted NICER count levels.

The numbers were digitised from Miller et al. (2019) Fig. 10 and Miller et al.
(2021) Fig. 6 and are pinned here so a change to the CSVs or the module is
caught. Targets are the one-cycle totals and ratios recorded in
``notes/revision-plan.md`` (Phase A) and the response letter, Point 2(iii).
"""

import sys
from pathlib import Path

import numpy as np
import pytest

SCRIPTS = Path(__file__).resolve().parents[1] / "scripts"
if str(SCRIPTS) not in sys.path:
    sys.path.insert(0, str(SCRIPTS))

import nicer_counts as nc


@pytest.fixture(params=["J0030", "J0740"])
def source(request):
    return nc.SOURCES[request.param]


def test_profile_lengths_match_published_binning(source):
    assert source.spot_counts.shape == (source.n_bins,)
    assert source.phase.shape == (source.n_bins,)
    assert np.all(np.diff(source.phase) > 0)


def test_total_is_spot_plus_background(source):
    np.testing.assert_allclose(
        source.total_counts, source.spot_counts + source.background_per_bin
    )


@pytest.mark.parametrize(
    "name, cycle_total, cycle_spot, cycle_bg",
    [("J0030", 1_341_407, 301_298, 1_040_109), ("J0740", 549_209, 52_188, 497_021)],
)
def test_one_cycle_totals_reproduce_phase_a_table(name, cycle_total, cycle_spot, cycle_bg):
    s = nc.SOURCES[name]
    assert s.total_counts.sum() == pytest.approx(cycle_total, rel=1e-3)
    assert s.spot_counts.sum() == pytest.approx(cycle_spot, rel=1e-3)
    assert s.background_per_bin * s.n_bins == pytest.approx(cycle_bg, rel=1e-3)


@pytest.mark.parametrize(
    "name, b_over_s, pf_total, pf_spot",
    [("J0030", 3.45, 0.200, 0.86), ("J0740", 9.5, 0.034, 0.37)],
)
def test_background_ratio_and_pulsed_fractions(name, b_over_s, pf_total, pf_spot):
    s = nc.SOURCES[name]
    assert s.background_to_mean_spot == pytest.approx(b_over_s, abs=0.05)
    assert nc.pulsed_fraction(s.total_counts) == pytest.approx(pf_total, abs=0.005)
    assert nc.pulsed_fraction(s.spot_counts) == pytest.approx(pf_spot, abs=0.01)


def test_j0030_count_rate_matches_miller_quoted_rate():
    s = nc.SOURCES["J0030"]
    # Miller et al. (2019) §2.2: long-term average ~0.7 counts/s.
    assert s.total_counts.sum() / s.exposure_s == pytest.approx(0.69, abs=0.02)


def test_j0740_3c50_background_is_the_high_bound():
    # Salmi et al. (2022) 3C50 model: B/<S> ≈ 12, above Miller's free-background 9.5.
    assert nc.J0740_3C50_BACKGROUND_TO_MEAN_SPOT > nc.SOURCES["J0740"].background_to_mean_spot
    assert nc.J0740_3C50_BACKGROUND_TO_MEAN_SPOT == pytest.approx(12, abs=1)


def test_sources_are_immutable(source):
    with pytest.raises(AttributeError):
        source.n_bins = 1  # frozen dataclass
    assert not source.spot_counts.flags.writeable
