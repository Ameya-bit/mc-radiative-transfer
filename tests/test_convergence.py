import numpy as np
import pytest

from mcrt.convergence import (
    statistical_error,
    loglog_slope,
    local_loglog_slopes,
    find_knee,
    n_for_target_error,
)


def test_statistical_error_matches_sample_std():
    """Error estimate is the sample std (ddof=1) across the seed axis."""
    values = np.array([[1.0, 2.0, 3.0], [10.0, 10.0, 10.0]])  # (n_N=2, n_seeds=3)
    err = statistical_error(values, axis=-1)
    assert np.allclose(err, [np.std([1, 2, 3], ddof=1), 0.0])


def test_statistical_error_reduces_per_bin_axis():
    """Per-bin observable (n_N, n_seeds, n_bins) reduces over the seed axis."""
    values = np.zeros((4, 5, 7))
    values[:, 0, :] = 1.0  # one seed offset → nonzero spread
    err = statistical_error(values, axis=1)
    assert err.shape == (4, 7)
    assert np.all(err > 0)


def test_loglog_slope_recovers_minus_half():
    """A pure N^(-1/2) statistical error gives a log-log slope of -0.5."""
    n = np.array([1e3, 3e3, 1e4, 3e4, 1e5, 3e5, 1e6])
    error = n ** (-0.5)
    assert loglog_slope(n, error) == pytest.approx(-0.5, abs=1e-9)


def test_loglog_slope_nan_when_insufficient_points():
    assert np.isnan(loglog_slope([1e3], [0.1]))
    assert np.isnan(loglog_slope([1e3, 1e4], [0.0, 0.0]))  # zero errors masked out


def test_local_loglog_slopes_length_and_values():
    n = np.array([1.0, 10.0, 100.0])
    error = np.array([1.0, 0.1, 0.01])  # exact -1 slope in log-log
    slopes = local_loglog_slopes(n, error)
    assert slopes.shape == (2,)
    assert np.allclose(slopes, [-1.0, -1.0])


def test_find_knee_detects_persistent_floor():
    """Error that falls as N^(-1/2) then stays on a floor → knee at the floor's start."""
    n = np.array([1e3, 1e4, 1e5, 1e6, 1e7])
    error = np.array([3e-2, 1e-2, 1e-3, 1e-3, 1e-3])  # at floor from 1e5 onward
    knee = find_knee(n, error)
    assert knee == pytest.approx(1e5)  # first N already on the floor


def test_find_knee_none_for_pure_power_law():
    """A tail bin still riding the -1/2 line never flattens → no knee."""
    n = np.array([1e3, 3e3, 1e4, 3e4, 1e5, 3e5, 1e6])
    error = n ** (-0.5)
    assert find_knee(n, error) is None


def test_find_knee_ignores_lone_noisy_flat_pair():
    """A single flat segment amid an otherwise falling curve is not a knee."""
    # falls, one noisy flat pair (1e4→3e4 unchanged), then resumes falling
    n = np.array([1e3, 1e4, 3e4, 1e5, 1e6])
    error = np.array([3e-2, 1e-2, 1e-2, 3e-3, 1e-3])
    assert find_knee(n, error) is None


def test_find_knee_ignores_lone_final_uptick():
    """A single noisy uptick at the largest N (one flat/positive pair) is not a knee."""
    n = np.array([1e3, 1e4, 1e5, 1e6])
    error = np.array([1e-2, 3e-3, 1e-3, 1.3e-3])  # rises at the end from noise
    assert find_knee(n, error) is None


def test_find_knee_ignores_rising_tail():
    """Two consecutive *rising* segments are noise, not a floor — must not be a knee."""
    n = np.array([1e3, 1e4, 1e5, 1e6, 2e6])
    error = np.array([1e-2, 3e-3, 1e-3, 1.3e-3, 1.6e-3])  # rises over the last two pairs
    assert find_knee(n, error) is None


def test_n_for_target_error_met_in_range():
    """Target reachable within the swept range → not extrapolated."""
    n = np.array([1e3, 1e4, 1e5, 1e6])
    error = 1.0 * n ** (-0.5)  # clean -1/2 line, C = 1
    result = n_for_target_error(n, error, target=2e-3)
    assert result.extrapolated is False
    assert result.n == pytest.approx(2.5e5, rel=1e-6)  # (1 / 2e-3)^2


def test_n_for_target_error_extrapolates_beyond_range():
    """Tail bin still above target at N_max → extrapolate along the -1/2 line."""
    n = np.array([1e3, 1e4, 1e5, 1e6])
    error = 1.0 * n ** (-0.5)  # error_max = 1e-3 at 1e6
    result = n_for_target_error(n, error, target=5e-4)
    assert result.extrapolated is True
    assert result.n == pytest.approx(4e6, rel=1e-6)  # (1 / 5e-4)^2
