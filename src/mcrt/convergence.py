"""Pure helpers for the photon-count convergence study.

Monte Carlo statistical error scales as ~N^(-1/2). With no closed-form "truth"
for the noisy tail bins, we estimate the error of an observable as its spread
(standard deviation) across several independent-seed runs at fixed N. Plotting
that error vs N on log-log axes gives a clean -1/2 slope while statistics-limited;
it flattens once the systematic floor (finite bin width, tolerances) is hit. The
bend is the knee.

Within the swept range (N <= 1e6) none of our observables actually reach a
systematic floor — they all ride the -1/2 line — so ``find_knee`` honestly
returns None and the production recommendation is set by an absolute error
tolerance instead (``n_for_target_error``), extrapolating along the -1/2 line
when the target is not met in range (the low-mu tail-bin case).

These functions are deliberately pure (arrays in, numbers out) so the knee /
target-N / per-bin-error logic is unit-tested instead of buried in a script's
__main__.
"""

from typing import NamedTuple, Optional

import numpy as np

# A flat (floored) segment has a local log-log slope near 0. We accept |slope|
# below this band: shallower than the statistics-limited -1/2 on the low side,
# and not strongly *rising* on the high side (a rising error is noise, not a
# floor). Halfway between -1/2 and 0.
FLAT_SLOPE_BAND = 0.25

# A floor must be confirmed by at least this many consecutive flat segments — one
# flat pair can be 5-seed noise (e.g. a lone error point that ticks up at the
# largest N), not a real systematic floor.
MIN_FLAT_SEGMENTS = 2


class TargetN(NamedTuple):
    """Photons needed to reach an error target. ``extrapolated`` is True when the
    target was not met within the swept range and N was estimated along the
    fitted -1/2 line (so it is an estimate, not a measured N)."""
    n: float
    extrapolated: bool


def statistical_error(values, axis: int = 1) -> np.ndarray:
    """Sample standard deviation across the seed axis.

    This is the spread used as the statistical-error estimate when there is no
    closed-form truth. The study lays out arrays with seeds on **axis 1**:
    scalar observables ``(n_N, n_seeds)`` → ``(n_N,)`` and per-bin observables
    ``(n_N, n_seeds, n_bins)`` → ``(n_N, n_bins)``, so axis 1 is the default.
    """
    values = np.asarray(values, dtype=float)
    return np.std(values, axis=axis, ddof=1)


def loglog_slope(n, error) -> float:
    """Least-squares slope of log(error) vs log(n).

    Returns NaN if fewer than two finite, positive points are available.
    A purely statistics-limited observable yields a slope near -0.5.
    """
    n = np.asarray(n, dtype=float)
    error = np.asarray(error, dtype=float)
    mask = np.isfinite(n) & np.isfinite(error) & (n > 0) & (error > 0)
    if mask.sum() < 2:
        return float("nan")
    slope = np.polyfit(np.log(n[mask]), np.log(error[mask]), 1)[0]
    return float(slope)


def local_loglog_slopes(n, error) -> np.ndarray:
    """Slope between each consecutive (n, error) pair in log-log space.

    Returns an array of length ``len(n) - 1``. Pairs with a non-positive or
    non-finite error are reported as NaN (cannot take a log there).
    """
    n = np.asarray(n, dtype=float)
    error = np.asarray(error, dtype=float)
    with np.errstate(divide="ignore", invalid="ignore"):
        logn = np.log(n)
        loge = np.where(error > 0, np.log(error), np.nan)
        return np.diff(loge) / np.diff(logn)


def find_knee(n, error, flat_band: float = FLAT_SLOPE_BAND) -> Optional[float]:
    """Smallest N at which the log-log error curve reaches a *persistent* floor.

    While statistics-limited the local slope sits near -1/2; at the systematic
    floor it sits near 0 and *stays* there. We require, from the knee onward: a
    contiguous run of finite local slopes (no NaN gap), at least
    ``MIN_FLAT_SEGMENTS`` of them, every one within ``±flat_band`` of zero — so
    neither a falling segment, a lone flat segment, a single noisy uptick, nor a
    rising tail counts as a floor. The returned N is the *start* of that run (the
    first photon count already at the floor).

    Returns None if the curve never persistently flattens within the swept range
    (e.g. a low-μ tail bin still riding its -1/2 line at the largest N). That
    "no knee found" is itself the honest finding for the binding tail bins.
    """
    n = np.asarray(n, dtype=float)
    slopes = local_loglog_slopes(n, error)
    for i in range(len(slopes)):
        tail = slopes[i:]
        if np.any(~np.isfinite(tail)):
            continue  # a NaN gap after i → not a clean persistent floor from here
        if len(tail) >= MIN_FLAT_SEGMENTS and np.all(np.abs(tail) < flat_band):
            return float(n[i])  # error is already at floor at n[i]
    return None


def n_for_target_error(n, error, target: float) -> Optional[TargetN]:
    """Photons needed to bring an observable's error to ``target``.

    Solved on the *fitted* log-log line (error ∝ N**slope), not raw per-point
    errors — the fit is robust to 5-seed scatter and ignores degenerate
    zero-variance points (a tail bin with no escapers at small N has zero spread,
    which is undersampling, not convergence). ``extrapolated`` is True when the
    required N lies beyond the largest swept N (the target is not met in range).

    Returns None when the fit is unusable (fewer than two positive-error points,
    or a non-decreasing slope).
    """
    n = np.asarray(n, dtype=float)
    error = np.asarray(error, dtype=float)
    mask = np.isfinite(n) & np.isfinite(error) & (n > 0) & (error > 0)
    if mask.sum() < 2:
        return None
    slope, intercept = np.polyfit(np.log(n[mask]), np.log(error[mask]), 1)
    if not np.isfinite(slope) or slope >= 0:
        return None
    n_req = float(np.exp((np.log(target) - intercept) / slope))
    return TargetN(n_req, extrapolated=n_req > float(n[mask].max()))
