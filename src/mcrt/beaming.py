"""Turn Monte Carlo escape angles into the emergent specific intensity I(μ)."""

from typing import Callable

import numpy as np

# A beaming law maps emission cosine μ = cos α (array) to specific intensity I(μ).
BeamingFunc = Callable[[np.ndarray], np.ndarray]


def extract_intensity(escaped_mu, n_bins=20):
    """Convert escape angles into the emergent specific intensity I(μ).

    The raw count of photons escaping per μ-bin is the emergent *flux*, which is
    proportional to I(μ)·μ (a photon escaping at angle θ contributes μ = cos θ to
    the normal flux). Dividing the binned counts by the bin-center μ recovers the
    specific intensity I(μ) that the Eddington and Chandrasekhar laws describe.

    Returns:
        (mu_centers, intensity) — intensity normalized to 1 at the μ≈1 bin.
    """
    counts, edges = np.histogram(escaped_mu, bins=n_bins, range=(0.0, 1.0))
    centers = 0.5 * (edges[:-1] + edges[1:])
    intensity = counts / centers          # flux → specific intensity: divide by μ
    intensity = intensity / intensity[-1]  # normalize to I(μ≈1) = 1
    return centers, intensity


def fit_limb_darkening_slope(mu_centers, intensity, mu_floor=0.1):
    """Best-fit linear law I(μ) = a(1 + bμ); return the normalized slope b.

    The μ→0 bins are excluded (mu_floor): grazing escapes are rare, so dividing
    their tiny counts by a tiny μ is dominated by noise. Eddington predicts b = 1.5.
    """
    mask = mu_centers > mu_floor
    A = np.vstack([np.ones(mask.sum()), mu_centers[mask]]).T
    a, slope = np.linalg.lstsq(A, intensity[mask], rcond=None)[0]
    return slope / a


def beaming_lookup(mu_centers, intensity) -> BeamingFunc:
    """Build a callable ``I(μ)`` that linearly interpolates a tabulated curve.

    ``mu_centers`` (ascending) and ``intensity`` are one optical-depth row of the
    beaming library (``data/beaming_library.npz``: ``mu_centers`` and a row of
    ``intensity_by_tau``). The returned function evaluates ``I`` at an arbitrary
    emission cosine ``μ = cos α`` — exactly the ``beaming`` callable that
    :func:`mcrt.pulse.point_spot_flux` expects for the isotropic-vs-realistic
    comparison, where only the brightness term changes and the geometry is held fixed.

    Outside the tabulated μ range the curve is **held flat at its end values**
    (``np.interp`` clamping). That is the honest choice: the grazing μ→0 tail is
    the noisiest part of the library — the very bins :func:`fit_limb_darkening_slope`
    excludes below μ≈0.1 — so extrapolating it would amplify Monte Carlo noise,
    and μ never exceeds the μ≈1 normalization bin in practice.
    """
    mu_centers = np.asarray(mu_centers, dtype=float)
    intensity = np.asarray(intensity, dtype=float)

    def beaming(mu):
        return np.interp(np.asarray(mu, dtype=float), mu_centers, intensity)

    return beaming
