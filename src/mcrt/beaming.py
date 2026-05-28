"""Turn Monte Carlo escape angles into the emergent specific intensity I(μ)."""

import numpy as np


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
