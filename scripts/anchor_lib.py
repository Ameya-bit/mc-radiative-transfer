"""Shared two-spot machinery for real-star beaming-swap anchors.

Both real-star anchor scripts — :mod:`j0030_anchor` (eclipsing, same-hemisphere
spots) and :mod:`j0740_anchor` (non-eclipsing, anti-phased spots) — reduce a
published fit to point spots and synthesize the star's pulse the same way:

    * each spot is one :func:`mcrt.compute_profile` call (the verified core),
    * a spot at longitude φ₀ is that profile rolled in phase (``np.roll``;
      separations land on the 1024-point grid), and
    * light is additive, so the star's flux is the weighted sum.

Swapping the (shared) ``beaming`` is therefore the *only* change between the
isotropic and realistic runs. Nothing here touches ``mcrt.pulse``.

The two scripts differ only in their anchors, headline metric, and plots; the
mechanics live here so they cannot drift apart.
"""

from typing import NamedTuple, Sequence

import numpy as np

from mcrt import compute_profile

LIBRARY_PATH = "data/beaming_library.npz"

N_PHASE = 1024
SHAPE_TAU = 10.0  # τ where the limb-darkening slope b(τ) peaks (Rung C result)
BACKGROUND_FRACS = np.array([0.0, 0.05, 0.10, 0.20, 0.30])


class Spot(NamedTuple):
    """One hot spot reduced to a point at its center colatitude.

    ``azimuth`` is the spot's longitude in cycles (only *relative* azimuth
    between a star's spots matters for the pulse; the absolute value is an
    arbitrary phase zero). ``weight`` is the relative bolometric amplitude
    ∝ (emitting area) × T_eff⁴.
    """

    colatitude: float  # θ_s, radians
    azimuth: float     # longitude, cycles
    weight: float      # relative amplitude ∝ area × T⁴


class Anchor(NamedTuple):
    """A published operating point: spacetime + multi-spot geometry."""

    label: str
    inclination: float    # i, radians
    compactness: float    # u = 2GM/Rc² (median)
    spots: Sequence[Spot]
    note: str             # weighting provenance / caveat


def multi_spot_flux(inclination, compactness, spots, beaming, n_phase=N_PHASE,
                    bending=None, rotation=None):
    """Weighted sum of point-spot fluxes; azimuth = integer phase roll.

    Light is additive, so the observed flux is the weighted sum of each spot's
    single-spot profile. A spot at longitude φ₀ produces the φ₀-shifted profile,
    which on a uniform full-cycle grid is ``np.roll`` by round(φ₀·n_phase). The
    same ``beaming`` is handed to every spot, so swapping it is the only change.

    ``bending`` selects the light-bending map, matching ``mcrt.compute_profile``:
    ``None`` (default) is Beloborodov's linear map — **bit-for-bit unchanged** from
    every prior anchor run — and an :class:`mcrt.bending.ExactBending` instance uses
    the exact Schwarzschild map (Track B2). It is handed identically to every spot,
    so exact-vs-linear is the only difference for a fixed geometry and ``beaming``.

    ``rotation`` (a :class:`mcrt.rotating.Rotation`, Track C) adds the Doppler +
    aberration layer to every spot; ``None`` (default) is the frozen slow-rotation
    flux, bit-for-bit unchanged. Each spot's β is set from its *own* colatitude
    (``spot_speed`` uses θ_s), and the roll places its whole — now fore-aft
    asymmetric — pulse at the spot's longitude, so the roll-and-add stays exact.
    Passed identically to every spot, so rotation-on-vs-off is the only difference
    for a fixed geometry, ``beaming``, and ``bending`` (Gate G2 coupling).
    """
    total = np.zeros(n_phase)
    for spot in spots:
        prof = compute_profile(inclination, spot.colatitude, compactness,
                               beaming=beaming, n_phase=n_phase, bending=bending,
                               rotation=rotation)
        shift = int(round(spot.azimuth * n_phase)) % n_phase
        total += spot.weight * np.roll(prof.flux, shift)
    return total


def single_spot_eclipsed_fraction(inclination, colatitude, compactness,
                                  n_phase=N_PHASE, bending=None):
    """Fraction of the rotation a *single* point spot spends below the horizon.

    Zero means the spot stays visible all rotation (J0740): then F_min > 0 and
    the pulsed fraction is unsaturated, so it can register the beaming swap. A
    large value (J0030) forces F_min = 0, pinning PF = 1 for every beaming and
    flattening the pulsed-fraction systematic to ΔPF ≈ 0.

    ``bending`` selects the light-bending map (``None`` = linear default; an
    :class:`mcrt.bending.ExactBending` uses the exact Schwarzschild horizon) so
    the visibility diagnostics match whatever map the flux is computed with.
    """
    prof = compute_profile(inclination, colatitude, compactness, n_phase=n_phase,
                           bending=bending)
    return float(np.mean(~prof.visible))


def single_spot_min_visible_mu(inclination, colatitude, compactness,
                               n_phase=N_PHASE, bending=None):
    """Smallest emission angle μ = cos α a *visible* single spot reaches.

    For a non-eclipsing but grazing geometry (J0740) this is the tiny but
    positive μ at the faint phase — the spot skims the limb without setting,
    which is why F_min is small yet non-zero. Returns ``nan`` if the spot is
    fully eclipsed at every phase (it never happens for the anchors used here).

    ``bending`` selects the light-bending map (``None`` = linear default; an
    :class:`mcrt.bending.ExactBending` uses the exact Schwarzschild map), keeping
    the grazing μ_min consistent with the flux's bending choice.
    """
    prof = compute_profile(inclination, colatitude, compactness, n_phase=n_phase,
                           bending=bending)
    if not np.any(prof.visible):
        return float("nan")
    return float(prof.cos_alpha[prof.visible].min())


def waveform_shape_change(iso_flux, real_flux):
    """Change between the two peak-normalized waveforms — the eclipse-immune metric.

    Returns (rms, max_local). Because each profile is divided by its own peak,
    this measures a *shape* difference and is insensitive to the absolute PF — it
    is the metric that still works when PF saturates (J0030), and a useful
    secondary read where PF does not (J0740).
    """
    ni = iso_flux / iso_flux.max()
    nr = real_flux / real_flux.max()
    diff = nr - ni
    return float(np.sqrt(np.mean(diff ** 2))), float(np.max(np.abs(diff)))


def delta_pf_vs_background(iso_flux, real_flux, bg_fracs=BACKGROUND_FRACS):
    """ΔPF as a common unpulsed background lifts the minimum of both curves.

    Adds a common absolute background B (a fraction of the isotropic peak) to
    both curves and reports ΔPF(B). Used as a *caveat* panel for the eclipsing
    star (where B un-saturates the pinned PF) and as a *robustness* panel for the
    non-eclipsing star (where a realistic background only dilutes the already-live
    ΔPF). We do not pick a B — this only shows the dependence.
    """
    f_scale = iso_flux.max()

    def pf(f):
        lo, hi = float(f.min()), float(f.max())
        return (hi - lo) / (hi + lo) if (hi + lo) > 0 else 0.0

    out = np.zeros(len(bg_fracs))
    for k, b in enumerate(bg_fracs):
        out[k] = pf(real_flux + b * f_scale) - pf(iso_flux + b * f_scale)
    return out


def load_library(path=LIBRARY_PATH):
    """Load the beaming library; return (tau_values, mu_centers, intensity_by_tau)."""
    d = np.load(path)
    return d["tau_values"], d["mu_centers"], d["intensity_by_tau"]


def shape_tau_index(tau_values, shape_tau=SHAPE_TAU):
    """Index of the profile-shape τ on the grid, raising if it is absent.

    Both anchor scripts snapshot a realistic waveform at ``SHAPE_TAU`` for their figures.
    Selecting it by an explicit index (rather than capturing whatever the loop last saw)
    keeps the two scripts identical and fails loudly if the library grid ever changes so
    that ``SHAPE_TAU`` is no longer sampled — no silent fallback to the wrong τ.
    """
    matches = np.where(np.asarray(tau_values) == shape_tau)[0]
    if matches.size == 0:
        raise ValueError(f"SHAPE_TAU={shape_tau} is not on the τ grid {tau_values!r}")
    return int(matches[0])


# Backwards-compatible alias: the J0030 script named the additive sum
# ``two_spot_flux``. It generalizes to any spot count, so the canonical name is
# ``multi_spot_flux``; keep the old name pointing at it.
two_spot_flux = multi_spot_flux
