"""Point-spot pulse-profile synthesis in the slow-rotation Schwarzschild regime.

A neutron star has a hot spot; the star spins; an observer measures brightness
vs. rotation **phase** φ — the **pulse profile**. This module builds that profile
for a single point-like spot and reduces it to the **pulsed fraction**
``PF = (F_max − F_min)/(F_max + F_min)``.

Three deterministic ingredients (no Monte Carlo here — the stochastic part lives
upstream in the beaming library):

1. **Geometry.** The angle ψ between the spot normal and the line of sight as the
   star rotates:  ``cos ψ(φ) = cos i cos θ_s + sin i sin θ_s cos φ``.
2. **Light bending — Beloborodov (2002).** Gravity curves rays, so the emission
   angle α (from the local radial normal) maps to ψ by the linear approximation
   ``cos α = u + (1 − u) cos ψ``, with compactness ``u ≡ R_s/R = 2GM/(Rc²)``.
   The Jacobian ``d(cos α)/d(cos ψ) = (1 − u)`` is constant — that is what makes
   the map cheap and the verification clean. Because the map shifts cos ψ upward,
   the spot stays visible (cos α ≥ 0) for some ψ > 90°: we "see around the back".
3. **Beaming.** The local surface brightness ``I(μ)`` at emission cosine
   ``μ = cos α``. The default here is **isotropic** (``I ≡ 1``); the Rung C result
   swaps in the scattering beaming function ``I(μ; τ)`` from
   ``data/beaming_library.npz`` by passing a ``beaming`` callable. Only that
   brightness term changes — the geometry above is shared, which is what makes the
   isotropic-vs-realistic comparison a controlled experiment.

The proportional (slow-rotation) bolometric flux from a point spot is

    F(φ) ∝ (1 − u) · I(cos α) · cos α      where visible, else 0

(see Poutanen & Beloborodov 2006 for the full expression). Fluxes are returned in
arbitrary units — only ratios (the pulse shape and PF) are physical.

References:
    Beloborodov (2002), ApJ 566, L85 — the bending approximation and Rung A benchmark.
    Poutanen & Beloborodov (2006), MNRAS 373, 836 — the bolometric point-spot flux.
"""

from typing import Callable, NamedTuple, Optional

import numpy as np

# A beaming law maps emission cosine μ = cos α (array) to specific intensity I(μ).
BeamingFunc = Callable[[np.ndarray], np.ndarray]


class PulseProfile(NamedTuple):
    """One rotation of a point spot. Arrays are aligned on ``phase``.

    ``flux`` is proportional (arbitrary units); ``cos_alpha`` is the local
    emission cosine μ at each phase (only meaningful where ``visible``).
    """
    phase: np.ndarray      # φ ∈ [0, 2π), radians
    flux: np.ndarray       # F(φ), proportional, zero where the spot has set
    cos_alpha: np.ndarray  # μ = cos α at each phase
    visible: np.ndarray    # bool mask: spot above the bending horizon


def cos_psi(phase, inclination: float, colatitude: float) -> np.ndarray:
    """cos of the angle between spot normal and line of sight, vs. phase.

    ``cos ψ(φ) = cos i · cos θ_s + sin i · sin θ_s · cos φ``. Angles in radians:
    ``inclination`` = observer i, ``colatitude`` = spot θ_s. The extremes are
    cos(i − θ_s) at φ = 0 (spot nearest) and cos(i + θ_s) at φ = π (farthest).
    """
    phase = np.asarray(phase, dtype=float)
    return (
        np.cos(inclination) * np.cos(colatitude)
        + np.sin(inclination) * np.sin(colatitude) * np.cos(phase)
    )


def bend(cos_psi_vals, compactness: float) -> np.ndarray:
    """Beloborodov (2002) light bending: cos α = u + (1 − u) cos ψ.

    Maps the geometric angle ψ to the local emission angle α. ``compactness`` is
    ``u = R_s/R ∈ [0, 1)``; u = 0 is flat space (cos α = cos ψ). The slope in
    cos ψ is the constant Jacobian (1 − u).
    """
    cos_psi_vals = np.asarray(cos_psi_vals, dtype=float)
    return compactness + (1.0 - compactness) * cos_psi_vals


def visibility_threshold(compactness: float) -> float:
    """Smallest cos ψ that is still visible: −u/(1 − u).

    The spot is visible when cos α ≥ 0, i.e. cos ψ ≥ −u/(1 − u). This threshold is
    negative, so ψ > 90° (the near-far hemisphere boundary) can remain visible —
    the gravitational "seeing around the star" effect.
    """
    return -compactness / (1.0 - compactness)


def point_spot_flux(
    phase,
    inclination: float,
    colatitude: float,
    compactness: float,
    beaming: Optional[BeamingFunc] = None,
) -> np.ndarray:
    """Proportional bolometric flux F(φ) from a point spot, zero where it has set.

    ``F(φ) ∝ (1 − u) · I(cos α) · cos α`` for visible phases, else 0. ``beaming``
    is the surface brightness law ``I(μ)`` evaluated at ``μ = cos α``; ``None``
    means isotropic (``I ≡ 1``). The geometry is identical regardless of
    ``beaming`` — passing the library's ``I(μ; τ)`` is the only change Rung C makes.
    """
    cos_a = bend(cos_psi(phase, inclination, colatitude), compactness)
    visible = cos_a >= 0.0

    intensity = np.ones_like(cos_a) if beaming is None else np.asarray(beaming(cos_a), dtype=float)
    flux = (1.0 - compactness) * intensity * cos_a
    return np.where(visible, flux, 0.0)


def compute_profile(
    inclination: float,
    colatitude: float,
    compactness: float,
    beaming: Optional[BeamingFunc] = None,
    n_phase: int = 1024,
) -> PulseProfile:
    """Sample one full rotation onto a uniform phase grid.

    ``n_phase`` phases on [0, 2π) (endpoint excluded so φ = 0 and, for even
    ``n_phase``, φ = π land exactly on grid points — the flux extremes for an
    isotropic spot). Returns a :class:`PulseProfile` bundling phase, flux,
    cos α, and the visibility mask for plotting and downstream sweeps.
    """
    phase = np.linspace(0.0, 2.0 * np.pi, n_phase, endpoint=False)
    cos_a = bend(cos_psi(phase, inclination, colatitude), compactness)
    visible = cos_a >= 0.0
    flux = point_spot_flux(phase, inclination, colatitude, compactness, beaming)
    return PulseProfile(phase=phase, flux=flux, cos_alpha=cos_a, visible=visible)


def pulsed_fraction(flux) -> float:
    """``(F_max − F_min)/(F_max + F_min)``; 0 when the flux is identically zero.

    The standard bolometric pulse-strength summary. The all-zero guard covers a
    spot that never rises (no divide-by-zero); any visible spot has F_max > 0.
    """
    flux = np.asarray(flux, dtype=float)
    f_max, f_min = float(flux.max()), float(flux.min())
    total = f_max + f_min
    if total == 0.0:
        return 0.0
    return (f_max - f_min) / total


def analytic_isotropic_pf(inclination: float, colatitude: float, compactness: float) -> float:
    """Closed-form pulsed fraction for an **always-visible** isotropic point spot.

    With isotropic emission the flux is ``F ∝ (1−u)(u + (1−u) cos ψ)``, monotonic
    in cos ψ, so the extremes sit at φ = 0 and φ = π:

        PF = (1−u)(cosψ_max − cosψ_min) / (2u + (1−u)(cosψ_max + cosψ_min)),

    with cosψ_max = cos(i − θ_s), cosψ_min = cos(i + θ_s). This is the Rung A
    benchmark. It is only valid when the spot never sets — if cosψ_min falls below
    the bending visibility threshold the minimum is an eclipse (F_min = 0), the
    extremes are no longer at φ = 0, π, and this returns nothing useful, so it
    raises instead.

    Raises:
        ValueError: if the geometry eclipses the spot (cosψ_min < −u/(1−u)).
    """
    cpsi_max = np.cos(inclination - colatitude)
    cpsi_min = np.cos(inclination + colatitude)
    if cpsi_min < visibility_threshold(compactness):
        raise ValueError(
            "spot is eclipsed for this geometry; the always-visible closed form "
            "does not apply — use compute_profile + pulsed_fraction instead"
        )
    numerator = (1.0 - compactness) * (cpsi_max - cpsi_min)
    denominator = 2.0 * compactness + (1.0 - compactness) * (cpsi_max + cpsi_min)
    return float(numerator / denominator)
