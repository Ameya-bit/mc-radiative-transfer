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
   Passing ``bending=`` an :class:`mcrt.bending.ExactBending` swaps this linear
   map for the exact Schwarzschild ray integral (map *and* Jacobian); the linear
   default is the untouched, verified path.
3. **Beaming.** The local surface brightness ``I(μ)`` at emission cosine
   ``μ = cos α``. The default here is **isotropic** (``I ≡ 1``); the
   isotropic-vs-realistic comparison swaps in the scattering beaming function ``I(μ; τ)`` from
   ``data/beaming_library.npz`` by passing a ``beaming`` callable. Only that
   brightness term changes — the geometry above is shared, which is what makes the
   isotropic-vs-realistic comparison a controlled experiment.

The proportional (slow-rotation) bolometric flux from a point spot is

    F(φ) ∝ D · I(cos α) · cos α      where visible, else 0

with lensing Jacobian ``D = d(cos α)/d(cos ψ)``: the constant ``(1 − u)`` for the
linear map, or the phase-dependent ``ExactBending.jacobian`` for the exact one
(see Poutanen & Beloborodov 2006 for the full expression). Fluxes are returned in
arbitrary units — only ratios (the pulse shape and PF) are physical.

An optional **fourth ingredient** — the star's rotational *velocity* — enters via
``bending=``'s sibling ``rotation=`` (a :class:`mcrt.rotating.Rotation`): the
special-relativistic Doppler boost and aberration of the S+D approximation
(Bogdanov et al. 2019), which multiply the flux by ``γ(θ) δⁿ`` and sample the
beaming at the aberrated cosine ``μ' = δ cos α``. It is off by default and the
ν → 0 limit is bit-for-bit the frozen flux above; see `mcrt.rotating`.

References:
    Beloborodov (2002), ApJ 566, L85 — the bending approximation and the analytic-check benchmark.
    Poutanen & Beloborodov (2006), MNRAS 373, 836 — the bolometric point-spot flux.
"""

from typing import NamedTuple, Optional

import numpy as np

# A beaming law maps emission cosine μ = cos α (array) to specific intensity I(μ).
# Defined in `beaming` (the dependency-free module that builds these curves) and
# reused here so the type has one home.
from .beaming import BeamingFunc

# The exact light-bending map (optional; the linear `bend` below is the default).
from .bending import ExactBending

# The optional rotational (Doppler + aberration) layer. Off unless a `Rotation` is
# passed; ν → 0 recovers the frozen slow-rotation flux below bit-for-bit. A
# `Rotation.band` (BandSpectrum) swaps the bolometric δⁿ boost for the band-limited
# blackbody weight (Track C3).
from .rotating import Rotation, band_boost, cos_xi, doppler_factor, spot_speed


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


def _bend_and_jacobian(cos_psi_vals, compactness: float, bending: Optional[ExactBending]):
    """(cos α, lensing Jacobian D) at each cos ψ, for the linear or exact map.

    ``bending=None`` uses Beloborodov's linear map (constant ``D = 1 − u``); an
    :class:`mcrt.bending.ExactBending` supplies the exact map and its phase-varying
    ``D``. Guards a compactness mismatch — an exact map built for a different ``u``
    than the flux is being computed at would silently corrupt the geometry.
    """
    if bending is None:
        cos_a = bend(cos_psi_vals, compactness)
        return cos_a, 1.0 - compactness
    if not np.isclose(bending.compactness, compactness):
        raise ValueError(
            f"bending map built for u = {bending.compactness} but flux requested "
            f"at u = {compactness}; build ExactBending with the matching compactness"
        )
    cos_a = bending.cos_alpha(cos_psi_vals)
    return cos_a, bending.jacobian(cos_psi_vals)


def point_spot_flux(
    phase,
    inclination: float,
    colatitude: float,
    compactness: float,
    beaming: Optional[BeamingFunc] = None,
    bending: Optional[ExactBending] = None,
    rotation: Optional[Rotation] = None,
) -> np.ndarray:
    """Proportional flux F(φ) from a point spot, zero where it has set.

    Frozen (slow-rotation) form: ``F(φ) ∝ D · I(cos α) · cos α`` for visible phases,
    else 0. ``beaming`` is the surface brightness law ``I(μ)`` evaluated at
    ``μ = cos α``; ``None`` means isotropic (``I ≡ 1``). ``bending`` selects the
    light-bending map: ``None`` is Beloborodov's linear default (``D = 1 − u``); an
    :class:`mcrt.bending.ExactBending` instance uses the exact Schwarzschild map and
    its Jacobian. The geometry is identical across ``beaming`` for a fixed
    ``bending`` — passing the library's ``I(μ; τ)`` is the only change the
    realistic-beaming comparison makes.

    ``rotation`` (a :class:`mcrt.rotating.Rotation`) adds the special-relativistic
    Doppler + aberration layer (Bogdanov et al. 2019 S+D recipe): the flux becomes
    ``F ∝ γ(θ) · δⁿ · I(μ') · μ' · D`` with Doppler factor δ(φ), aberrated cosine
    ``μ' = δ cos α`` (used for both the beaming lookup and the projection, via the
    Lorentz invariant dS cosα = dS' cosα'), and n = 4 (energy) or 3 (photon) flux
    boost. ``None`` (the default) leaves the frozen form untouched, bit-for-bit; the
    visibility mask ``cos α ≥ 0`` is a static-frame geometric condition and is
    unchanged by rotation (δ > 0 ⇒ μ' and cos α share sign).
    """
    cos_psi_vals = cos_psi(phase, inclination, colatitude)
    cos_a, jacobian = _bend_and_jacobian(cos_psi_vals, compactness, bending)
    visible = cos_a >= 0.0

    if rotation is None:
        mu, boost, area_weight = cos_a, 1.0, 1.0
    else:
        beta = spot_speed(rotation.spin_hz, rotation.radius_km, colatitude, compactness)
        delta = doppler_factor(beta, cos_xi(phase, cos_a, cos_psi_vals, inclination))
        mu = np.clip(delta * cos_a, -1.0, 1.0)              # μ' = cos α' (aberrated)
        if rotation.band is None:
            boost = delta ** rotation.flux_exponent         # δ⁴ energy / δ³ photon
        else:                                               # band-limited blackbody
            boost = band_boost(delta, compactness, rotation.band, rotation.photon_flux)
        area_weight = 1.0 / np.sqrt(1.0 - beta**2)          # γ(θ): comoving-area factor

    intensity = np.ones_like(cos_a) if beaming is None else np.asarray(beaming(mu), dtype=float)
    flux = area_weight * boost * jacobian * intensity * mu
    return np.where(visible, flux, 0.0)


def compute_profile(
    inclination: float,
    colatitude: float,
    compactness: float,
    beaming: Optional[BeamingFunc] = None,
    n_phase: int = 1024,
    bending: Optional[ExactBending] = None,
    rotation: Optional[Rotation] = None,
) -> PulseProfile:
    """Sample one full rotation onto a uniform phase grid.

    ``n_phase`` phases on [0, 2π) (endpoint excluded so φ = 0 and, for even
    ``n_phase``, φ = π land exactly on grid points — the flux extremes for an
    isotropic, non-rotating spot). ``bending`` selects the linear (``None``) or exact
    (:class:`mcrt.bending.ExactBending`) light-bending map; ``rotation`` (a
    :class:`mcrt.rotating.Rotation`) adds the Doppler + aberration layer to the flux
    (``None`` = frozen slow-rotation). Returns a :class:`PulseProfile` bundling
    phase, flux, cos α, and the visibility mask. ``cos_alpha`` is the *static-frame*
    geometric emission cosine (unchanged by rotation — the aberrated μ' lives inside
    the flux only), so the geometry stays a controlled constant across beaming and
    rotation swaps.
    """
    phase = np.linspace(0.0, 2.0 * np.pi, n_phase, endpoint=False)
    cos_a, _ = _bend_and_jacobian(
        cos_psi(phase, inclination, colatitude), compactness, bending)
    visible = cos_a >= 0.0
    flux = point_spot_flux(phase, inclination, colatitude, compactness, beaming, bending, rotation)
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

    with cosψ_max = cos(i − θ_s), cosψ_min = cos(i + θ_s). This is the analytic-check
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
