"""Special-relativistic Doppler + aberration layer for a rotating hot spot.

The slow-rotation pulse pipeline (`pulse`) freezes the star: the spot's
*position* sweeps with phase, but its *velocity* is ignored. At the spin of a
real NICER target — J0740 at 346.5 Hz — the surface moves at β ≈ 0.12c, and that
motion (i) Doppler-boosts the observed flux and (ii) **aberrates** the emission
angle. The aberration is the term that *couples to the beaming swap*, because it
shifts the μ = cos α at which the beaming law I(μ) is sampled: the same rotation
that brightens the approaching phase also drags its emission angle, so the
isotropic→realistic systematic no longer factors out of the geometry. This module
adds that layer following the Schwarzschild + Doppler (S+D) recipe of Bogdanov
et al. (2019). It is **off by default** — pass a :class:`Rotation` to
:func:`mcrt.pulse.point_spot_flux` / :func:`mcrt.pulse.compute_profile` to opt in;
ν → 0 recovers the frozen pipeline bit-for-bit.

Physics (Bogdanov et al. 2019, ApJL 887, L26, §2; equation numbers theirs):

  - **Surface speed** (eq. 11):  ``β(θ) = 2πν R sinθ / (c √(1 − u))``, with the
    ``√(1 − u)`` the gravitational-redshift enhancement of the *locally measured*
    velocity (Ω = 2πν is the spin as seen at infinity). Lorentz ``γ = 1/√(1−β²)``.
  - **Ray/velocity angle** (static frame):  ``cos ξ = − sinα sin i sinφ / sinψ``,
    ξ between the photon direction and the azimuthal spot velocity. At φ = 0, π the
    motion is transverse (cos ξ = 0).
  - **Doppler factor** (eq. 12):  ``δ = 1 / [γ (1 − β cos ξ)]``.
  - **Aberration** (eq. 13):  ``cos α' = δ cos α`` — sample the beaming at the
    *comoving* cosine ``μ' = δ cos α``.
  - **Observed bolometric flux** (eq. 20 integrated over energy, with the eq. 8/14
    comoving-area factor):  ``F ∝ γ(θ) · δⁿ · I(μ') · μ' · D`` where ``D`` is the
    lensing Jacobian d cosα/d cosψ (from `pulse`/`bending`) and the projection
    cosine is the comoving ``μ'`` (via the Lorentz invariant dS cosα = dS' cosα',
    eq. 14). The exponent is ``n = 4`` for bolometric **energy** flux — the natural
    match for the grey, energy-integrated library I(μ) — and ``n = 3`` for the
    **photon**-number flux; :attr:`Rotation.photon_flux` selects the latter.

**Deliberately out of scope of the production flux** (see next-steps Track C /
Gate G2): photon light-travel-time delay across the star (eq. 18–19; ~2% of phase)
and rotational oblateness (~3% at 346 Hz, AlGendy & Morsink 2014). Both are
*beaming-independent* geometry terms, so they do not couple to the isotropic→
realistic swap at first order. The delay quadrature is provided here as a
driver-level utility (:func:`travel_time_delay`) — used by the SD1c validation and
the C4 caveat audit, never inside :func:`mcrt.pulse.point_spot_flux`.

References:
    Bogdanov et al. (2019), ApJL 887, L26 — the S+D recipe and the SD1 code-comparison suite.
    Poutanen & Beloborodov (2006), MNRAS 373, 836 — the underlying analytic S+D flux.
"""

from dataclasses import dataclass
from typing import Optional

import numpy as np

# Speed of light in km/s: with ν in Hz and R in km, 2πνR is a km/s speed directly.
SPEED_OF_LIGHT_KM_S = 299792.458

# Cap on E/kT exponents inside expm1: beyond this the Planck tail underflows to a
# clean 0 contribution anyway, and the cap avoids overflow warnings.
MAX_PLANCK_EXPONENT = 700.0


@dataclass(frozen=True)
class BandSpectrum:
    """A comoving blackbody spectrum observed through a fixed instrument band.

    Opting a :class:`Rotation` into this replaces the *bolometric* δⁿ boost with the
    exact in-band weight :func:`band_boost` — the band-limited generalization whose
    wide-band limits recover δ³ (photon) / δ⁴ (energy). Blackbody is a stand-in for
    the hydrogen-atmosphere spectra of the published fits; the angular (beaming)
    dependence stays whatever ``I(μ)`` the flux call passes, i.e. the spectrum is
    assumed separable, I′(E′, μ′) = Planck(E′; kT) × I(μ′) — the same separability
    the grey library already assumes.

    Attributes:
        kt_kev: comoving effective temperature kT of the spot (keV).
        e_min_kev / e_max_kev: instrument band edges in *observed* energy at
            infinity (keV) — e.g. NICER's calibrated 0.3–1.5 keV for J0740.
        n_quad: Gauss–Legendre nodes across the band (64 is converged to ≪1e-8).

    Raises:
        ValueError: on a non-positive temperature or an empty/inverted band.
    """

    kt_kev: float
    e_min_kev: float
    e_max_kev: float
    n_quad: int = 64

    def __post_init__(self) -> None:
        if not (self.kt_kev > 0.0 and 0.0 <= self.e_min_kev < self.e_max_kev):
            raise ValueError(
                f"invalid band: kT = {self.kt_kev} keV, "
                f"[{self.e_min_kev}, {self.e_max_kev}] keV — need kT > 0 and "
                "0 ≤ e_min < e_max"
            )


@dataclass(frozen=True)
class Rotation:
    """Rigid-rotation parameters for the Doppler + aberration layer.

    Attributes:
        spin_hz: rotation frequency ν as seen at infinity (Hz).
        radius_km: stellar circumferential radius R (km). With the compactness u
            (passed through the flux call) this fixes β via Bogdanov eq. 11.
        photon_flux: if ``True`` use the δ³ photon-number boost instead of the
            default δ⁴ bolometric-energy boost. The grey beaming library I(μ) is an
            energy intensity, so δ⁴ is the consistent default; δ³ is exposed for
            photon-flux cross-checks (e.g. NICER count-rate conventions).
        band: optional :class:`BandSpectrum`. When set, the bolometric δⁿ boost is
            replaced by the band-limited blackbody weight :func:`band_boost`
            (``photon_flux`` then selects in-band counts vs in-band energy;
            counts + a NICER band is the observable-faithful pairing). ``None``
            keeps the bolometric grey treatment.
    """

    spin_hz: float
    radius_km: float
    photon_flux: bool = False
    band: Optional[BandSpectrum] = None

    @property
    def flux_exponent(self) -> int:
        """δ-power in the observed flux: 3 (photon number) or 4 (bolometric energy)."""
        return 3 if self.photon_flux else 4


def band_boost(delta, compactness: float, band: BandSpectrum,
               photon_flux: bool = True):
    """In-band Doppler weight Φ_k(δ)/Φ_k(1) — the band-limited δⁿ replacement.

    For a comoving Planck spectrum, every Doppler and gravitational factor in the
    observed spectral flux (Bogdanov eq. 20 with E′ = E/(δ√(1−u)), eq. 15) collapses
    into a single *effective observed temperature*: the δ³ intensity boost exactly
    cancels the E′³ Planck numerator (the same collapse the SD1c validation uses),
    leaving the observed photon spectrum ∝ E² / (exp[E/(δ√(1−u)·kT)] − 1). So the
    in-band photon count (k = 2) and in-band energy flux (k = 3) are

        Φ_k(δ) = ∫_{E1}^{E2} E^k dE / (exp[E/(δ√(1−u)·kT)] − 1),

    and the phase-dependent flux weight replacing δⁿ is Φ_k(δ)/Φ_k(1) — normalized
    at δ = 1 so that ν → 0 (δ ≡ 1) stays bit-for-bit the frozen flux, and so equal-
    temperature multi-spot sums keep their relative weights. Wide-band limit:
    Φ_k ∝ (δ√(1−u)kT)^{k+1} ⇒ the ratio → δ³ (photon) / δ⁴ (energy). On a Wien-tail
    band (E ≫ δ√(1−u)kT, the J0740 regime) the weight is *exponentially* steeper
    than any power — the effective exponent d ln Φ/d ln δ exceeds 4.

    ``delta`` is scalar or array (from :func:`doppler_factor`); returns the same
    shape. ``photon_flux=True`` (counts) is the NICER-faithful convention.
    """
    exponent_k = 2 if photon_flux else 3
    nodes, weights = np.polynomial.legendre.leggauss(band.n_quad)
    energy = (0.5 * (band.e_max_kev - band.e_min_kev) * nodes
              + 0.5 * (band.e_max_kev + band.e_min_kev))
    t_eff_static = np.sqrt(1.0 - compactness) * band.kt_kev   # δ = 1 effective temp

    def phi(d: np.ndarray) -> np.ndarray:
        x = np.clip(energy / (d[..., None] * t_eff_static), None, MAX_PLANCK_EXPONENT)
        return (weights * energy**exponent_k / np.expm1(x)).sum(axis=-1)

    delta_arr = np.atleast_1d(np.asarray(delta, dtype=float))
    ratio = phi(delta_arr) / phi(np.ones(1))
    return ratio if np.ndim(delta) else float(ratio[0])


def spot_speed(spin_hz: float, radius_km: float, colatitude: float,
               compactness: float) -> float:
    """Local static-frame surface speed β = v/c at a spot (Bogdanov 2019 eq. 11).

    ``β = 2πν R sinθ / (c √(1 − u))``. The ``√(1 − u)`` is the gravitational-
    redshift enhancement of the velocity a local static observer measures (Ω = 2πν
    is the spin as seen at infinity). ``colatitude`` is θ_s (radians);
    ``compactness`` is ``u = R_s/R``. Returns dimensionless β.

    Raises:
        ValueError: if the resulting β is not in [0, 1) (super-luminal input).
    """
    beta = (2.0 * np.pi * spin_hz * radius_km * np.sin(colatitude)) / (
        SPEED_OF_LIGHT_KM_S * np.sqrt(1.0 - compactness))
    if not 0.0 <= beta < 1.0:
        raise ValueError(
            f"spot speed β = {beta:.4f} is not in [0, 1) for ν = {spin_hz} Hz, "
            f"R = {radius_km} km, θ = {colatitude} rad, u = {compactness}; check inputs"
        )
    return float(beta)


def lorentz_gamma(beta: float) -> float:
    """Lorentz factor γ = 1/√(1 − β²)."""
    return 1.0 / np.sqrt(1.0 - beta**2)


def cos_xi(phase, cos_alpha, cos_psi_vals, inclination: float) -> np.ndarray:
    """cos ξ between the photon direction and the spot velocity (static frame).

    ``cos ξ = − sinα sin i sinφ / sinψ`` (see the module derivation): the azimuthal
    (φ̂) velocity is perpendicular to the radial normal, so only the in-plane
    component sinα of the ray projects onto it, weighted by the observer's
    azimuthal projection sin i sinφ / sinψ. ``sinα`` and ``sinψ`` are taken from the
    cosines (both angles lie in [0, π], so their sines are ≥ 0).

    At φ = 0, π the motion is transverse and cos ξ = 0; those phases are also where
    sinψ → 0, and the numerator vanishes there too (sinφ = 0), so the guarded
    divide returns 0 rather than 0/0.
    """
    phase = np.asarray(phase, dtype=float)
    sin_alpha = np.sqrt(np.clip(1.0 - np.asarray(cos_alpha, dtype=float) ** 2, 0.0, 1.0))
    sin_psi = np.sqrt(np.clip(1.0 - np.asarray(cos_psi_vals, dtype=float) ** 2, 0.0, 1.0))
    numerator = -sin_alpha * np.sin(inclination) * np.sin(phase)

    out = np.zeros_like(numerator)
    nonzero = sin_psi > 1e-12
    out[nonzero] = numerator[nonzero] / sin_psi[nonzero]
    return np.clip(out, -1.0, 1.0)


def doppler_factor(beta: float, cos_xi_vals) -> np.ndarray:
    """Doppler factor δ = 1 / [γ (1 − β cos ξ)] (Bogdanov 2019 eq. 12).

    ``beta`` is the spot speed (scalar, from :func:`spot_speed`); ``cos_xi_vals`` is
    the ray/velocity cosine (array, from :func:`cos_xi`). δ > 1 blueshifts
    (approaching), δ < 1 redshifts; at cos ξ = 0 it is the transverse value 1/γ.
    """
    gamma = lorentz_gamma(beta)
    return 1.0 / (gamma * (1.0 - beta * np.asarray(cos_xi_vals, dtype=float)))


def travel_time_delay(cos_alpha, compactness: float, radius_km: float,
                      n_quad: int = 512) -> np.ndarray:
    """Photon travel-time delay Δt(α) vs. a radial ray, seconds (Bogdanov eq. 18).

    ``Δt = (R/c) ∫₀¹ dx / [x²(1 − u x)] · { [1 − β_b² x²(1 − u x)]^(−1/2) − 1 }``
    with ``β_b = sinα/√(1−u)`` the impact parameter b/R (eq. 17) and ``x = R/r``.
    The integrand is finite at x → 0 (→ β_b²/2); the only singularity is the
    integrable √ at x = 1 for grazing α = π/2, which the interior Gauss–Legendre
    nodes never hit. Flat-space limit (u = 0): Δt = (R/c)(1 − cos α) exactly.
    Eclipsed samples (cos α < 0) get Δt = 0; they carry no flux.

    **Driver-level utility, not production flux physics** — the observed-phase warp
    Δφ[cycles] = ν·Δt is geometry-only (beaming-independent): exactly PF-invariant
    for a single spot (a pure phase reparametrization preserves max/min) and
    second-order for multi-spot sums, which is what the C4 caveat audit bounds.
    """
    cos_alpha = np.asarray(cos_alpha, dtype=float)
    sin_alpha = np.sqrt(np.clip(1.0 - cos_alpha**2, 0.0, 1.0))
    beta_b = sin_alpha / np.sqrt(1.0 - compactness)             # b/R (eq. 17)

    t, w = np.polynomial.legendre.leggauss(n_quad)
    x = 0.5 * (t + 1.0)
    weights = 0.5 * w
    one_minus_ux = 1.0 - compactness * x
    radicand = 1.0 - (beta_b[..., None] ** 2) * x**2 * one_minus_ux
    integrand = ((1.0 / np.sqrt(np.clip(radicand, 1e-300, None))) - 1.0) / (x**2 * one_minus_ux)

    radius_m = radius_km * 1.0e3
    delay = (radius_m / (SPEED_OF_LIGHT_KM_S * 1.0e3)) * (integrand * weights).sum(axis=-1)
    return np.where(cos_alpha >= 0.0, delay, 0.0)
