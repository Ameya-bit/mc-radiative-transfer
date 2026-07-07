"""Verification tests for the Doppler + aberration layer (`mcrt.rotating`).

Two load-bearing groups mirror the repo's testing philosophy — pin the physics to
something *algebraically independent* of the module under test:

1. **The explicit Lorentz-boost cross-check.** `cos_xi`, `doppler_factor`, and the
   aberration ``μ' = δ cos α`` are re-derived from first principles by building the
   photon and velocity 3-vectors and boosting the photon's null 4-vector into the
   comoving frame — sharing no code with the module's trig shortcuts.
2. **ν → 0 bit-for-bit.** Passing a zero-spin :class:`Rotation` through the flux must
   return the frozen slow-rotation flux *exactly*, so the layer is provably additive.

Plus an external anchor: the 200 Hz SD1c code-comparison waveform (Bogdanov et al.
2019), which the layer must reproduce (with the geometry-only travel-time delay
applied) to the same ~1% linear-bending accuracy as the static SD1a case.
"""

from pathlib import Path

import numpy as np
import pytest

from mcrt import (
    BandSpectrum,
    Rotation,
    band_boost,
    bend,
    compute_profile,
    cos_psi,
    cos_xi,
    doppler_factor,
    lorentz_gamma,
    point_spot_flux,
    spot_speed,
    travel_time_delay,
)

SPEED_OF_LIGHT_KM_S = 299792.458

# An always-visible geometry (spot never sets) so flux symmetry is clean to reason
# about: non-rotating flux is even in φ, and rotation is the only thing that breaks it.
ALWAYS_VISIBLE = dict(inclination=np.deg2rad(30.0), colatitude=np.deg2rad(20.0), compactness=0.3)

# J0740 Riley compactness/geometry, for realistic β magnitudes.
J0740_U = 0.494
J0740_R_KM = 12.39
J0740_SPIN_HZ = 346.5


# --- spot_speed (Bogdanov eq. 11) --------------------------------------------

def test_spot_speed_matches_closed_form_inline():
    """β = 2πν R sinθ / (c √(1−u)), recomputed independently of the module."""
    nu, r, theta, u = 346.5, 12.39, np.deg2rad(77.0), 0.494
    expected = (2.0 * np.pi * nu * r * np.sin(theta)) / (SPEED_OF_LIGHT_KM_S * np.sqrt(1.0 - u))
    assert spot_speed(nu, r, theta, u) == pytest.approx(expected)


def test_spot_speed_j0740_is_order_tenth_c():
    """J0740's equatorial surface moves at ≈ 0.12–0.13 c — the regime that matters.

    (Guards against silently dropping the √(1−u) gravitational enhancement, which
    would give a spuriously low ≈ 0.09.)
    """
    beta_eq = spot_speed(J0740_SPIN_HZ, J0740_R_KM, np.deg2rad(90.0), J0740_U)
    assert 0.11 < beta_eq < 0.14
    # The √(1−u) factor raises β above the naive 2πνR sinθ/c:
    naive = 2.0 * np.pi * J0740_SPIN_HZ * J0740_R_KM / SPEED_OF_LIGHT_KM_S
    assert beta_eq > naive


def test_spot_speed_grows_with_sin_colatitude():
    """β ∝ sinθ: an equatorial spot is faster than a polar one."""
    eq = spot_speed(600.0, 12.0, np.deg2rad(90.0), 0.3)
    polar = spot_speed(600.0, 12.0, np.deg2rad(10.0), 0.3)
    assert eq > polar > 0.0


def test_spot_speed_rejects_superluminal():
    """An absurd spin that drives β ≥ 1 raises rather than returning nonsense."""
    with pytest.raises(ValueError):
        spot_speed(5000.0, 30.0, np.deg2rad(90.0), 0.5)


# --- lorentz_gamma / doppler_factor ------------------------------------------

def test_lorentz_gamma_basic():
    assert lorentz_gamma(0.0) == pytest.approx(1.0)
    assert lorentz_gamma(0.6) == pytest.approx(1.25)   # 1/√(1−0.36)


def test_doppler_factor_matches_closed_form():
    """δ = 1/[γ(1 − β cos ξ)], recomputed inline."""
    beta = 0.1
    cxi = np.array([-0.5, 0.0, 0.3, 0.9])
    expected = 1.0 / ((1.0 / np.sqrt(1.0 - beta**2)) * (1.0 - beta * cxi))
    assert np.allclose(doppler_factor(beta, cxi), expected)


def test_doppler_is_transverse_at_phase_0_and_pi():
    """At φ = 0, π the velocity is transverse (cos ξ = 0), so δ = 1/γ (redshift)."""
    g = ALWAYS_VISIBLE
    phase = np.array([0.0, np.pi])
    cps = cos_psi(phase, g["inclination"], g["colatitude"])
    cosa = bend(cps, g["compactness"])
    cxi = cos_xi(phase, cosa, cps, g["inclination"])
    assert np.allclose(cxi, 0.0, atol=1e-12)
    beta = spot_speed(400.0, 12.0, g["colatitude"], g["compactness"])
    assert np.allclose(doppler_factor(beta, cxi), 1.0 / lorentz_gamma(beta))


# --- the explicit Lorentz-boost cross-check (algebraically independent) -------

def _independent_delta_and_mu_prime(phase, inclination, colatitude, compactness, beta):
    """δ and μ' = cos α' from explicit 3-vectors + a null-4-vector Lorentz boost.

    Builds the radial normal n̂, azimuthal velocity direction v̂, observer direction
    k̂_obs, and the photon's initial direction k̂₀ = cosα n̂ + sinα t̂ (α from linear
    bending), then boosts the photon 4-vector (1, k̂₀) into the frame moving at β v̂.
    δ = 1/E' and cos α' = n̂·k̂₀'. Shares no code with `mcrt.rotating`.
    """
    i = inclination
    k_obs = np.array([np.sin(i), 0.0, np.cos(i)])
    deltas, mus = [], []
    for phi in np.atleast_1d(phase):
        n_hat = np.array([np.sin(colatitude) * np.cos(phi),
                          np.sin(colatitude) * np.sin(phi),
                          np.cos(colatitude)])
        v_hat = np.array([-np.sin(phi), np.cos(phi), 0.0])          # +φ̂ (rotation sense)
        cos_ps = float(n_hat @ k_obs)
        cos_a = compactness + (1.0 - compactness) * cos_ps          # linear bending
        sin_ps = np.sqrt(max(1.0 - cos_ps**2, 0.0))
        if sin_ps < 1e-12:
            t_hat = np.zeros(3)
        else:
            t_hat = (k_obs - cos_ps * n_hat) / sin_ps
        sin_a = np.sqrt(max(1.0 - cos_a**2, 0.0))
        k0 = cos_a * n_hat + sin_a * t_hat                          # photon direction

        beta_vec = beta * v_hat
        gamma = 1.0 / np.sqrt(1.0 - beta**2)
        e_prime = gamma * (1.0 - beta_vec @ k0)                     # comoving energy (E=1)
        # Spatial part of the boosted null 4-vector, then normalize by E' to a unit dir.
        b2 = beta**2
        k_prime = k0 + ((gamma - 1.0) * (beta_vec @ k0) / b2 - gamma) * beta_vec
        k_prime_hat = k_prime / e_prime
        deltas.append(1.0 / e_prime)                               # δ = E/E' = 1/E'
        mus.append(float(n_hat @ k_prime_hat))                     # cos α'
    return np.array(deltas), np.array(mus)


@pytest.mark.parametrize("i_deg, th_deg", [(90.0, 90.0), (30.0, 20.0), (60.0, 110.0), (80.0, 45.0)])
def test_doppler_and_aberration_match_explicit_boost(i_deg, th_deg):
    """cos_xi/doppler_factor and μ' = δ cos α reproduce the explicit Lorentz boost.

    The module's closed-form trig must equal a from-scratch 4-vector boost to machine
    precision — for a spread of geometries, including grazing/eclipsing colatitudes.
    """
    i, theta, u = np.deg2rad(i_deg), np.deg2rad(th_deg), 0.4
    phase = np.linspace(0.0, 2.0 * np.pi, 37, endpoint=False)
    beta = spot_speed(500.0, 12.0, theta, u)

    cps = cos_psi(phase, i, theta)
    cosa = bend(cps, u)
    delta_mod = doppler_factor(beta, cos_xi(phase, cosa, cps, i))
    mu_mod = delta_mod * cosa                                       # μ' the flux uses

    delta_ref, mu_ref = _independent_delta_and_mu_prime(phase, i, theta, u, beta)
    assert np.allclose(delta_mod, delta_ref, atol=1e-12)
    assert np.allclose(mu_mod, mu_ref, atol=1e-12)


def test_aberrated_cosine_stays_physical():
    """μ' = δ cos α ∈ [−1, 1] for all phases (aberration is a real Lorentz map).

    In particular μ' never exceeds 1 even at the blueshifted phases, where δ > 1 but
    cos α is correspondingly small — the guard that lets the beaming lookup stay valid.
    """
    i, theta, u = np.deg2rad(90.0), np.deg2rad(90.0), 0.3
    phase = np.linspace(0.0, 2.0 * np.pi, 512, endpoint=False)
    beta = spot_speed(600.0, 12.0, theta, u)
    cps = cos_psi(phase, i, theta)
    cosa = bend(cps, u)
    mu_prime = doppler_factor(beta, cos_xi(phase, cosa, cps, i)) * cosa
    assert mu_prime.max() <= 1.0 + 1e-12
    assert mu_prime.min() >= -1.0 - 1e-12
    # where the spot is visible (cos α ≥ 0), the aberrated cosine is also ≥ 0:
    assert np.all(mu_prime[cosa >= 0.0] >= -1e-12)


# --- ν → 0 bit-for-bit + the layer actually does something --------------------

def test_zero_spin_is_bit_for_bit_frozen_flux():
    """A zero-spin Rotation reproduces the non-rotating flux exactly (array_equal)."""
    g = ALWAYS_VISIBLE
    phase = np.linspace(0.0, 2.0 * np.pi, 128, endpoint=False)
    frozen = point_spot_flux(phase, g["inclination"], g["colatitude"], g["compactness"])
    zero_spin = point_spot_flux(phase, g["inclination"], g["colatitude"], g["compactness"],
                                rotation=Rotation(spin_hz=0.0, radius_km=12.0))
    assert np.array_equal(frozen, zero_spin)


def test_zero_spin_bit_for_bit_through_compute_profile_and_beaming():
    """The bit-for-bit identity also holds through compute_profile with a beaming law."""
    from mcrt import eddington_limb_darkening
    g = ALWAYS_VISIBLE
    frozen = compute_profile(g["inclination"], g["colatitude"], g["compactness"],
                             beaming=eddington_limb_darkening)
    zero = compute_profile(g["inclination"], g["colatitude"], g["compactness"],
                           beaming=eddington_limb_darkening,
                           rotation=Rotation(spin_hz=0.0, radius_km=12.0))
    assert np.array_equal(frozen.flux, zero.flux)
    assert np.array_equal(frozen.cos_alpha, zero.cos_alpha)


def test_finite_spin_actually_changes_the_flux():
    """At a real spin the flux departs from the frozen one — the layer is not a no-op."""
    g = ALWAYS_VISIBLE
    phase = np.linspace(0.0, 2.0 * np.pi, 256, endpoint=False)
    frozen = point_spot_flux(phase, g["inclination"], g["colatitude"], g["compactness"])
    spun = point_spot_flux(phase, g["inclination"], g["colatitude"], g["compactness"],
                           rotation=Rotation(spin_hz=J0740_SPIN_HZ, radius_km=J0740_R_KM))
    assert not np.allclose(frozen, spun)


def test_small_spin_limit_approaches_frozen():
    """As ν → 0 the rotating flux converges to the frozen one (continuity of the layer)."""
    g = ALWAYS_VISIBLE
    phase = np.linspace(0.0, 2.0 * np.pi, 128, endpoint=False)
    frozen = point_spot_flux(phase, g["inclination"], g["colatitude"], g["compactness"])
    frozen_n = frozen / frozen.max()
    prev = np.inf
    for nu in (100.0, 10.0, 1.0):
        spun = point_spot_flux(phase, g["inclination"], g["colatitude"], g["compactness"],
                               rotation=Rotation(spin_hz=nu, radius_km=12.0))
        dev = float(np.abs(spun / spun.max() - frozen_n).max())
        assert dev < prev            # monotone shrink toward the frozen shape
        prev = dev
    assert prev < 5e-4               # at 1 Hz the departure is negligible (≪ 1% bending floor)


# --- rotation breaks the fore-aft symmetry (the Doppler signature) -----------

def test_rotation_breaks_phase_symmetry_that_frozen_flux_preserves():
    """Frozen flux is even in φ; rotation makes it asymmetric (approaching ≠ receding)."""
    g = ALWAYS_VISIBLE
    n = 512
    frozen = compute_profile(g["inclination"], g["colatitude"], g["compactness"], n_phase=n).flux
    spun = compute_profile(g["inclination"], g["colatitude"], g["compactness"], n_phase=n,
                           rotation=Rotation(spin_hz=J0740_SPIN_HZ, radius_km=J0740_R_KM)).flux
    # φ_k ↔ −φ_k on the periodic grid maps index k → (n − k) mod n.
    mirror = (-np.arange(n)) % n
    assert np.allclose(frozen, frozen[mirror])        # frozen: even in φ
    assert not np.allclose(spun, spun[mirror])        # rotating: fore-aft asymmetric


# --- the δ⁴ energy vs δ³ photon flux flag ------------------------------------

def test_photon_flux_flag_removes_one_doppler_power():
    """Energy (δ⁴) and photon (δ³) fluxes differ by exactly one factor of δ per phase.

    Their pointwise ratio (where visible) must equal the Doppler factor δ(φ).
    """
    i, theta, u = np.deg2rad(90.0), np.deg2rad(90.0), 0.3
    phase = np.linspace(0.0, 2.0 * np.pi, 256, endpoint=False)
    energy = point_spot_flux(phase, i, theta, u,
                             rotation=Rotation(J0740_SPIN_HZ, J0740_R_KM, photon_flux=False))
    photon = point_spot_flux(phase, i, theta, u,
                             rotation=Rotation(J0740_SPIN_HZ, J0740_R_KM, photon_flux=True))
    cps = cos_psi(phase, i, theta)
    cosa = bend(cps, u)
    delta = doppler_factor(spot_speed(J0740_SPIN_HZ, J0740_R_KM, theta, u),
                           cos_xi(phase, cosa, cps, i))
    vis = cosa > 1e-6
    assert np.allclose(energy[vis] / photon[vis], delta[vis])
    assert Rotation(1.0, 1.0).flux_exponent == 4
    assert Rotation(1.0, 1.0, photon_flux=True).flux_exponent == 3


# --- band-limited boost (BandSpectrum / band_boost, Track C3) -----------------

# J0740 Riley regime: kT from log₁₀T = 5.99 (Riley 2021 ST-U), NICER band
# 0.3–1.5 keV (channels [30, 150)). √(1−u)·kT ≈ 0.060 keV, so the band sits at
# E/kT_obs ≈ 5–25 — deep on the Wien tail, the exponential-sensitivity regime.
RILEY_BAND = BandSpectrum(kt_kev=0.0842, e_min_kev=0.3, e_max_kev=1.5)


def test_band_boost_is_unity_at_delta_one():
    """δ = 1 must give weight exactly 1 in both conventions (the normalization)."""
    assert band_boost(1.0, J0740_U, RILEY_BAND, photon_flux=True) == pytest.approx(1.0)
    arr = band_boost(np.ones(5), J0740_U, RILEY_BAND, photon_flux=False)
    assert np.allclose(arr, 1.0)


def test_band_boost_wide_band_recovers_bolometric_powers():
    """A band covering (essentially) the whole spectrum reduces to δ³ / δ⁴.

    Φ_k(δ) ∝ (δ√(1−u)kT)^{k+1} as [E1, E2] → (0, ∞), so the normalized weight is a
    pure power: photon (k = 2) → δ³, energy (k = 3) → δ⁴ — i.e. the band machinery
    *contains* both bolometric conventions as limits, pinning it to the verified δⁿ.
    """
    wide = BandSpectrum(kt_kev=0.0842, e_min_kev=1e-4, e_max_kev=8.0, n_quad=512)
    for delta in (0.88, 0.95, 1.05, 1.14):
        assert band_boost(delta, J0740_U, wide, photon_flux=True) == pytest.approx(
            delta**3, rel=1e-3)
        assert band_boost(delta, J0740_U, wide, photon_flux=False) == pytest.approx(
            delta**4, rel=1e-3)


def test_band_boost_wien_tail_is_steeper_than_bolometric():
    """On the J0740 NICER band (Wien tail) the weight beats δ⁴ in both directions.

    Flux there goes like exp(−E/δkT_obs): blueshifted phases gain *more* than δ⁴,
    redshifted phases lose more — the deep-dive §8.2 direction argument, as a test.
    """
    for delta in (1.05, 1.12):
        assert band_boost(delta, J0740_U, RILEY_BAND) > delta**4
    for delta in (0.88, 0.95):
        assert band_boost(delta, J0740_U, RILEY_BAND) < delta**4


def test_band_boost_monotone_in_delta():
    """More blueshift, more in-band flux: Φ(δ)/Φ(1) strictly increasing in δ."""
    deltas = np.linspace(0.85, 1.15, 31)
    assert np.all(np.diff(band_boost(deltas, J0740_U, RILEY_BAND)) > 0.0)


def test_band_boost_quadrature_converged():
    """64 Gauss–Legendre nodes match 512 to ≪1e-8 across the physical δ range."""
    fine = BandSpectrum(kt_kev=0.0842, e_min_kev=0.3, e_max_kev=1.5, n_quad=512)
    d = np.linspace(0.85, 1.15, 7)
    assert np.allclose(band_boost(d, J0740_U, RILEY_BAND),
                       band_boost(d, J0740_U, fine), rtol=1e-10)


def test_band_zero_spin_is_bit_for_bit_frozen_flux():
    """ν → 0 stays bit-for-bit frozen with a band attached (δ ≡ 1 ⇒ weight ≡ 1)."""
    g = ALWAYS_VISIBLE
    phase = np.linspace(0.0, 2.0 * np.pi, 128, endpoint=False)
    frozen = point_spot_flux(phase, g["inclination"], g["colatitude"], g["compactness"])
    banded = point_spot_flux(phase, g["inclination"], g["colatitude"], g["compactness"],
                             rotation=Rotation(spin_hz=0.0, radius_km=12.0,
                                               photon_flux=True, band=RILEY_BAND))
    assert np.array_equal(frozen, banded)


def test_band_spectrum_validates_inputs():
    """Non-positive kT or an inverted/empty band raises rather than mis-integrating."""
    with pytest.raises(ValueError):
        BandSpectrum(kt_kev=-0.1, e_min_kev=0.3, e_max_kev=1.5)
    with pytest.raises(ValueError):
        BandSpectrum(kt_kev=0.1, e_min_kev=1.5, e_max_kev=0.3)


# --- travel-time delay utility (Bogdanov eq. 18; driver-level, C4) ------------

def test_travel_time_delay_flat_space_closed_form():
    """u → 0: Δt = (R/c)(1 − cos α) exactly — elementary flat-space geometry.

    (The eq. 18 integral evaluates in closed form at u = 0; pins the quadrature to
    something derived with no shared code.)
    """
    cosa = np.linspace(0.05, 1.0, 12)
    delay = travel_time_delay(cosa, 0.0, 12.0)
    expected = (12.0 / SPEED_OF_LIGHT_KM_S) * (1.0 - cosa)
    assert np.allclose(delay, expected, rtol=1e-6)


def test_travel_time_delay_radial_zero_monotone_eclipse_masked():
    """Δt(α=0) = 0 exactly; Δt grows toward grazing; eclipsed samples get 0."""
    cosa = np.linspace(1.0, 0.0, 21)
    delay = travel_time_delay(cosa, J0740_U, J0740_R_KM)
    assert delay[0] == 0.0
    assert np.all(np.diff(delay) > 0.0)          # larger α (smaller μ) → longer path
    assert np.all(travel_time_delay(np.array([-0.3, -0.01]), J0740_U, J0740_R_KM) == 0.0)


def test_travel_time_delay_j0740_warp_is_percent_scale():
    """At J0740 the grazing warp ν·Δt is a few % of a cycle — the stated caveat scale.

    (≈ 0.036 cyc at cos α = 0: the u = 0.494 compactness roughly doubles the flat-space
    R/c scale for grazing paths. Guards order of magnitude, not the exact value.)
    """
    warp = J0740_SPIN_HZ * travel_time_delay(np.array([0.0]), J0740_U, J0740_R_KM)[0]
    assert 0.01 < warp < 0.06


# --- external anchor: SD1c 200 Hz code-comparison waveform -------------------

# Bogdanov et al. 2019 (ApJL 887 L26) Table 1, test SD1c: point spot, i = θ_s = 90°,
# M = 1.4 M_sun, R = 12 km, ν = 200 Hz, blackbody kT = 0.35 keV (comoving), isotropic,
# monochromatic photon flux at 1 keV. Same geometry as SD1a (used in code_comparison.py)
# but rotating — so it tests exactly the Doppler + aberration layer. See
# scripts/c1_doppler_validate.py for the full driver and figure.
SD1C_REFERENCE = Path(__file__).resolve().parents[1] / "data" / "l26_reference" / "SD1c_test_IM.txt"
SD1_M_SUN_KM, SD1_M, SD1_R, SD1_KT, SD1_EOBS = 1.47662, 1.4, 12.0, 0.35, 1.0
SD1_U = 2.0 * SD1_M_SUN_KM * SD1_M / SD1_R
SD1C_TOLERANCE = 0.015   # ~1% Beloborodov-bending floor (as SD1a) + a little headroom


def _sd1c_flux(phase_cycles):
    """SD1c 1 keV photon-flux shape via the collapsed eq. (20) + travel-time warp.

    Independent inline re-implementation (uses only mcrt primitives), so agreement
    with the IM reference is a genuine end-to-end check of the Doppler layer.
    """
    incl = colat = np.deg2rad(90.0)
    n = 4096
    phi = np.linspace(0.0, 2.0 * np.pi, n, endpoint=False)
    cps = cos_psi(phi, incl, colat)
    cosa = bend(cps, SD1_U)
    delta = doppler_factor(spot_speed(200.0, SD1_R, colat, SD1_U), cos_xi(phi, cosa, cps, incl))
    e_prime = SD1_EOBS / (delta * np.sqrt(1.0 - SD1_U))
    flux = np.where(cosa >= 0.0, (delta * cosa) * (1.0 - SD1_U) / np.expm1(e_prime / SD1_KT), 0.0)

    # Light-travel-time delay Δφ[cyc] = ν Δt(α) (Bogdanov eq. 18), applied to match the
    # reference; geometry-only, so it is excluded from the production layer.
    sin_a = np.sqrt(np.clip(1.0 - cosa**2, 0.0, 1.0))
    beta_b = sin_a / np.sqrt(1.0 - SD1_U)
    t, w = np.polynomial.legendre.leggauss(512)
    x = 0.5 * (t + 1.0)
    one_minus_ux = 1.0 - SD1_U * x
    integrand = ((1.0 / np.sqrt(np.clip(1.0 - beta_b[..., None] ** 2 * x**2 * one_minus_ux,
                                        1e-300, None))) - 1.0) / (x**2 * one_minus_ux)
    dt = (SD1_R * 1.0e3 / 299792458.0) * (integrand * (0.5 * w)).sum(axis=-1)
    dt = np.where(cosa >= 0.0, dt, 0.0)
    obs = (phi / (2.0 * np.pi) + 200.0 * dt) % 1.0

    order = np.argsort(obs)
    xs, ys = obs[order], flux[order]
    xs_ext = np.concatenate([xs - 1.0, xs, xs + 1.0])
    ys_ext = np.concatenate([ys, ys, ys])
    return np.interp(phase_cycles % 1.0, xs_ext, ys_ext)


@pytest.mark.skipif(not SD1C_REFERENCE.exists(),
                    reason="L26 supplementary reference (data/l26_reference/SD1c_test_IM.txt) not present")
def test_sd1c_rotating_waveform_matches_im_reference():
    """Our 200 Hz Doppler waveform reproduces the IM SD1c reference to ~1%.

    The same geometry non-rotating (SD1a) already matches to ~0.8% (code_comparison.py);
    this shows the Doppler + aberration layer carries that agreement to 200 Hz, where the
    pulse is visibly skewed. Residual is the linear-bending approximation, worst at the
    grazing eclipse edge — identical in character to SD1a.
    """
    ref = np.loadtxt(SD1C_REFERENCE)
    phase, f_ref = ref[:, 0], ref[:, 1]
    ours = _sd1c_flux(phase)
    ours_n, ref_n = ours / ours.max(), f_ref / f_ref.max()

    assert np.abs(ours_n - ref_n).max() < SD1C_TOLERANCE

    # The reference is genuinely rotating: it carries a fore-aft asymmetry (~14% of
    # peak) that a static pulse would not, and ours reproduces it.
    def fore_aft(curve):
        return float(np.sqrt(np.mean((curve - curve[(-np.arange(curve.size)) % curve.size]) ** 2)))
    assert fore_aft(ref_n) > 0.05
    assert fore_aft(ours_n) == pytest.approx(fore_aft(ref_n), abs=0.02)
