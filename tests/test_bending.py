"""Validation tests for exact Schwarzschild bending (`mcrt.bending`).

The load-bearing checks are *external references* the module must reproduce, not
internal self-consistency: the flat-space limit ψ(α) = α (algebraically exact),
an independent high-resolution quadrature of the same integral, and an
independently-derived analytic Jacobian. On top of those, the tests pin the
relationship to the linear Beloborodov map the exact integral is meant to
supersede — agreement to ~1% at moderate compactness, a gap that grows toward the
grazing angles where J0740's faint phase lives — and confirm the linear default
path in `mcrt.pulse` is left bit-for-bit unchanged.
"""

from pathlib import Path

import numpy as np
import pytest

from mcrt.bending import (
    DEFAULT_N_QUAD,
    ExactBending,
    bend_exact,
    deflection_angle,
)
from mcrt.pulse import (
    bend,
    compute_profile,
    cos_psi,
    point_spot_flux,
    pulsed_fraction,
    visibility_threshold,
)
from mcrt.theory import eddington_limb_darkening

# Neutron-star compactnesses that anchor the whole project (see test_pulse.py):
# SD1a code-comparison (u ≈ 0.3445), J0740 Riley (0.494) and Miller (0.444).
SD1A_U = 0.3445
J0740_U = 0.494


# --- the deflection integral against external references -----------------------

def test_deflection_reduces_to_flat_space_at_zero_compactness():
    """u → 0: ψ(α) = α exactly (no bending). The integral becomes arcsin(sin α)."""
    alpha = np.linspace(0.0, 1.5, 40)          # up to ~86°; the grazing endpoint
    psi = deflection_angle(alpha, 1e-6)        # itself is an integrable singularity
    assert np.allclose(psi, alpha, atol=1e-5)


def test_deflection_matches_independent_high_resolution_quadrature():
    """Gauss–Legendre ψ(α) agrees with a fine, independent trapezoid of the same
    integrand — certifies the quadrature value, not just its internal convergence.
    """
    u = J0740_U
    for alpha in (0.5, 1.0, 1.4):
        beta = np.sin(alpha) / np.sqrt(1.0 - u)
        x = np.linspace(0.0, 1.0, 4_000_001)[1:-1]   # exclude the singular endpoints
        trap = np.trapezoid(beta / np.sqrt(1.0 - beta**2 * x**2 * (1.0 - u * x)), x)
        gl = deflection_angle(np.array([alpha]), u)[0]
        assert gl == pytest.approx(trap, abs=1e-4)


def test_deflection_quadrature_is_converged():
    """Doubling the node count leaves ψ unchanged to ~1e-8 — the node default is safe."""
    u = J0740_U
    alpha = np.array([0.3, 0.9, 1.4])
    coarse = deflection_angle(alpha, u, n_quad=DEFAULT_N_QUAD)
    fine = deflection_angle(alpha, u, n_quad=2 * DEFAULT_N_QUAD)
    assert np.allclose(coarse, fine, atol=1e-8)


def test_deflection_exceeds_emission_angle_and_is_monotone():
    """Gravity bends rays outward (ψ ≥ α) and ψ(α) increases — the invertibility premise."""
    u = J0740_U
    alpha = np.linspace(0.0, 0.5 * np.pi, 200)
    psi = deflection_angle(alpha, u)
    assert np.all(psi >= alpha - 1e-9)
    assert np.all(np.diff(psi) > 0.0)


# --- the inverted map: cos α(cos ψ) -------------------------------------------

def test_map_is_identity_at_zero_compactness():
    """With no bending the map is cos α = cos ψ (and a facing spot always emits along n).

    Sampled on the non-grazing interior (cos ψ ≥ 0.05): the exact grazing horizon at
    ψ = π/2 sits on the integrand's √-singularity, where fixed quadrature is
    ~1e-3-limited — a zero-flux edge, and the tight ψ = α check lives in
    ``test_deflection_reduces_to_flat_space_at_zero_compactness``.
    """
    eb = ExactBending(1e-6)
    cps = np.linspace(0.05, 1.0, 50)
    assert np.allclose(eb.cos_alpha(cps), cps, atol=1e-4)
    assert eb.cos_alpha(1.0) == pytest.approx(1.0, abs=1e-6)


def test_map_is_monotone_increasing_in_cos_psi():
    """cos α rises with cos ψ across the whole visible range (needed for a single-valued map)."""
    eb = ExactBending(J0740_U)
    cps = np.linspace(eb.cos_psi_visible, 1.0, 500)
    ca = eb.cos_alpha(cps)
    assert np.all(np.diff(ca) >= -1e-12)


def test_eclipsed_cos_psi_maps_below_zero():
    """Below the visibility threshold the map returns cos α < 0, so the spot reads invisible."""
    eb = ExactBending(J0740_U)
    below = eb.cos_psi_visible - 0.05
    assert eb.cos_alpha(below) < 0.0
    # and exactly at the threshold it grazes (cos α ≈ 0).
    assert eb.cos_alpha(eb.cos_psi_visible) == pytest.approx(0.0, abs=1e-6)


def test_bend_exact_convenience_matches_the_class():
    """The module-level convenience wrapper equals building the class and evaluating it."""
    cps = np.linspace(-0.5, 1.0, 20)
    assert np.allclose(bend_exact(cps, SD1A_U), ExactBending(SD1A_U).cos_alpha(cps))


# --- the lensing Jacobian D(ψ) = d cos α / d cos ψ -----------------------------

def _analytic_jacobian(alpha, u, n_quad=4000):
    """D(α) = sin α / (sin ψ · dψ/dα) with dψ/dα derived in closed form:

    differentiating ψ = ∫ β(1 − c)^(−1/2) dx  (c = β² x² (1 − u x), β = sin α/√(1−u))
    gives ∂ψ/∂β = ∫ (1 − c)^(−3/2) dx, and dβ/dα = cos α/√(1−u). Independent of the
    module's np.gradient path, so it validates the stored Jacobian.
    """
    t, w = np.polynomial.legendre.leggauss(n_quad)
    x = 0.5 * (t + 1.0)
    weights = 0.5 * w
    beta = np.sin(alpha) / np.sqrt(1.0 - u)
    c = beta**2 * x**2 * (1.0 - u * x)
    dpsi_dalpha = (np.cos(alpha) / np.sqrt(1.0 - u)) * ((1.0 - c) ** -1.5 * weights).sum()
    psi = deflection_angle(np.array([alpha]), u, n_quad=n_quad)[0]
    return np.sin(alpha) / (np.sin(psi) * dpsi_dalpha)


def test_jacobian_at_facing_spot_equals_one_minus_u():
    """At cos ψ = 1 (α = 0) the exact Jacobian recovers the linear map's constant (1 − u)."""
    for u in (SD1A_U, 0.444, J0740_U):
        assert ExactBending(u).jacobian(1.0) == pytest.approx(1.0 - u, abs=1e-4)


def test_jacobian_matches_independent_analytic_derivative():
    """The stored D(cos ψ) reproduces the closed-form d cos α/d cos ψ, including the
    grazing region (cos α = 0.005) where J0740's faint phase reads the library.
    """
    u = J0740_U
    eb = ExactBending(u)
    for cos_a in (0.005, 0.01, 0.05, 0.2, 0.5, 0.9):
        alpha = np.arccos(cos_a)
        cps = float(np.cos(deflection_angle(np.array([alpha]), u)[0]))
        assert eb.jacobian(cps) == pytest.approx(_analytic_jacobian(alpha, u), abs=2e-3)


# --- relationship to the linear Beloborodov map --------------------------------

def test_agrees_with_beloborodov_to_one_percent_at_sd1a_compactness():
    """At the validated SD1a compactness (u ≈ 0.34) the linear map is within ~1% in cos α."""
    eb = ExactBending(SD1A_U)
    cps = np.linspace(eb.cos_psi_visible, 1.0, 400)
    max_gap = float(np.max(np.abs(eb.cos_alpha(cps) - bend(cps, SD1A_U))))
    assert max_gap < 0.015


def test_linear_map_error_grows_with_compactness():
    """The exact–linear gap is larger at J0740's u = 0.494 than at SD1a's 0.34 — which is
    why the exact re-check matters precisely for the most compact anchor.
    """
    def max_gap(u):
        eb = ExactBending(u)
        cps = np.linspace(eb.cos_psi_visible, 1.0, 400)
        return float(np.max(np.abs(eb.cos_alpha(cps) - bend(cps, u))))

    assert max_gap(J0740_U) > max_gap(SD1A_U) > 0.0


def test_exact_visibility_threshold_differs_from_linear_but_stays_around_the_back():
    """Exact bending still lets the spot be seen past ψ = 90° (cos ψ < 0), but the
    grazing threshold is not the linear −u/(1−u): a measurable geometric difference.
    """
    eb = ExactBending(J0740_U)
    assert eb.cos_psi_visible < 0.0                                   # sees around the back
    assert eb.cos_psi_visible != pytest.approx(visibility_threshold(J0740_U), abs=1e-3)


# --- integration with pulse.py: default untouched, exact opt-in works ----------

def test_linear_default_flux_is_unchanged_bit_for_bit():
    """point_spot_flux with bending=None reproduces the pre-existing (1−u)·I·cos α path."""
    i, th, u = np.deg2rad(30.0), np.deg2rad(20.0), 0.3
    phase = np.linspace(0.0, 2.0 * np.pi, 256, endpoint=False)
    flux = point_spot_flux(phase, i, th, u)
    cps = cos_psi(phase, i, th)
    cos_a = bend(cps, u)
    expected = np.where(cos_a >= 0.0, (1.0 - u) * cos_a, 0.0)
    assert np.array_equal(flux, expected)


def test_compactness_mismatch_between_flux_and_bending_map_is_rejected():
    """A bending map built for the wrong u must not silently corrupt the geometry."""
    with pytest.raises(ValueError):
        point_spot_flux(np.array([0.0]), 1.5, 1.35, J0740_U,
                        bending=ExactBending(0.3))


def test_exact_bending_changes_flux_but_keeps_a_valid_profile():
    """Swapping in the exact map perturbs the waveform yet leaves cos α and visibility sane."""
    i, th, u = 1.5284, 1.35, J0740_U
    lin = compute_profile(i, th, u, n_phase=1024)
    exact = compute_profile(i, th, u, n_phase=1024, bending=ExactBending(u))

    assert not np.allclose(lin.flux, exact.flux)                # the map genuinely moved
    assert np.all(exact.cos_alpha[exact.visible] >= 0.0)        # visible ⇒ cos α ≥ 0
    assert exact.flux.min() >= 0.0                              # no negative flux leaks


def test_exact_bending_shifts_j0740_delta_pf_but_keeps_it_live():
    """The measured B-track result in miniature: exact bending moves J0740's two-spot ΔPF
    by more than the seed error (~0.003) yet the systematic stays live (ΔPF > 0.1).

    Uses the Eddington limb-darkening stand-in so the claim rests on the geometry, not on
    the production library numbers; Track B2 repeats it against the real I(μ; τ).
    """
    incl, u = 1.5284, J0740_U
    spots = [(1.35, 0.0), (1.89, 0.442)]
    n = 1024

    def two_spot_delta_pf(bending):
        def flux(beaming):
            total = np.zeros(n)
            for colat, azim in spots:
                base = compute_profile(incl, colat, u, beaming=beaming,
                                       n_phase=n, bending=bending).flux
                total += np.roll(base, int(round(azim * n)) % n)
            return total
        return pulsed_fraction(flux(eddington_limb_darkening)) - pulsed_fraction(flux(None))

    d_linear = two_spot_delta_pf(None)
    d_exact = two_spot_delta_pf(ExactBending(u))
    assert d_exact > 0.1                        # still a live PF systematic
    assert abs(d_exact - d_linear) > 0.003      # but shifted beyond seed error


# --- the SD1a self-check: exact bending shrinks the code-comparison residual ----

SD1A_REFERENCE = (Path(__file__).resolve().parents[1]
                  / "data" / "l26_reference" / "SD1a_test_IM.txt")


@pytest.mark.skipif(not SD1A_REFERENCE.exists(),
                    reason="L26 supplementary reference (data/l26_reference/SD1a_test_IM.txt) not present")
def test_exact_bending_shrinks_sd1a_residual_vs_im_reference():
    """Against the independent IM waveform, exact bending should not be worse than the
    linear map — the "free self-check" of Track B1 (the ~1% gap comes from bending).
    """
    ref = np.loadtxt(SD1A_REFERENCE)
    phase_cycles, f_ref = ref[:, 0], ref[:, 1]
    phi = 2.0 * np.pi * phase_cycles
    i = th = np.deg2rad(90.0)

    ref_n = f_ref / f_ref.max()
    lin = point_spot_flux(phi, i, th, SD1A_U)
    exact = point_spot_flux(phi, i, th, SD1A_U, bending=ExactBending(SD1A_U))
    lin_res = np.abs(lin / lin.max() - ref_n).max()
    exact_res = np.abs(exact / exact.max() - ref_n).max()

    assert exact_res <= lin_res + 1e-3          # exact is no worse (and generally better)


# --- constructor validation ----------------------------------------------------

@pytest.mark.parametrize("bad_u", [-0.1, 2.0 / 3.0, 0.7, 1.0, 1.5])
def test_constructor_rejects_out_of_range_compactness(bad_u):
    """Rejects negatives and any u ≥ 2/3 (surface at/inside the photon sphere, where
    the deflection integral diverges) — with a physics message, not a resolution one.
    """
    with pytest.raises(ValueError):
        ExactBending(bad_u)
