"""Rung A verification for the point-spot pulse-profile machinery.

The load-bearing test is ``test_rung_a_*``: the full numerical pipeline
(phase grid → cos ψ → Beloborodov bending → visibility clip → flux → pulsed
fraction) must reproduce an *algebraically independent* closed form for the
isotropic point spot. The closed forms are derived inline here from the
emission-angle extremes, so they share no code with the module under test.
"""

import numpy as np
import pytest

from mcrt.pulse import (
    analytic_isotropic_pf,
    bend,
    compute_profile,
    cos_psi,
    point_spot_flux,
    pulsed_fraction,
    visibility_threshold,
)

# A geometry the spot never dips behind the star: cos ψ_min = cos(i+θ_s) stays
# above the bending visibility threshold −u/(1−u). Verified in the test body.
ALWAYS_VISIBLE = dict(inclination=np.deg2rad(30.0), colatitude=np.deg2rad(20.0), compactness=0.3)

# A geometry with ψ_max > 90°: the spot is eclipsed near φ = π but is partly
# visible "around the back" (cos ψ < 0 yet cos α ≥ 0) thanks to light bending.
ECLIPSING = dict(inclination=np.deg2rad(80.0), colatitude=np.deg2rad(80.0), compactness=0.2)


# --- geometry helpers --------------------------------------------------------

def test_cos_psi_extremes_are_at_phase_0_and_pi():
    """cos ψ(0) = cos(i − θ_s) (spot nearest), cos ψ(π) = cos(i + θ_s) (farthest)."""
    i, theta = np.deg2rad(60.0), np.deg2rad(20.0)
    assert cos_psi(0.0, i, theta) == pytest.approx(np.cos(i - theta))
    assert cos_psi(np.pi, i, theta) == pytest.approx(np.cos(i + theta))


def test_bend_is_identity_on_a_pole_on_spot():
    """A spot facing us (cos ψ = 1) emits along the normal (cos α = 1) for any u."""
    for u in (0.0, 0.2, 0.4):
        assert bend(1.0, u) == pytest.approx(1.0)


def test_bend_jacobian_is_constant_one_minus_u():
    """d(cos α)/d(cos ψ) = (1 − u): the clean constant that makes the map cheap."""
    u = 0.35
    psi_grid = np.linspace(-1.0, 1.0, 50)
    d = np.diff(bend(psi_grid, u)) / np.diff(psi_grid)
    assert np.allclose(d, 1.0 - u)


def test_visibility_threshold_admits_seeing_around_the_back():
    """Threshold is −u/(1−u) < 0, so cos ψ < 0 (ψ > 90°) can still be visible."""
    u = 0.3
    assert visibility_threshold(u) == pytest.approx(-u / (1.0 - u))
    assert visibility_threshold(u) < 0.0


# --- pulsed_fraction ----------------------------------------------------------

def test_pulsed_fraction_basic_ratio():
    assert pulsed_fraction(np.array([1.0, 3.0])) == pytest.approx(0.5)  # (3-1)/(3+1)


def test_pulsed_fraction_zero_when_flux_identically_zero():
    """An invisible spot (all-zero flux) yields PF = 0, not a divide-by-zero."""
    assert pulsed_fraction(np.zeros(8)) == 0.0


# --- Rung A: isotropic point spot vs. closed form ----------------------------

def test_rung_a_flux_shape_matches_analytic_isotropic_form():
    """Where visible, isotropic flux is exactly F ∝ (1−u)(u + (1−u) cos ψ).

    This pins the flux *pipeline* (bending + projection) to machine precision,
    independent of any pulsed-fraction bookkeeping.
    """
    g = ALWAYS_VISIBLE
    phase = np.linspace(0.0, 2.0 * np.pi, 256, endpoint=False)
    u = g["compactness"]
    cps = cos_psi(phase, g["inclination"], g["colatitude"])
    expected = (1.0 - u) * (u + (1.0 - u) * cps)  # isotropic: I ≡ 1, cos α = u+(1−u)cosψ

    flux = point_spot_flux(phase, g["inclination"], g["colatitude"], u)  # beaming=None → isotropic
    assert np.allclose(flux, expected)


def test_rung_a_pulsed_fraction_matches_closed_form_always_visible():
    """Numerical PF reproduces the closed-form isotropic PF to < 1% (here far tighter).

    Closed form (always visible, so extremes sit at φ = 0, π):
        PF = (1−u)(cosψ_max − cosψ_min) / (2u + (1−u)(cosψ_max + cosψ_min)).
    Computed inline from cos(i∓θ_s) — no shared code with the module.
    """
    g = ALWAYS_VISIBLE
    i, theta, u = g["inclination"], g["colatitude"], g["compactness"]
    cpsi_max, cpsi_min = np.cos(i - theta), np.cos(i + theta)
    assert cpsi_min >= visibility_threshold(u)  # confirm the spot truly never sets

    pf_closed = ((1.0 - u) * (cpsi_max - cpsi_min)) / (2.0 * u + (1.0 - u) * (cpsi_max + cpsi_min))

    prof = compute_profile(i, theta, u, n_phase=1024)
    pf_numeric = pulsed_fraction(prof.flux)

    assert pf_numeric == pytest.approx(pf_closed, rel=1e-3)            # the < 1% acceptance target
    assert analytic_isotropic_pf(i, theta, u) == pytest.approx(pf_closed)  # module closed form agrees too


def test_rung_a_eclipse_drives_flux_to_zero_and_pf_to_one():
    """When ψ_max > 90° past the bending threshold, the spot sets: F_min = 0, PF = 1."""
    g = ECLIPSING
    prof = compute_profile(g["inclination"], g["colatitude"], g["compactness"], n_phase=1024)
    assert prof.flux.min() == 0.0
    assert prof.flux.max() > 0.0
    assert pulsed_fraction(prof.flux) == pytest.approx(1.0)


def test_rung_a_spot_is_visible_around_the_back():
    """Light bending makes part of the far hemisphere (cos ψ < 0) still visible."""
    g = ECLIPSING
    prof = compute_profile(g["inclination"], g["colatitude"], g["compactness"], n_phase=1024)
    cps = cos_psi(prof.phase, g["inclination"], g["colatitude"])
    assert np.any(prof.visible & (cps < 0.0))   # seeing around the star
    assert np.any(~prof.visible)                # but it does fully set somewhere


def test_analytic_isotropic_pf_rejects_eclipsing_geometry():
    """The closed form only holds when the spot never sets; otherwise it must refuse."""
    g = ECLIPSING
    with pytest.raises(ValueError):
        analytic_isotropic_pf(g["inclination"], g["colatitude"], g["compactness"])


# --- beaming swap (the seam Rung C will use) ---------------------------------

def test_beaming_callable_modulates_flux_but_not_geometry():
    """Passing I(μ) changes only the brightness term; isotropic ≡ constant beaming."""
    g = ALWAYS_VISIBLE
    phase = np.linspace(0.0, 2.0 * np.pi, 128, endpoint=False)
    iso = point_spot_flux(phase, g["inclination"], g["colatitude"], g["compactness"])
    const = point_spot_flux(phase, g["inclination"], g["colatitude"], g["compactness"],
                            beaming=lambda mu: np.full_like(np.asarray(mu, float), 1.0))
    assert np.allclose(iso, const)
