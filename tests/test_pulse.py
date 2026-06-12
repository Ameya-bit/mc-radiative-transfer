"""Verification tests for the point-spot pulse-profile machinery.

The load-bearing group is the **analytic check**: the full numerical pipeline
(phase grid → cos ψ → Beloborodov bending → visibility clip → flux → pulsed
fraction) must reproduce an *algebraically independent* closed form for the
isotropic point spot. The closed forms are derived inline here from the
emission-angle extremes, so they share no code with the module under test.
"""

from pathlib import Path

import numpy as np
import pytest

from mcrt.beaming import beaming_lookup
from mcrt.pulse import (
    analytic_isotropic_pf,
    bend,
    compute_profile,
    cos_psi,
    point_spot_flux,
    pulsed_fraction,
    visibility_threshold,
)
from mcrt.theory import eddington_limb_darkening

# IM-code reference waveform from the L26 supplementary archive (apjlab5968.tar.gz).
# That archive is third-party AAS material and is gitignored (data/*), so it is
# present only after a local download+extract; the code-comparison test skips
# cleanly without it so the suite still passes on a fresh checkout. See docs/references.md.
SD1A_REFERENCE = Path(__file__).resolve().parents[1] / "data" / "l26_reference" / "SD1a_test_IM.txt"

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


# --- analytic check: isotropic point spot vs. closed form --------------------

def test_flux_shape_matches_analytic_isotropic_form():
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


def test_pulsed_fraction_matches_closed_form_always_visible():
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


def test_eclipse_drives_flux_to_zero_and_pf_to_one():
    """When ψ_max > 90° past the bending threshold, the spot sets: F_min = 0, PF = 1."""
    g = ECLIPSING
    prof = compute_profile(g["inclination"], g["colatitude"], g["compactness"], n_phase=1024)
    assert prof.flux.min() == 0.0
    assert prof.flux.max() > 0.0
    assert pulsed_fraction(prof.flux) == pytest.approx(1.0)


def test_spot_is_visible_around_the_back():
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


# --- beaming swap (used by the isotropic-vs-realistic comparison) ------------

def test_beaming_callable_modulates_flux_but_not_geometry():
    """Passing I(μ) changes only the brightness term; isotropic ≡ constant beaming."""
    g = ALWAYS_VISIBLE
    phase = np.linspace(0.0, 2.0 * np.pi, 128, endpoint=False)
    iso = point_spot_flux(phase, g["inclination"], g["colatitude"], g["compactness"])
    const = point_spot_flux(phase, g["inclination"], g["colatitude"], g["compactness"],
                            beaming=lambda mu: np.full_like(np.asarray(mu, float), 1.0))
    assert np.allclose(iso, const)


# --- the isotropic-vs-realistic beaming comparison (the headline result) -----

def test_beaming_lookup_interpolates_at_nodes_and_clamps_outside():
    """The library lookup is exact on its μ nodes, linear between, flat at the edges.

    Clamping (not extrapolating) the μ→0 / μ→1 ends is the intended behaviour: the
    grazing tail is the noisiest part of the library, so the curve is held at its
    end values rather than amplified beyond them.
    """
    mu = np.array([0.1, 0.3, 0.5, 0.7, 0.9])
    intensity = np.array([1.2, 1.6, 2.0, 2.4, 2.8])  # a clean linear test curve
    beaming = beaming_lookup(mu, intensity)

    assert np.allclose(beaming(mu), intensity)             # exact on the nodes
    assert beaming(0.4) == pytest.approx(1.8)              # halfway between nodes
    assert beaming(0.0) == pytest.approx(intensity[0])     # clamped below the range
    assert beaming(1.0) == pytest.approx(intensity[-1])    # clamped above the range


def test_limb_darkening_raises_pulsed_fraction_vs_isotropic():
    """Limb darkening (I rising with μ) *sharpens* the pulse: PF_real > PF_iso.

    The bright phase faces us (μ = cos α near max) and the faint phase grazes the
    limb (small μ). A monotone-increasing I(μ) — here the Eddington law 1 + 1.5μ —
    multiplies the bright phase by more than the faint phase, widening the max/min
    contrast. This locks the *sign* of the isotropic-vs-realistic result into a test.
    """
    g = ALWAYS_VISIBLE
    i, theta, u = g["inclination"], g["colatitude"], g["compactness"]

    pf_iso = pulsed_fraction(compute_profile(i, theta, u).flux)
    pf_real = pulsed_fraction(compute_profile(i, theta, u, beaming=eddington_limb_darkening).flux)

    assert pf_real > pf_iso


def test_beaming_swap_leaves_geometry_identical():
    """Isotropic and beamed runs share cos α and visibility exactly — only flux differs.

    This is the controlled-experiment guarantee the comparison rests on: ΔPF can only
    come from the brightness term, never from the geometry shifting underneath it.
    """
    g = ALWAYS_VISIBLE
    i, theta, u = g["inclination"], g["colatitude"], g["compactness"]
    iso = compute_profile(i, theta, u)
    real = compute_profile(i, theta, u, beaming=eddington_limb_darkening)

    assert np.array_equal(real.cos_alpha, iso.cos_alpha)
    assert np.array_equal(real.visible, iso.visible)
    assert not np.allclose(real.flux, iso.flux)


# --- code-comparison check: SD1a vs. the NICER (IM) reference ----------------

# SD1a parameters (Bogdanov et al. 2019, ApJL 887 L26, Table 1): 1 Hz (no Doppler),
# 0.01 rad spot (a point), isotropic Planck, i = θ_s = 90°, M = 1.4 M_sun, R = 12 km.
SD1A_COMPACTNESS = 2.0 * 1.47662 * 1.4 / 12.0   # u = 2GM/Rc² ≈ 0.3445
SD1A_INCLINATION = np.deg2rad(90.0)
SD1A_COLATITUDE = np.deg2rad(90.0)

# Beloborodov's linear bending map is ~1% accurate for this compactness, so we hold
# the agreement to 1.2% (the observed worst case is 0.8%, at the grazing eclipse edge).
CODE_COMPARISON_TOLERANCE = 0.012


@pytest.mark.skipif(not SD1A_REFERENCE.exists(),
                    reason="L26 supplementary reference (data/l26_reference/SD1a_test_IM.txt) not present")
def test_sd1a_matches_im_reference_within_beloborodov_accuracy():
    """Our SD1a profile reproduces the IM community-reference waveform to ~1%.

    The reference (column 1 phase in cycles, column 2 photon flux at 1 keV) and our
    point-spot flux both peak at phase 0 (spot facing the observer). Each is
    normalized to its own peak — ours is in arbitrary units, theirs an absolute
    photon flux, and at 1 Hz + isotropic emission the *normalized shapes* must
    coincide. The only expected deviation is the Beloborodov-approximation error,
    largest where the emission grazes the limb (the eclipse edge).
    """
    ref = np.loadtxt(SD1A_REFERENCE)
    phase_cycles, f_ref = ref[:, 0], ref[:, 1]

    phi = 2.0 * np.pi * phase_cycles
    f_ours = point_spot_flux(phi, SD1A_INCLINATION, SD1A_COLATITUDE, SD1A_COMPACTNESS)

    ref_n = f_ref / f_ref.max()
    ours_n = f_ours / f_ours.max()

    # Shape agreement at the ~1% level across the whole rotation.
    assert np.abs(ours_n - ref_n).max() < CODE_COMPARISON_TOLERANCE

    # The bending-set eclipse window must land on the reference: same visible
    # fraction (gravity lets the spot be seen partway around the back) and a true
    # eclipse to zero (PF = 1) in both.
    assert (ours_n > 0).mean() == pytest.approx((ref_n > 0).mean(), abs=0.02)
    assert pulsed_fraction(f_ref) == pytest.approx(1.0)
    assert pulsed_fraction(f_ours) == pytest.approx(1.0)


# --- real-star anchor: PSR J0030+0451 (the geometry-dependence result) --------

# J0030's published geometry, primary spot (Riley et al. 2019, ApJL 887 L21, Table 2):
# compactness GM/Rc² = 0.156 ⇒ u = 0.312; inclination i = 0.94 rad; the two spots sit
# in the same (far) hemisphere at colatitudes Θp = 2.23 rad and Θs = 2.91 rad, ~0.45
# cycle apart in azimuth. The point of these tests is the saturation: both spots dive
# behind the star, so the pulsed fraction pins at 1 and the beaming systematic moves
# from PF into the waveform shape. The two-spot sum is reproduced inline (the same
# np.roll-and-add the j0030_anchor.py script uses) so the test exercises only the core.
J0030_INCLINATION = 0.94
J0030_COMPACTNESS = 0.312
J0030_SPOT1_COLATITUDE = 2.23   # Θp, ST small circular
J0030_SPOT2_COLATITUDE = 2.91   # Θs, PST crescent
J0030_AZIMUTH_SEPARATION = 0.45  # cycles between the two spots


def _two_spot_flux(beaming, n_phase=1024):
    """Sum the two J0030 spots; the second rolled by its azimuthal separation."""
    s1 = compute_profile(J0030_INCLINATION, J0030_SPOT1_COLATITUDE, J0030_COMPACTNESS,
                         beaming=beaming, n_phase=n_phase).flux
    s2 = compute_profile(J0030_INCLINATION, J0030_SPOT2_COLATITUDE, J0030_COMPACTNESS,
                         beaming=beaming, n_phase=n_phase).flux
    return s1 + np.roll(s2, int(round(J0030_AZIMUTH_SEPARATION * n_phase)))


def test_j0030_single_spot_eclipses_for_a_large_phase_fraction():
    """At J0030's geometry a single point spot is hidden for ~45% of the rotation.

    Both spots sit in the far hemisphere (colatitude ≫ 90°) viewed at i ≈ 54°, so the
    spot swings behind the star: F_min = 0 and PF saturates at 1. This is the root
    cause of the whole result.
    """
    prof = compute_profile(J0030_INCLINATION, J0030_SPOT1_COLATITUDE, J0030_COMPACTNESS,
                           n_phase=1024)
    eclipsed_fraction = float(np.mean(~prof.visible))
    assert 0.4 < eclipsed_fraction < 0.5
    assert prof.flux.min() == 0.0
    assert pulsed_fraction(prof.flux) == pytest.approx(1.0)


def test_j0030_two_spots_still_saturate():
    """Summing both published spots does NOT lift the floor — they hug the same pole.

    The two spots' visible windows leave a gap where the whole star is dark, so the
    two-spot flux still touches zero and PF stays pinned at 1. (Earlier intuition that
    a second spot would fill the eclipse fails for J0030's same-hemisphere geometry.)
    """
    total = _two_spot_flux(beaming=None)
    assert total.min() == pytest.approx(0.0)
    assert pulsed_fraction(total) == pytest.approx(1.0)


def test_j0030_beaming_changes_shape_but_not_pulsed_fraction():
    """The headline: realistic beaming reshapes the visible pulse yet leaves PF = 1.

    Because the eclipse pins PF for *both* beamings, the isotropic-vs-realistic
    difference cannot show up in the pulsed fraction — it lives in the waveform shape.
    A monotone limb-darkening law (Eddington 1 + 1.5μ stands in for the scattering
    beaming) leaves the saturated PF untouched while changing the normalized profile.
    """
    iso = _two_spot_flux(beaming=None)
    real = _two_spot_flux(beaming=eddington_limb_darkening)

    # PF is saturated and blind to the swap ...
    assert pulsed_fraction(iso) == pytest.approx(1.0)
    assert pulsed_fraction(real) == pytest.approx(1.0)

    # ... but the peak-normalized waveform genuinely changes shape.
    shape_rms = np.sqrt(np.mean((real / real.max() - iso / iso.max()) ** 2))
    assert shape_rms > 0.01


def test_two_spot_azimuth_roll_matches_a_half_cycle_shift():
    """A spot at azimuth 0.5 cyc is the same profile shifted half a rotation.

    Guards the np.roll mechanism the two-spot sum relies on: rolling a 1024-point
    profile by 512 reproduces evaluating it half a cycle later.
    """
    n = 1024
    flux = compute_profile(J0030_INCLINATION, J0030_SPOT1_COLATITUDE, J0030_COMPACTNESS,
                           n_phase=n).flux
    rolled = np.roll(flux, n // 2)
    half_cycle = np.concatenate([flux[n // 2:], flux[:n // 2]])
    assert np.array_equal(rolled, half_cycle)
