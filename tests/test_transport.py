"""Track D — exact transport validation against Chandrasekhar H(μ).

Covers the opt-in isotropic-scattering phase function (D1) and certifies that
the transport machinery, run with the phase function for which H(μ) is exact,
reproduces H(μ) within Monte Carlo error (D3). The default Thomson path must
stay bit-for-bit unchanged.
"""

import numpy as np
import pytest

from mcrt import (
    Simulation,
    chandrasekhar_h,
    sample_isotropic_angle,
    sample_thomson_angle,
)

# H(μ)·μ is the emergent photon-flux law of a conservative semi-infinite
# isotropic atmosphere; integrated per bin it is the expected escape
# distribution the transport must reproduce.
BIN_EDGES = np.linspace(0.0, 1.0, 21)
BIN_CENTERS = 0.5 * (BIN_EDGES[:-1] + BIN_EDGES[1:])
_TRAPZ = np.trapezoid if hasattr(np, "trapezoid") else np.trapz


def _h_flux_per_bin():
    return np.array(
        [
            _TRAPZ(chandrasekhar_h(x := np.linspace(lo, hi, 400)) * x, x)
            for lo, hi in zip(BIN_EDGES[:-1], BIN_EDGES[1:])
        ]
    )


# --- D1: the isotropic sampler ---------------------------------------------

def test_isotropic_sampler_in_range():
    rng = np.random.default_rng(0)
    samples = np.array([sample_isotropic_angle(rng) for _ in range(5000)])
    assert np.all(samples >= -1.0) and np.all(samples <= 1.0)


def test_isotropic_sampler_is_uniform():
    """μ uniform on [-1, 1]: mean ≈ 0, variance ≈ 1/3."""
    rng = np.random.default_rng(1)
    samples = np.array([sample_isotropic_angle(rng) for _ in range(200_000)])
    assert np.isclose(samples.mean(), 0.0, atol=0.01)
    assert np.isclose(samples.var(), 1.0 / 3.0, rtol=0.02)


def test_isotropic_flatter_than_thomson():
    """Thomson P(μ)=(3/4)(1+μ²) is peaked at |μ|=1; isotropic is flat.

    A histogram of |μ| is therefore more forward/back-weighted for Thomson.
    """
    rng = np.random.default_rng(2)
    iso = np.abs([sample_isotropic_angle(rng) for _ in range(100_000)])
    thom = np.abs([sample_thomson_angle(rng) for _ in range(100_000)])
    # Thomson puts more weight near |μ|=1.
    assert (thom > 0.8).mean() > (iso > 0.8).mean()


# --- D1: opt-in wiring, default unchanged ----------------------------------

def test_default_phase_function_is_thomson_bit_for_bit():
    """The default engine path is byte-identical to explicit Thomson."""
    a = Simulation(tau_total=5.0, num_photons=800, rng=np.random.default_rng(7))
    a.run()
    b = Simulation(
        tau_total=5.0, num_photons=800, rng=np.random.default_rng(7),
        phase_function="thomson",
    )
    b.run()
    assert np.array_equal(a.results["escaped_mu"], b.results["escaped_mu"])
    assert a.results["absorbed_count"] == b.results["absorbed_count"]


def test_isotropic_differs_from_thomson():
    a = Simulation(
        tau_total=5.0, num_photons=800, rng=np.random.default_rng(7),
        phase_function="thomson",
    )
    a.run()
    c = Simulation(
        tau_total=5.0, num_photons=800, rng=np.random.default_rng(7),
        phase_function="isotropic",
    )
    c.run()
    assert not np.array_equal(a.results["escaped_mu"], c.results["escaped_mu"])


def test_phase_function_attribute_recorded():
    sim = Simulation(tau_total=1.0, num_photons=10, phase_function="isotropic")
    assert sim.phase_function == "isotropic"


def test_invalid_phase_function_raises():
    with pytest.raises(ValueError, match="phase_function must be"):
        Simulation(tau_total=1.0, num_photons=10, phase_function="rayleigh")


# --- D3: isotropic transport reproduces H(μ) -------------------------------

def test_isotropic_transport_matches_chandrasekhar_h():
    """A thick (τ=10, semi-infinite) isotropic slab emits ∝ H(μ)·μ within error.

    Flux-space goodness-of-fit: compare the raw escape distribution to the
    per-bin H(μ)·μ prediction (no I=flux/μ reconstruction, which would carry the
    bin-center division bias). Deterministic seed; the full 5-seed / 4×10⁵-escaped
    driver (scripts/d_isotropic_validate.py) tightens this to reduced χ² = 0.70.
    """
    sim = Simulation(
        tau_total=10.0, num_photons=60_000, rng=np.random.default_rng(0),
        phase_function="isotropic",
    )
    sim.run()
    counts = np.histogram(sim.results["escaped_mu"], bins=BIN_EDGES)[0].astype(float)

    interior = BIN_CENTERS > 0.1
    h_flux = _h_flux_per_bin()
    exp = h_flux / h_flux[interior].sum() * counts[interior].sum()
    chi2 = np.sum((counts[interior] - exp[interior]) ** 2 / exp[interior])
    reduced_chi2 = chi2 / (interior.sum() - 1)

    # Centered on ~1 at this budget (measured 1.23 at seed 0); 3.0 guards against
    # gross transport breakage without being flaky.
    assert reduced_chi2 < 3.0
