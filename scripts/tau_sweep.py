"""Build the beaming-function library I(μ; τ) across a range of optical depths.

The validated engine extracts the emergent specific intensity I(μ) for a single
atmosphere thickness. Pulse-profile synthesis needs I(μ) as a *function of* the
scattering optical depth τ, because the amount of limb darkening is set by how
many times a photon scatters before escaping:

    τ ≪ 1 (thin)  → little scattering → nearly isotropic emission, slope b → 0
    τ ≫ 1 (thick) → many scatterings  → approaches Eddington / Chandrasekhar, b → ~1.5

This script sweeps τ, extracts I(μ) at each (reusing mcrt.beaming), and saves the
family of curves as a lookup table (data/beaming_library.npz) so downstream
pulse-profile code can interpolate I(μ) for a chosen τ instead of re-running the
Monte Carlo every time. Two figures visualise the trend.

Run from the repository root:  python scripts/tau_sweep.py
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from mcrt import (
    Simulation,
    chandrasekhar_h,
    eddington_limb_darkening,
    extract_intensity,
    fit_limb_darkening_slope,
)

# Fixed seed so the library is reproducible when re-checked. Threaded through the
# engine as an explicit Generator (the per-run seed offsets the convergence study
# needs are built on the same Simulation(rng=...) hook).
SEED = 20260601

TAU_VALUES = [0.1, 0.3, 1.0, 3.0, 10.0, 30.0]
NUM_PHOTONS = 200_000
N_BINS = 20

LIBRARY_PATH = 'data/beaming_library.npz'
CURVES_PATH = 'data/beaming_tau_curves.png'
SLOPE_PATH = 'data/beaming_slope_vs_tau.png'


def build_library(tau_values=TAU_VALUES, num_photons=NUM_PHOTONS, n_bins=N_BINS):
    """Run the engine at each τ and tabulate I(μ; τ).

    Returns:
        (tau_values, mu_centers, intensity_by_tau, b_of_tau)
        - mu_centers:       float64[n_bins], shared μ grid (same bins for every τ)
        - intensity_by_tau: float64[n_tau, n_bins], I(μ) normalized to 1 at μ≈1
        - b_of_tau:         float64[n_tau], fitted limb-darkening slope per τ
    """
    rng = np.random.default_rng(SEED)

    mu_centers = None
    intensity_by_tau = []
    b_of_tau = []

    for tau in tau_values:
        sim = Simulation(tau_total=tau, num_photons=num_photons, rng=rng)
        sim.run()

        centers, intensity = extract_intensity(sim.results['escaped_mu'], n_bins=n_bins)
        slope = fit_limb_darkening_slope(centers, intensity)

        mu_centers = centers  # identical across τ (range (0,1), fixed n_bins)
        intensity_by_tau.append(intensity)
        b_of_tau.append(slope)

        escaped = len(sim.results['escaped_mu'])
        print(f"  τ={tau:>5}: escaped {escaped:>7}/{num_photons}  →  slope b = {slope:.3f}")

    return (
        np.asarray(tau_values, dtype=float),
        np.asarray(mu_centers, dtype=float),
        np.asarray(intensity_by_tau, dtype=float),
        np.asarray(b_of_tau, dtype=float),
    )


def save_library(tau_values, mu_centers, intensity_by_tau, b_of_tau, path=LIBRARY_PATH):
    np.savez(
        path,
        tau_values=tau_values,
        mu_centers=mu_centers,
        intensity_by_tau=intensity_by_tau,
        b_of_tau=b_of_tau,
    )
    print(f"✓ Beaming-function library saved to {path}")


def plot_curves(tau_values, mu_centers, intensity_by_tau, path=CURVES_PATH):
    """Overlay the I(μ) curve for each τ against the theory references."""
    mu_th = np.linspace(1e-3, 1.0, 200)
    edd = eddington_limb_darkening(mu_th) / eddington_limb_darkening(1.0)
    H = chandrasekhar_h(mu_th)
    H = H / chandrasekhar_h(1.0)[0]

    colors = cm.viridis(np.linspace(0.0, 0.9, len(tau_values)))

    plt.figure(figsize=(8, 6))
    for tau, intensity, color in zip(tau_values, intensity_by_tau, colors):
        plt.plot(mu_centers, intensity, 'o-', color=color, ms=4, label=f'τ = {tau}')
    plt.plot(mu_th, edd, 'r--', lw=1.5, label='Eddington (1 + 1.5μ)')
    plt.plot(mu_th, H, 'k-', lw=1.5, label='Chandrasekhar H(μ)')
    plt.xlabel(r'$\mu = \cos(\theta)$')
    plt.ylabel(r'Specific intensity $I(\mu)$ (normalized to $\mu=1$)')
    plt.title('Beaming-function library: I(μ) vs. optical depth τ')
    plt.legend(fontsize=9)
    plt.grid(True, alpha=0.3)
    plt.savefig(path)
    print(f"✓ I(μ; τ) curve family saved to {path}")


def plot_slope_vs_tau(tau_values, b_of_tau, path=SLOPE_PATH):
    """Show the limb-darkening slope rising from ~0 (thin) toward Eddington (thick)."""
    plt.figure(figsize=(8, 5))
    plt.semilogx(tau_values, b_of_tau, 'o-', color="#2c7fb8", label='MC best-fit slope b(τ)')
    plt.axhline(1.5, color="#c0392b", ls='--', label='Eddington (b = 1.5)')
    plt.axhline(0.0, color="#7f8c8d", ls=':', label='Isotropic (b = 0)')
    plt.xlabel('Optical depth τ')
    plt.ylabel('Best-fit limb-darkening slope b')
    plt.title('Limb darkening vs. optical depth')
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.savefig(path)
    print(f"✓ slope-vs-τ plot saved to {path}")


if __name__ == "__main__":
    print(f"Building beaming-function library over τ = {TAU_VALUES} (N={NUM_PHOTONS} each)...")
    tau_values, mu_centers, intensity_by_tau, b_of_tau = build_library()
    save_library(tau_values, mu_centers, intensity_by_tau, b_of_tau)
    plot_curves(tau_values, mu_centers, intensity_by_tau)
    plot_slope_vs_tau(tau_values, b_of_tau)
