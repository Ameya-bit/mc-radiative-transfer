"""Photon-count parameter study for the beaming function.

Re-runs the simulation at increasing photon counts and tracks the best-fit
limb-darkening slope b (Eddington predicts 1.5). The shape should converge as N
grows while the Monte Carlo noise shrinks — answering the reviewer's question
"does the fit change with more photons?" (the trend is convergence, not drift).
"""

import numpy as np
import matplotlib.pyplot as plt
from mcrt import Simulation, extract_intensity, fit_limb_darkening_slope

PHOTON_COUNTS = [2_000, 10_000, 50_000, 200_000]
TAU_TOTAL = 10.0
EDDINGTON_SLOPE = 1.5


def run_study(photon_counts=PHOTON_COUNTS, tau_total=TAU_TOTAL):
    slopes = []
    for n in photon_counts:
        sim = Simulation(tau_total=tau_total, num_photons=n)
        sim.run()
        mu_centers, intensity = extract_intensity(sim.results['escaped_mu'])
        slope = fit_limb_darkening_slope(mu_centers, intensity)
        slopes.append(slope)
        print(f"  N={n:>8}: best-fit slope b = {slope:.3f}")
    return np.array(photon_counts), np.array(slopes)


def plot_convergence(photon_counts, slopes, output_path='data/beaming_convergence.png'):
    plt.figure(figsize=(8, 5))
    plt.semilogx(photon_counts, slopes, 'o-', color="#2c7fb8", label='MC best-fit slope b')
    plt.axhline(EDDINGTON_SLOPE, color="#c0392b", ls='--', label='Eddington (b = 1.5)')
    plt.xlabel('Number of photons')
    plt.ylabel('Best-fit limb-darkening slope b')
    plt.title(f'Beaming-function convergence (τ = {TAU_TOTAL})')
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.savefig(output_path)
    print(f"✓ Convergence plot saved to {output_path}")


if __name__ == "__main__":
    print("Running photon-count parameter study...")
    counts, slopes = run_study()
    plot_convergence(counts, slopes)
