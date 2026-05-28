import numpy as np
import matplotlib.pyplot as plt
from mcrt import (
    Simulation,
    chandrasekhar_h,
    eddington_limb_darkening,
    extract_intensity,
    fit_limb_darkening_slope,
)

def validate_energy_conservation(tau_total=1.0, num_photons=5000):
    """Verify that every photon is either escaped or absorbed."""
    sim = Simulation(tau_total=tau_total, num_photons=num_photons)
    sim.run()
    
    total_tracked = len(sim.results['escaped_mu']) + sim.results['absorbed_count']
    assert total_tracked == num_photons, f"Energy conservation failed: {total_tracked} tracked vs {num_photons} injected"
    print(f"✓ Energy Conservation: {total_tracked}/{num_photons} photons accounted for.")

def validate_mean_free_path(tau_total=10.0, num_photons=1000):
    """Verify that the average path length between scatterings is ~1.0 in tau units."""
    sim = Simulation(tau_total=tau_total, num_photons=num_photons)
    sim.run()
    
    all_scatterings = np.sum(sim.results['num_scatterings'])
    all_path_lengths = np.sum(sim.results['total_path_lengths'])
    
    if all_scatterings > 0:
        mfp = all_path_lengths / all_scatterings
        print(f"✓ Mean Free Path: {mfp:.4f} (Expected ~1.0)")
        # Note: In a finite slab, there's a slight bias because the last step to escape
        # doesn't end in a scattering. But for large tau, it should be close to 1.
        assert 0.9 <= mfp <= 1.1, f"MFP {mfp} is too far from 1.0"
    else:
        print("! Mean Free Path: No scatterings occurred (check tau_total).")

def plot_beaming_function(tau_total=10.0, num_photons=200000, n_bins=20):
    """Extract the beaming function I(μ) and compare to Eddington and Chandrasekhar."""
    print(f"Running simulation for beaming function (tau={tau_total}, N={num_photons})...")
    sim = Simulation(tau_total=tau_total, num_photons=num_photons)
    sim.run()

    mu_centers, intensity = extract_intensity(sim.results['escaped_mu'], n_bins=n_bins)
    slope = fit_limb_darkening_slope(mu_centers, intensity)
    print(f"  Best-fit limb-darkening slope b = {slope:.3f} (Eddington predicts 1.5)")

    # Theory curves, normalized to 1 at μ=1 to match the MC normalization.
    mu_th = np.linspace(1e-3, 1.0, 200)
    edd = eddington_limb_darkening(mu_th) / eddington_limb_darkening(1.0)
    H = chandrasekhar_h(mu_th)
    H = H / chandrasekhar_h(1.0)[0]

    plt.figure(figsize=(8, 6))
    plt.plot(mu_centers, intensity, 'o', color="#2c7fb8", label=f'MC intensity I(μ) (τ={tau_total})')
    plt.plot(mu_th, edd, 'r--', label='Eddington (1 + 1.5μ)')
    plt.plot(mu_th, H, 'g-', label='Chandrasekhar H(μ)')
    plt.xlabel(r'$\mu = \cos(\theta)$')
    plt.ylabel(r'Specific intensity $I(\mu)$ (normalized to $\mu=1$)')
    plt.title('Beaming Function: specific intensity vs. theory')
    plt.legend()
    plt.grid(True, alpha=0.3)

    output_path = 'data/beaming_function.png'
    plt.savefig(output_path)
    print(f"✓ Beaming function plot saved to {output_path}")

if __name__ == "__main__":
    validate_energy_conservation()
    validate_mean_free_path()
    plot_beaming_function()
