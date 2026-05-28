import numpy as np
import matplotlib.pyplot as plt
from mcrt import Simulation

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

def plot_beaming_function(tau_total=10.0, num_photons=20000):
    """Generate and plot the beaming function I(mu)."""
    print(f"Running simulation for beaming function (tau={tau_total})...")
    sim = Simulation(tau_total=tau_total, num_photons=num_photons)
    sim.run()
    
    mus = sim.results['escaped_mu']
    
    plt.figure(figsize=(8, 6))
    plt.hist(mus, bins=25, density=True, alpha=0.7, label=f'MC Result (tau={tau_total})')
    
    # Theory: For a semi-infinite atmosphere (Eddington limit), I(mu) ~ 1 + 1.5*mu
    # We normalize it for comparison
    mu_range = np.linspace(0, 1, 100)
    theory = (1 + 1.5 * mu_range) 
    # Use np.trapezoid (NumPy 2.0+) or manual integration
    theory /= np.trapezoid(theory, mu_range) 
    
    plt.plot(mu_range, theory, 'r--', label='Eddington Approximation (1 + 1.5μ)')
    
    plt.xlabel(r'$\mu = \cos(\theta)$')
    plt.ylabel('Intensity (Normalized)')
    plt.title('Beaming Function Extraction')
    plt.legend()
    plt.grid(True, alpha=0.3)
    
    output_path = 'data/beaming_function.png'
    plt.savefig(output_path)
    print(f"✓ Beaming function plot saved to {output_path}")

if __name__ == "__main__":
    validate_energy_conservation()
    validate_mean_free_path()
    plot_beaming_function()
