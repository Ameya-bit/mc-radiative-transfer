import numpy as np
import matplotlib.pyplot as plt
from monte_carlo import Simulation

def plot_random_walks(tau_total=5.0, num_photons=10):
    """Visualize 3D paths of a few photons."""
    sim = Simulation(tau_total=tau_total, num_photons=num_photons)
    sim.run()
    
    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(111, projection='3d')
    
    for i, p in enumerate(sim.photons):
        # We need to modify the Photon class to store the history of steps
        # or just run a mini-simulation here that stores them.
        # Since the current Photon class doesn't store history, 
        # let's write a small wrapper to track paths for this visualization.
        pass

    # Actually, let's just make a separate path tracker for this demo
    for _ in range(num_photons):
        path = [[0, 0, tau_total]]
        pos = np.array([0.0, 0.0, tau_total])
        
        # Initial upward direction
        phi = np.random.uniform(0, 2 * np.pi)
        costheta = np.random.uniform(0, 1)
        sintheta = np.sqrt(1 - costheta**2)
        direction = np.array([sintheta * np.cos(phi), sintheta * np.sin(phi), -costheta])
        
        alive = True
        while alive:
            d_tau = -np.log(np.random.random())
            pos = pos + direction * d_tau
            path.append(pos.copy())
            
            if pos[2] <= 0 or pos[2] >= tau_total:
                alive = False
            else:
                # Scatter (Thomson)
                while True:
                    mu = np.random.uniform(-1, 1)
                    if np.random.random() < (0.75 * (1 + mu**2)) / 1.5:
                        break
                
                # Rotate (using simplified rotation for visualization script)
                phi_s = np.random.uniform(0, 2 * np.pi)
                sin_s = np.sqrt(1 - mu**2)
                
                if abs(direction[2]) > 0.999:
                    direction = np.array([sin_s * np.cos(phi_s), sin_s * np.sin(phi_s), mu if direction[2] > 0 else -mu])
                else:
                    v = direction
                    sqrt_1_vz2 = np.sqrt(1 - v[2]**2)
                    dx = (sin_s * (v[0] * v[2] * np.cos(phi_s) - v[1] * np.sin(phi_s)) / sqrt_1_vz2) + v[0] * mu
                    dy = (sin_s * (v[1] * v[2] * np.cos(phi_s) + v[0] * np.sin(phi_s)) / sqrt_1_vz2) + v[1] * mu
                    dz = -sin_s * sqrt_1_vz2 * np.cos(phi_s) + v[2] * mu
                    direction = np.array([dx, dy, dz])
                    direction /= np.linalg.norm(direction)

        path = np.array(path)
        ax.plot(path[:,0], path[:,1], path[:,2], alpha=0.6)

    # Plot boundaries
    ax.set_zlim(0, tau_total)
    ax.set_xlabel('X (Optical Depth)')
    ax.set_ylabel('Y (Optical Depth)')
    ax.set_zlabel('Tau (Vertical Depth)')
    ax.set_title(f'Photon Random Walks (tau_total={tau_total})')
    ax.invert_zaxis() # Tau=0 is top
    
    output_path = 'figures/photon_paths.png'
    plt.savefig(output_path)
    print(f"✓ Path visualization saved to {output_path}")

if __name__ == "__main__":
    plot_random_walks()
