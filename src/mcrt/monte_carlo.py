import numpy as np
from .utils import (
    sample_step_size,
    sample_thomson_angle,
    sample_isotropic_angle,
    rotate_vector,
    get_random_direction,
)

# Opt-in scattering phase functions. "thomson" (the default) is the dipole
# P(mu) = (3/4)(1 + mu^2) the atmosphere physically uses; "isotropic" is the
# P(mu) = 1/2 phase function for which Chandrasekhar's H(mu) is the exact
# emergent intensity, used to certify the transport machinery (Track D).
SCATTER_SAMPLERS = {
    "thomson": sample_thomson_angle,
    "isotropic": sample_isotropic_angle,
}

class Photon:
    """
    Tracks the state of a single photon packet in the atmosphere.
    Coordinates are in optical depth units.
    """
    def __init__(self, tau_injection, direction=None):
        # Position: (x, y, tau)
        # We assume x, y start at 0. tau is the vertical optical depth.
        self.pos = np.array([0.0, 0.0, tau_injection])
        
        # Direction: unit vector (dx, dy, d_tau)
        if direction is None:
            # Default: inject vertically upward (-tau direction)
            self.dir = np.array([0.0, 0.0, -1.0])
        else:
            self.dir = direction / np.linalg.norm(direction)
            
        self.alive = True
        self.escaped = False
        self.num_scatterings = 0
        self.total_path_length = 0.0

    def propagate(self, step_size):
        """Move the photon along its current direction."""
        self.pos += self.dir * step_size
        self.total_path_length += step_size

    def scatter(self, cos_theta, rng=None):
        """Update direction based on scattering angle."""
        self.dir = rotate_vector(self.dir, cos_theta, rng=rng)
        self.num_scatterings += 1

class Simulation:
    """
    Manages the Monte Carlo simulation of many photons.
    """
    def __init__(self, tau_total, num_photons=1000, rng=None, phase_function="thomson"):
        self.tau_total = tau_total
        self.num_photons = num_photons
        # Scattering phase function. Default "thomson" leaves every existing run
        # bit-for-bit unchanged (same sampler, same RNG draws); "isotropic" is
        # the opt-in Track D mode validated against Chandrasekhar H(μ).
        if phase_function not in SCATTER_SAMPLERS:
            raise ValueError(
                f"phase_function must be one of {sorted(SCATTER_SAMPLERS)}, "
                f"got {phase_function!r}"
            )
        self.phase_function = phase_function
        self._scatter_sampler = SCATTER_SAMPLERS[phase_function]
        # An explicit Generator makes a run reproducible; convergence studies
        # pass per-run seed offsets. Default to a fresh Generator so each
        # Simulation is independent without touching global np.random state.
        self.rng = rng if rng is not None else np.random.default_rng()
        self.photons = []
        self.results = {
            'escaped_mu': [],
            'num_scatterings': [],
            'total_path_lengths': [],
            'absorbed_count': 0
        }

    def run(self):
        """Run the simulation for all photons."""
        for _ in range(self.num_photons):
            # Inject photon at the bottom boundary (tau = tau_total), moving upward.
            # The thermal source is isotropic in *specific intensity*, so the photon
            # number crossing the boundary per unit mu goes as N(mu) ∝ mu (the cos θ
            # projection). Sampling costheta = sqrt(U) reproduces that ∝ mu law;
            # uniform(0,1) would instead over-produce grazing photons (isotropic in
            # solid angle, not in intensity) — see docs/deep-dives/v0.6.1.
            phi = self.rng.uniform(0, 2 * np.pi)
            costheta = np.sqrt(self.rng.uniform(0, 1))  # 0 to 1 is upward (-z in tau coords)
            sintheta = np.sqrt(1 - costheta**2)
            
            # Note: In our coordinate system, tau decreases upwards.
            # So d_tau = -costheta
            direction = np.array([
                sintheta * np.cos(phi),
                sintheta * np.sin(phi),
                -costheta
            ])
            
            p = Photon(self.tau_total, direction)
            
            while p.alive:
                # 1. Sample step size
                d_tau = sample_step_size(rng=self.rng)
                
                # 2. Propagate
                p.propagate(d_tau)
                
                # 3. Check boundaries
                if p.pos[2] <= 0:
                    # Escaped from the top
                    p.alive = False
                    p.escaped = True
                    # mu = cos(theta) = dz / norm (but here dz is d_tau)
                    # The escape angle theta is measured relative to the normal (vertical)
                    # Our d_tau is negative for upward. So mu = -p.dir[2]
                    self.results['escaped_mu'].append(-p.dir[2])
                    self.results['num_scatterings'].append(p.num_scatterings)
                    self.results['total_path_lengths'].append(p.total_path_length)
                
                elif p.pos[2] >= self.tau_total:
                    # Re-entered the bottom (absorbed)
                    p.alive = False
                    p.escaped = False
                    self.results['absorbed_count'] += 1
                
                else:
                    # 4. Scatter (phase function selected at construction;
                    # "thomson" default preserves the original RNG stream)
                    mu_scatter = self._scatter_sampler(rng=self.rng)
                    p.scatter(mu_scatter, rng=self.rng)
            
            self.photons.append(p)

if __name__ == "__main__":
    # Quick test run
    sim = Simulation(tau_total=1.0, num_photons=1000)
    sim.run()
    print(f"Simulation finished.")
    print(f"Escaped: {len(sim.results['escaped_mu'])}")
    print(f"Absorbed: {sim.results['absorbed_count']}")
    print(f"Average scatterings: {np.mean(sim.results['num_scatterings']) if sim.results['num_scatterings'] else 0}")
