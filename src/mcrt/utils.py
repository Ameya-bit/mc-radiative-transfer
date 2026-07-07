"""Utility functions and constants for Monte Carlo Radiative Transfer.

(Using optical depth coordinates, so many constants are normalized to 1.)
"""

import numpy as np

# All samplers accept an optional `rng` (a numpy.random.Generator). When it is
# None they fall back to the legacy global `np.random` module, so direct callers
# (and the existing unit tests) keep their original behaviour; the Simulation
# threads an explicit Generator through for reproducible, independent runs.

def sample_step_size(rng=None):
    """
    Sample a step size delta_tau from the distribution P(tau) = exp(-tau).
    This corresponds to the distance a photon travels before scattering.

    Args:
        rng: Optional numpy.random.Generator. Defaults to global np.random.

    Returns:
        float: The sampled optical depth step size.
    """
    if rng is None:
        rng = np.random
    return -np.log(rng.random())

def sample_thomson_angle(rng=None):
    """
    Sample the scattering angle theta from the Thomson phase function:
    P(mu) = (3/4) * (1 + mu^2), where mu = cos(theta).

    We use rejection sampling for simplicity and clarity.

    Args:
        rng: Optional numpy.random.Generator. Defaults to global np.random.

    Returns:
        float: The sampled cosine of the scattering angle (mu).
    """
    if rng is None:
        rng = np.random
    while True:
        mu = rng.uniform(-1, 1)
        p_mu = 0.75 * (1 + mu**2)
        if rng.random() < p_mu / 1.5:  # max(p_mu) = 1.5
            return mu

def sample_isotropic_angle(rng=None):
    """
    Sample the scattering-angle cosine for *isotropic* scattering:
    P(mu) = 1/2 (constant), so mu = cos(theta) is uniform on [-1, 1].

    Isotropic scattering redistributes photons uniformly over solid angle. It is
    the phase function for which Chandrasekhar's H-function gives the *exact*
    emergent specific intensity of a conservative semi-infinite atmosphere.
    Swapping the engine onto it (from the default Thomson dipole
    P(mu) = (3/4)(1 + mu^2), see :func:`sample_thomson_angle`) lets the transport
    and geometry machinery be validated against a closed-form solution with the
    phase function controlled — attention item 2 / Track D.

    Args:
        rng: Optional numpy.random.Generator. Defaults to global np.random.

    Returns:
        float: The sampled cosine of the scattering angle (mu).
    """
    if rng is None:
        rng = np.random
    return rng.uniform(-1.0, 1.0)


def get_random_direction(rng=None):
    """
    Generate a random unit vector in 3D (isotropic distribution).

    Args:
        rng: Optional numpy.random.Generator. Defaults to global np.random.

    Returns:
        np.array: A 3D unit vector [dx, dy, dz].
    """
    if rng is None:
        rng = np.random
    phi = rng.uniform(0, 2 * np.pi)
    costheta = rng.uniform(-1, 1)
    sintheta = np.sqrt(1 - costheta**2)

    dx = sintheta * np.cos(phi)
    dy = sintheta * np.sin(phi)
    dz = costheta

    return np.array([dx, dy, dz])

def rotate_vector(vector, costheta_scatter, rng=None):
    """
    Rotate a direction vector by a scattering angle theta.
    The azimuthal angle of scattering is assumed to be uniform [0, 2pi].

    Args:
        vector (np.array): Current direction unit vector.
        costheta_scatter (float): Cosine of the scattering angle.
        rng: Optional numpy.random.Generator. Defaults to global np.random.

    Returns:
        np.array: The new direction unit vector.
    """
    if rng is None:
        rng = np.random
    # Sample azimuthal angle
    phi = rng.uniform(0, 2 * np.pi)
    sintheta_scatter = np.sqrt(1 - costheta_scatter**2)
    
    # If the vector is nearly vertical, handle it specially to avoid division by zero
    if abs(vector[2]) > 0.999:
        return np.array([
            sintheta_scatter * np.cos(phi),
            sintheta_scatter * np.sin(phi),
            costheta_scatter if vector[2] > 0 else -costheta_scatter
        ])
    
    # General Case: Rotate the vector
    # Using the standard transformation to a new coordinate system centered on the current direction
    v = vector
    sin_phi_scatter = np.sin(phi)
    cos_phi_scatter = np.cos(phi)
    
    sqrt_1_minus_vz2 = np.sqrt(1 - v[2]**2)
    
    new_dx = (sintheta_scatter * (v[0] * v[2] * cos_phi_scatter - v[1] * sin_phi_scatter) / sqrt_1_minus_vz2) + v[0] * costheta_scatter
    new_dy = (sintheta_scatter * (v[1] * v[2] * cos_phi_scatter + v[0] * sin_phi_scatter) / sqrt_1_minus_vz2) + v[1] * costheta_scatter
    new_dz = -sintheta_scatter * sqrt_1_minus_vz2 * cos_phi_scatter + v[2] * costheta_scatter
    
    new_dir = np.array([new_dx, new_dy, new_dz])
    return new_dir / np.linalg.norm(new_dir)
