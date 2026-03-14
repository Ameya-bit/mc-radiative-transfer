import numpy as np

"""
Utility functions and constants for Monte Carlo Radiative Transfer.
"""

# Physical Constants
# (Using optical depth coordinates, so many constants are normalized to 1)

def sample_step_size():
    """
    Sample a step size delta_tau from the distribution P(tau) = exp(-tau).
    This corresponds to the distance a photon travels before scattering.
    
    Returns:
        float: The sampled optical depth step size.
    """
    return -np.log(np.random.random())

def sample_thomson_angle():
    """
    Sample the scattering angle theta from the Thomson phase function:
    P(mu) = (3/4) * (1 + mu^2), where mu = cos(theta).
    
    We use rejection sampling for simplicity and clarity.
    
    Returns:
        float: The sampled cosine of the scattering angle (mu).
    """
    while True:
        mu = np.random.uniform(-1, 1)
        p_mu = 0.75 * (1 + mu**2)
        if np.random.random() < p_mu / 1.5:  # max(p_mu) = 1.5
            return mu

def get_random_direction():
    """
    Generate a random unit vector in 3D (isotropic distribution).
    
    Returns:
        np.array: A 3D unit vector [dx, dy, dz].
    """
    phi = np.random.uniform(0, 2 * np.pi)
    costheta = np.random.uniform(-1, 1)
    sintheta = np.sqrt(1 - costheta**2)
    
    dx = sintheta * np.cos(phi)
    dy = sintheta * np.sin(phi)
    dz = costheta
    
    return np.array([dx, dy, dz])

def rotate_vector(vector, costheta_scatter):
    """
    Rotate a direction vector by a scattering angle theta.
    The azimuthal angle of scattering is assumed to be uniform [0, 2pi].
    
    Args:
        vector (np.array): Current direction unit vector.
        costheta_scatter (float): Cosine of the scattering angle.
        
    Returns:
        np.array: The new direction unit vector.
    """
    # Sample azimuthal angle
    phi = np.random.uniform(0, 2 * np.pi)
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
