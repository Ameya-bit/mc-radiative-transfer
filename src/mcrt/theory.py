"""Analytic reference curves for the emergent specific intensity I(μ).

These are the theoretical limb-darkening laws a scattering atmosphere should
reproduce, used to validate the Monte Carlo beaming function.
"""

import numpy as np


def eddington_limb_darkening(mu):
    """Eddington-approximation limb darkening, I(μ) ∝ 1 + (3/2)μ (unnormalized)."""
    mu = np.asarray(mu, dtype=float)
    return 1.0 + 1.5 * mu


def chandrasekhar_h(mu, albedo=1.0, n_nodes=100, n_iter=2000, tol=1e-12):
    """Chandrasekhar H-function for isotropic scattering with single-scattering albedo.

    Solves the nonlinear integral equation

        H(μ) = 1 + (ω/2) μ H(μ) ∫₀¹ H(μ') / (μ + μ') dμ'

    by fixed-point iteration on Gauss-Legendre nodes, then evaluates H at the
    requested μ. For a conservative (ω = 1) atmosphere the emergent specific
    intensity is I(μ) ∝ H(μ); the Eddington law 1 + 1.5μ is its linear
    approximation.

    Args:
        mu: scalar or array of μ = cos(angle from surface normal), in [0, 1].
        albedo: single-scattering albedo ω (1.0 = conservative scattering).
        n_nodes: number of Gauss-Legendre quadrature nodes on [0, 1].
        n_iter: maximum fixed-point iterations.
        tol: convergence tolerance on max |ΔH| between iterations.

    Returns:
        np.ndarray of H(μ) values, same shape as the input μ.
    """
    mu_in = np.atleast_1d(np.asarray(mu, dtype=float))

    # Gauss-Legendre nodes/weights mapped from [-1, 1] to [0, 1].
    x, w = np.polynomial.legendre.leggauss(n_nodes)
    nodes = 0.5 * (x + 1.0)
    weights = 0.5 * w

    # Iterate H on the quadrature nodes.
    node_sum = nodes[:, None] + nodes[None, :]  # μ_i + μ_j, always > 0
    H = np.ones_like(nodes)
    for _ in range(n_iter):
        integral = ((weights * H)[None, :] / node_sum).sum(axis=1)
        H_new = 1.0 / (1.0 - 0.5 * albedo * nodes * integral)
        if np.max(np.abs(H_new - H)) < tol:
            H = H_new
            break
        H = H_new

    # Evaluate H at the requested μ using the converged node values.
    integral = ((weights * H)[None, :] / (mu_in[:, None] + nodes[None, :])).sum(axis=1)
    return 1.0 / (1.0 - 0.5 * albedo * mu_in * integral)
