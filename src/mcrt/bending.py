"""Exact Schwarzschild light bending — the reference the linear map approximates.

`pulse.bend` uses Beloborodov's (2002) linear map ``cos α = u + (1 − u) cos ψ``,
which is exact to ~1% for u ≲ 0.5 but degrades at grazing emission (the faint
phase of a compact star like J0740, u = 0.494). This module computes the *exact*
ray deflection instead, by numerically integrating the Schwarzschild null
geodesic, so the linear approximation can be checked against the truth and the
anchors re-run on the exact map where it matters (next-steps Track B).

The physics. A photon leaving the surface at radius R with emission angle α from
the local radial normal travels to a distant observer, arriving at the geometric
angle ψ between the spot normal and the line of sight. Gravity bends the path, so
ψ > α. With compactness ``u = R_s/R = 2GM/(Rc²)`` the impact parameter is
``b = R sin α / √(1 − u)`` and the total deflection is the standard quadrature

    ψ(α) = ∫_R^∞  (b / r²) dr / √(1 − (b/r)² (1 − R_s/r)).

Substituting ``x = R/r ∈ (0, 1]`` and ``β = b/R = sin α / √(1 − u)`` removes the
improper upper limit and gives the form integrated here,

    ψ(α) = ∫_0^1  β dx / √(1 − β² x² (1 − u x)).

For u → 0 the integrand collapses to ``sin α /√(1 − sin²α x²)`` and the integral
is ``arcsin(sin α) = α`` — flat space, no bending — which is the exactness check
in the tests. For α ≤ π/2 and u ≲ 0.5 the integrand has no interior turning point
and only an integrable √-singularity at the x = 1 endpoint when α = π/2 (grazing);
Gauss–Legendre nodes are strictly interior, so the quadrature never touches it.

ψ(α) is monotone increasing, so ``cos α`` is a monotone increasing function of
``cos ψ`` — the map `pulse.bend` supplies linearly. :class:`ExactBending`
tabulates ψ(α) once, inverts it, and exposes both ``cos α(cos ψ)`` and the exact
**lensing Jacobian** ``D(ψ) = d cos α / d cos ψ`` (the constant ``1 − u`` in the
linear map, now a function of phase). The flux weight ``F ∝ D · I(cos α) · cos α``
uses D in place of ``1 − u``; see `pulse.point_spot_flux`.

References:
    Pechenick, Ftaclas & Cohen (1983), ApJ 274, 846 — the exact deflection integral.
    Beloborodov (2002), ApJ 566, L85 — the linear approximation this validates.
    La Placa et al. (2019), RNAAS 3, 99 — modern analytic bending approximations.
"""

import numpy as np

# Table resolution: emission-angle samples and quadrature nodes. Both are cheap
# (one (n_alpha × n_quad) array, built once per compactness) and set well above
# where the inverted map or its numerical Jacobian stop changing.
DEFAULT_N_ALPHA = 4096
DEFAULT_N_QUAD = 512

# Photon-sphere limit: the surface must lie outside r = 3GM/c² = 1.5 R_s, i.e.
# u = R_s/R < 2/3. At or above it, grazing rays orbit the star, ψ(π/2) → ∞, the
# deflection integrand's radicand turns negative at interior nodes, and the map is
# not invertible — a physics boundary, not a resolution problem.
PHOTON_SPHERE_U = 2.0 / 3.0


def deflection_angle(alpha, compactness: float, n_quad: int = DEFAULT_N_QUAD) -> np.ndarray:
    """Exact Schwarzschild deflection ψ(α) by Gauss–Legendre quadrature.

    ``ψ(α) = ∫_0^1 β dx / √(1 − β² x² (1 − u x))`` with ``β = sin α / √(1 − u)``
    and ``u = compactness``. ``alpha`` is the emission angle from the radial
    normal (radians, in [0, π/2]); returns ψ (radians) of the same shape. At
    u = 0 this returns α exactly.
    """
    alpha = np.asarray(alpha, dtype=float)
    beta = np.sin(alpha) / np.sqrt(1.0 - compactness)

    # Gauss–Legendre nodes/weights mapped from [-1, 1] to (0, 1); nodes are
    # strictly interior, so the grazing endpoint singularity at x = 1 is never hit.
    t, w = np.polynomial.legendre.leggauss(n_quad)
    x = 0.5 * (t + 1.0)
    weights = 0.5 * w

    b2 = (beta**2)[..., None]                       # β², broadcast over nodes
    radicand = 1.0 - b2 * x**2 * (1.0 - compactness * x)
    integrand = beta[..., None] / np.sqrt(radicand)
    return (integrand * weights).sum(axis=-1)


class ExactBending:
    """Exact bending map for one compactness: ``cos α(cos ψ)`` and its Jacobian.

    Tabulates the deflection ψ(α) on a dense α grid, then inverts it so that a
    geometric ``cos ψ`` (from :func:`mcrt.pulse.cos_psi`) maps to the emission
    ``cos α`` and to the lensing Jacobian ``D = d cos α / d cos ψ``. Pass an
    instance as the ``bending=`` argument of :func:`mcrt.pulse.point_spot_flux`
    or :func:`mcrt.pulse.compute_profile` to swap the linear map for the exact one.

    Attributes:
        compactness: the u the map was built for.
        cos_psi_visible: smallest cos ψ that still escapes (cos α = 0 at grazing);
            the exact analogue of ``pulse.visibility_threshold``. Below it the spot
            is eclipsed and ``cos_alpha`` returns a negative value (invisible).
        psi_max: the deflection at grazing emission, ψ(π/2) — the largest angle
            still seen "around the back".
    """

    def __init__(self, compactness: float, n_alpha: int = DEFAULT_N_ALPHA,
                 n_quad: int = DEFAULT_N_QUAD):
        if not 0.0 <= compactness < PHOTON_SPHERE_U:
            raise ValueError(
                f"compactness u must be in [0, 2/3); got {compactness}. The surface "
                f"must lie outside the photon sphere (R > 3GM/c², u < 2/3) for the "
                f"deflection integral to converge — this is a physics limit, not a "
                f"quadrature-resolution one."
            )
        self.compactness = float(compactness)

        # Sample α ∈ [0, π/2] (normal → grazing) and integrate the deflection.
        alpha = np.linspace(0.0, 0.5 * np.pi, n_alpha)
        psi = deflection_angle(alpha, compactness, n_quad=n_quad)

        # ψ(α) must be strictly increasing for the inverse to exist; and ψ ≤ π so
        # that cos ψ stays monotone (true for u ≲ 0.5 — well outside the photon
        # sphere). Guard both rather than silently returning a garbled map.
        if not np.all(np.diff(psi) > 0.0):
            raise ValueError("deflection ψ(α) is not strictly monotonically increasing; "
                             "increase n_quad or check the compactness")
        if psi[-1] >= np.pi:
            raise ValueError(f"grazing deflection ψ_max = {psi[-1]:.3f} ≥ π at "
                             f"u = {compactness}; cos ψ is no longer invertible")

        self.psi_max = float(psi[-1])

        # Both cos ψ and cos α decrease as α runs 0 → π/2. Reverse to ascending
        # cos ψ so np.interp can invert the map: cos ψ ↑ ⇒ cos α ↑.
        self._cos_psi = np.cos(psi)[::-1]          # ascending: cos ψ_max … 1
        self._cos_alpha = np.cos(alpha)[::-1]      # ascending: 0 … 1
        self.cos_psi_visible = float(self._cos_psi[0])

        # Lensing Jacobian D(cos ψ) = d cos α / d cos ψ, from the tabulated map.
        # np.gradient handles the non-uniform cos ψ spacing; at cos ψ = 1 (α = 0)
        # this recovers the linear map's constant (1 − u), a self-consistency check.
        self._jacobian = np.gradient(self._cos_alpha, self._cos_psi)

        # Slope at the grazing end, used to extrapolate eclipsed phases smoothly
        # into cos α < 0 (so the visibility clip in pulse.point_spot_flux fires).
        self._grazing_slope = float(self._jacobian[0])

    def cos_alpha(self, cos_psi_vals) -> np.ndarray:
        """Emission cosine ``cos α`` for a geometric ``cos ψ`` (exact bending).

        Monotone interpolation of the inverted deflection map. For ``cos ψ`` below
        the visibility threshold the spot is eclipsed: the return is linearly
        extrapolated below zero so the ``cos α ≥ 0`` visibility test rejects it,
        matching how the linear ``pulse.bend`` goes negative there.
        """
        q = np.asarray(cos_psi_vals, dtype=float)
        out = np.interp(q, self._cos_psi, self._cos_alpha)
        # np.interp clamps to 0 below the table; replace with a negative linear
        # extrapolation so eclipsed phases read as invisible, not grazing.
        below = q < self.cos_psi_visible
        if np.any(below):
            out = np.where(below,
                           self._grazing_slope * (q - self.cos_psi_visible),
                           out)
        return out

    def jacobian(self, cos_psi_vals) -> np.ndarray:
        """Lensing Jacobian ``D(ψ) = d cos α / d cos ψ`` at a geometric ``cos ψ``.

        Replaces the linear map's constant ``1 − u`` in the point-spot flux
        weight. Held flat at the grazing slope below the visibility threshold
        (those phases are eclipsed and carry zero flux regardless).
        """
        q = np.asarray(cos_psi_vals, dtype=float)
        return np.interp(q, self._cos_psi, self._jacobian)

    def visibility_threshold(self) -> float:
        """Smallest visible ``cos ψ`` (exact analogue of ``pulse.visibility_threshold``)."""
        return self.cos_psi_visible


def bend_exact(cos_psi_vals, compactness: float, n_alpha: int = DEFAULT_N_ALPHA,
               n_quad: int = DEFAULT_N_QUAD) -> np.ndarray:
    """Exact ``cos α(cos ψ)`` — the drop-in exact analogue of ``pulse.bend``.

    Convenience wrapper that builds an :class:`ExactBending` and evaluates its
    map. Prefer constructing :class:`ExactBending` once when calling repeatedly
    (each call rebuilds the ψ(α) table). ``compactness`` is ``u = R_s/R``.
    """
    return ExactBending(compactness, n_alpha=n_alpha, n_quad=n_quad).cos_alpha(cos_psi_vals)
