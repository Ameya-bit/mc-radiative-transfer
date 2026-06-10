"""Demonstrate and verify the point-spot pulse-profile machinery (v0.8.0).

Two panels, saved to data/pulse_profile_analytic.png:

  (A) Analytic check — the numerical pipeline (compute_profile) plotted on
      top of the closed-form isotropic flux F ∝ (1−u)(u + (1−u) cos ψ) for an
      always-visible geometry. The points sit on the line; the residual is shown
      to be machine-precision, which is exactly what the unit test asserts.

  (B) Gravitational smoothing — pulsed fraction vs compactness u for the same
      geometry. As u grows, light bending first lifts the spot out of eclipse
      (PF pinned at 1 while part of the orbit is hidden) and then smooths the
      pulse (PF falls). The vertical line marks the eclipse → always-visible
      transition where the closed form takes over.

This is deterministic (no Monte Carlo) — the figure is identical on every run.

Run from the repository root:  python scripts/pulse_demo.py
"""

import numpy as np
import matplotlib.pyplot as plt

from mcrt import (
    analytic_isotropic_pf,
    compute_profile,
    cos_psi,
    pulsed_fraction,
    visibility_threshold,
)

# One geometry used for both panels (i, θ_s in degrees). At low compactness the
# spot is eclipsed near φ = π (cos ψ_min = cos 105° ≈ −0.26 sits below the bending
# horizon); by u ≈ 0.3 bending lifts it into permanent view, where the Beloborodov
# closed form applies — so we run the analytic check at u = 0.3.
INCLINATION_DEG = 45.0
COLATITUDE_DEG = 60.0
COMPACTNESS_ANALYTIC = 0.3

FIGURE_PATH = "data/pulse_profile_analytic.png"


def panel_analytic(ax, ax_resid):
    """Numerical profile vs. closed-form isotropic flux, plus the residual."""
    i, theta = np.deg2rad(INCLINATION_DEG), np.deg2rad(COLATITUDE_DEG)
    u = COMPACTNESS_ANALYTIC

    prof = compute_profile(i, theta, u, n_phase=400)
    phase_cycles = prof.phase / (2.0 * np.pi)

    # Closed form, derived independently of the pipeline: isotropic ⇒ I ≡ 1.
    analytic = (1.0 - u) * (u + (1.0 - u) * cos_psi(prof.phase, i, theta))

    ax.plot(phase_cycles, analytic, "-", lw=2, color="#1f77b4",
            label=r"closed form  $(1-u)(u+(1-u)\cos\psi)$")
    ax.plot(phase_cycles[::12], prof.flux[::12], "o", ms=5, color="#d62728",
            label="compute_profile (numerical)")
    ax.set_xlabel("rotational phase  φ / 2π")
    ax.set_ylabel("flux  F(φ)  [arb.]")
    ax.set_title(f"(A) analytic check: isotropic point spot  (i={INCLINATION_DEG:.0f}°, "
                 f"θ_s={COLATITUDE_DEG:.0f}°, u={u})")
    ax.legend(loc="lower center", fontsize=8)

    resid = prof.flux - analytic
    ax_resid.plot(phase_cycles, resid, "-", color="#555555", lw=1)
    ax_resid.set_xlabel("rotational phase  φ / 2π")
    ax_resid.set_ylabel("residual")
    ax_resid.set_title(f"numerical − analytic  (max |Δ| = {np.abs(resid).max():.1e})")
    ax_resid.ticklabel_format(axis="y", style="sci", scilimits=(0, 0))


def panel_smoothing(ax):
    """Pulsed fraction vs. compactness for the same geometry."""
    i, theta = np.deg2rad(INCLINATION_DEG), np.deg2rad(COLATITUDE_DEG)
    u_grid = np.linspace(0.0, 0.45, 46)

    pf = np.array([pulsed_fraction(compute_profile(i, theta, u).flux) for u in u_grid])

    # Compactness at which the spot stops being eclipsed: cos ψ_min == −u/(1−u).
    cpsi_min = np.cos(i + theta)
    u_no_eclipse = next((u for u in u_grid if cpsi_min >= visibility_threshold(u)), None)

    ax.plot(u_grid, pf, "-", lw=2, color="#2ca02c")
    if u_no_eclipse is not None and u_no_eclipse > u_grid[0]:
        ax.axvline(u_no_eclipse, ls="--", color="#888888", lw=1)
        ax.text(u_no_eclipse + 0.006, 0.6, "spot stops\nsetting (eclipse →\nalways visible)",
                fontsize=8, color="#555555")
        # Closed-form PF overlaid where it applies (always-visible branch).
        u_av = u_grid[u_grid >= u_no_eclipse]
        pf_av = [analytic_isotropic_pf(i, theta, u) for u in u_av]
        ax.plot(u_av, pf_av, "o", ms=4, color="#1f77b4", label="closed form (always visible)")
        ax.legend(loc="upper right", fontsize=8)
    ax.set_xlabel("compactness  u = R_s/R")
    ax.set_ylabel("pulsed fraction  PF")
    ax.set_title("(B) gravitational smoothing: stronger bending → smaller PF")


def main():
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.6), layout="constrained")
    panel_analytic(axes[0], axes[1])
    panel_smoothing(axes[2])
    fig.savefig(FIGURE_PATH, dpi=130)
    print(f"saved {FIGURE_PATH}")

    # Echo the numbers the deep dive quotes, so the figure and prose stay in sync.
    i, theta = np.deg2rad(INCLINATION_DEG), np.deg2rad(COLATITUDE_DEG)
    for u in (0.0, 0.2, 0.3, 0.4):
        pf = pulsed_fraction(compute_profile(i, theta, u).flux)
        print(f"  u={u:.1f}  PF={pf:.4f}")


if __name__ == "__main__":
    main()
