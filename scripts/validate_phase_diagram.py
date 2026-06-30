"""Validate the phase-diagram machinery against the canonical antipodal limit.

The geometry phase diagram (``scripts/phase_diagram.py``) extends the standard
two-spot visibility picture into a regime the literature has not mapped: spots at
the *same* colatitude, separated by a varying azimuth Δφ, labeled by where the
beaming systematic ΔPF lands. Before trusting that extension we show the same
engine reproduces the canonical map in the limit where it is already established —
**antipodal** spots (colatitudes θ and 180−θ, half-cycle apart) — exactly the
Poutanen & Beloborodov (2006) visibility classification (their Fig. 5).

Two tiers, from the project roadmap (`docs/next-steps.md`, step 2):

1. **Single-spot tier (Beloborodov).** A point spot sets at its farthest phase
   (φ = π, where cos ψ = cos(i + θ)) exactly when bending can no longer lift it
   above the horizon — i.e. when

       cos(i + θ)  <  −u / (1 − u)      [= mcrt.visibility_threshold(u)]

   We confirm the *engine's* eclipse (from ``single_spot_eclipsed_fraction`` via
   the simulated visibility mask) flips precisely at this analytic condition,
   across a grid of (i, θ, u). This is the same closed-form check as the v0.8.0
   Rung-A benchmark, now swept over the plane.

2. **Two-spot tier (Poutanen & Beloborodov 2006, Fig. 5).** For an antipodal pair
   the secondary's geometry mirrors the primary: cos ψ₂(φ) = −cos ψ₁(φ). With the
   symmetric bending threshold c ≡ u/(1−u) the two spots' visibility windows are
   reflections, so the rotation falls into four classes set by how the primary's
   cos ψ range [m, M] = [cos(i+θ), cos(i−θ)] sits relative to [−c, c]:

       IV  both spots visible all rotation            [m, M] ⊆ [−c, c]
       III each spot eclipses for part, never both     [−c, c] ⊂ (m, M)
       II  one spot always visible, the other eclipses  one endpoint inside [−c, c]
       I   a single spot visible all rotation           [m, M] ∩ [−c, c] = ∅

   (At least one spot is always visible for an antipodal pair — there is no total
   eclipse — which is the n_min = 1 floor PB06 note.) The class boundaries are the
   curves cos(i±θ) = ±c; the boundary m = −c is the *same* Beloborodov condition as
   tier 1. We classify each (i, θ) cell two independent ways — analytically from
   (m, M, c), and numerically by counting the engine's visible-spot mask per phase —
   and confirm they agree everywhere. Reproducing PB06's map in the antipodal limit
   is what earns trust in the equal-colatitude, swept-Δφ extension.

The Zhao positioning. Zhao, Psaltis & Özel (2024) take exactly this antipodal
two-spot setup and route the wrong-beaming systematic into inferred radius via the
waveform's harmonic amplitude — without ever using the pulsed fraction, varying the
azimuthal separation, or stating a tiling criterion. The third panel here is the
direct side-by-side: on the antipodal geometry, ΔPF as a function of Δφ, with the
tiling boundary marked — the three layers Zhao omits, made concrete.

Run from the repository root:  python3 scripts/validate_phase_diagram.py
"""

from typing import NamedTuple

import numpy as np
import matplotlib.pyplot as plt

from mcrt import (
    beaming_lookup,
    compute_profile,
    pulsed_fraction,
    visibility_threshold,
)
from anchor_lib import (
    N_PHASE,
    Spot,
    load_library,
    multi_spot_flux,
    shape_tau_index,
    single_spot_eclipsed_fraction,
)

FIGURE_PATH = "data/validate_phase_diagram.png"

# Class codes (integers so they tile a heatmap cleanly). Roman numerals follow the
# physical reading above; PB06 Fig. 5 partitions the (i, θ) plane into the same four
# regions by the visibility curves cos(i±θ) = ±u/(1−u).
CLASS_I, CLASS_II, CLASS_III, CLASS_IV = 1, 2, 3, 4
CLASS_NAMES = {
    CLASS_I: "I  · one spot, always",
    CLASS_II: "II · one always, one eclipses",
    CLASS_III: "III · both eclipse, never both",
    CLASS_IV: "IV · both, always",
}

# (i, θ) sweep for the class map and the Beloborodov grid check. Degrees on [0,180]
# for θ (full colatitude) and [0,90] for i (observer colatitude; the map is
# symmetric about i = 90°, so the upper half is redundant).
I_AXIS_DEG = np.linspace(1.0, 89.0, 89)
THETA_AXIS_DEG = np.linspace(1.0, 179.0, 90)

# Representative compact spacetime for the class map / Zhao panel (J0740-like).
VALID_U = 0.494


# --- Tier 1: Beloborodov single-spot eclipse condition --------------------------

def single_spot_sets_analytic(inclination, colatitude, compactness) -> bool:
    """Closed form: the spot dips below the bending horizon at its farthest phase.

    The minimum cos ψ over a rotation is cos(i + θ) (at φ = π). The spot is
    eclipsed for some phase exactly when that minimum is below the visibility
    threshold −u/(1−u).
    """
    return bool(np.cos(inclination + colatitude) < visibility_threshold(compactness))


def beloborodov_grid_mismatches(i_axis=I_AXIS_DEG, theta_axis=THETA_AXIS_DEG,
                                compactnesses=(0.1, 0.3, 0.494)) -> int:
    """Count (i, θ, u) cells where the engine and the analytic condition disagree.

    For each cell, the engine eclipses iff ``single_spot_eclipsed_fraction`` > 0;
    the analytic prediction is :func:`single_spot_sets_analytic`. A correct engine
    yields zero mismatches across the whole grid.
    """
    mismatches = 0
    for u in compactnesses:
        for i_deg in i_axis:
            for th_deg in theta_axis:
                incl, theta = np.deg2rad(i_deg), np.deg2rad(th_deg)
                engine_sets = single_spot_eclipsed_fraction(incl, theta, u) > 0.0
                if engine_sets != single_spot_sets_analytic(incl, theta, u):
                    mismatches += 1
    return mismatches


# --- Tier 2: Poutanen & Beloborodov (2006) antipodal visibility classes ---------

def antipodal_class_analytic(inclination, colatitude, compactness) -> int:
    """PB06 visibility class of an antipodal pair from the closed-form geometry.

    The primary's cos ψ spans [m, M] = [cos(i+θ), cos(i−θ)]; the antipode's spans
    [−M, −m]. With threshold c = u/(1−u) both spots are visible where the primary's
    cos ψ ∈ [−c, c]. The four classes follow from how [m, M] overlaps [−c, c].
    """
    c = compactness / (1.0 - compactness)
    m = float(np.cos(inclination + colatitude))      # min cos ψ of primary
    big_m = float(np.cos(inclination - colatitude))  # max cos ψ of primary
    # Primary is visible where cos ψ ≥ −c; the antipode where cos ψ ≤ c
    # (its cos ψ = −cos ψ_primary). Each spot is "always" / "never" / "partial".
    p_always, p_never = m >= -c, big_m < -c
    s_always, s_never = big_m <= c, m > c
    if p_always and s_always:
        return CLASS_IV                          # both up all rotation
    if (not p_always and not p_never) and (not s_always and not s_never):
        return CLASS_III                         # both eclipse for part, never both
    if (p_always and s_never) or (s_always and p_never):
        return CLASS_I                           # only one spot ever visible
    return CLASS_II                              # one always up, the other eclipses


def antipodal_class_numerical(inclination, colatitude, compactness,
                              n_phase=N_PHASE) -> int:
    """Same class, read off the engine: count visible spots per phase.

    Builds the primary (colatitude θ) and antipode (colatitude π−θ, rolled half a
    cycle) visibility masks from ``compute_profile`` and classifies by the per-phase
    visible-spot count (its min and max over the rotation, and whether a both-visible
    interval coexists with single-visible intervals).
    """
    prim = compute_profile(inclination, colatitude, compactness, n_phase=n_phase)
    anti = compute_profile(inclination, np.pi - colatitude, compactness, n_phase=n_phase)
    anti_vis = np.roll(anti.visible, n_phase // 2)   # antipode is half a cycle away
    n_vis = prim.visible.astype(int) + anti_vis.astype(int)
    n_min, n_max = int(n_vis.min()), int(n_vis.max())
    if n_min == 2:
        return CLASS_IV               # always both
    if n_max == 1:
        return CLASS_I                # never both → only one spot ever seen
    # n_max == 2 and n_min == 1: both-visible interval plus single-visible interval.
    # Class III iff each spot is individually eclipsed somewhere; else Class II.
    both_eclipse = (not prim.visible.all()) and (not anti_vis.all())
    return CLASS_III if both_eclipse else CLASS_II


def class_maps(compactness=VALID_U, i_axis=I_AXIS_DEG, theta_axis=THETA_AXIS_DEG):
    """Analytic and numerical class grids over (θ, i); plus their mismatch count."""
    n_i, n_t = len(i_axis), len(theta_axis)
    analytic = np.zeros((n_i, n_t), dtype=int)
    numerical = np.zeros((n_i, n_t), dtype=int)
    for ii, i_deg in enumerate(i_axis):
        for jt, th_deg in enumerate(theta_axis):
            incl, theta = np.deg2rad(i_deg), np.deg2rad(th_deg)
            analytic[ii, jt] = antipodal_class_analytic(incl, theta, compactness)
            numerical[ii, jt] = antipodal_class_numerical(incl, theta, compactness)
    mismatches = int(np.count_nonzero(analytic != numerical))
    return analytic, numerical, mismatches


# --- Zhao positioning: the three layers they omit, on their own geometry --------

def zhao_extension(beaming_shape, colatitude_deg=60.0, inclination_deg=60.0,
                   compactness=VALID_U, dphi_axis=None):
    """ΔPF vs azimuthal separation for an antipodal pair — Zhao's setup, extended.

    Zhao et al. (2024) use this two-antipodal-spot geometry but report a harmonic
    amplitude routed into radius, never the pulsed fraction and never a separation
    sweep. Here we vary Δφ across [0, 0.5] for the antipodal colatitude pair (θ,
    180−θ) and report PF_iso, PF_real, and ΔPF, marking where the pair tiles.
    """
    if dphi_axis is None:
        dphi_axis = np.linspace(0.0, 0.5, 51)
    incl = np.deg2rad(inclination_deg)
    theta = np.deg2rad(colatitude_deg)
    anti = np.pi - theta
    out = []
    for dphi in dphi_axis:
        spots = (Spot(theta, 0.0, 1.0), Spot(anti, float(dphi), 1.0))
        iso = multi_spot_flux(incl, compactness, spots, None)
        real = multi_spot_flux(incl, compactness, spots, beaming_shape)
        out.append((float(dphi), pulsed_fraction(iso), pulsed_fraction(real),
                    pulsed_fraction(real) - pulsed_fraction(iso)))
    return np.array(out)  # columns: dphi, pf_iso, pf_real, delta_pf


class Validation(NamedTuple):
    belo_mismatches: int
    analytic: np.ndarray
    numerical: np.ndarray
    class_mismatches: int
    zhao: np.ndarray


def run() -> Validation:
    tau_values, mu_centers, intensity_by_tau = load_library()
    beaming_shape = beaming_lookup(mu_centers, intensity_by_tau[shape_tau_index(tau_values)])
    belo = beloborodov_grid_mismatches()
    analytic, numerical, class_mm = class_maps()
    zhao = zhao_extension(beaming_shape)
    return Validation(belo, analytic, numerical, class_mm, zhao)


def plot_result(v: Validation, path=FIGURE_PATH, compactness=VALID_U):
    """Three panels: PB06 class map (analytic), the engine's class map, ΔPF(Δφ)."""
    fig, axes = plt.subplots(1, 3, figsize=(17, 5.0))
    fig.subplots_adjust(wspace=0.45)
    cmap = plt.cm.viridis
    extent = [THETA_AXIS_DEG[0], THETA_AXIS_DEG[-1], I_AXIS_DEG[0], I_AXIS_DEG[-1]]

    im = None
    for ax, grid, title in (
        (axes[0], v.analytic, "PB06 Fig. 5 classes — analytic"),
        (axes[1], v.numerical, "PB06 Fig. 5 classes — from the engine"),
    ):
        im = ax.imshow(grid, origin="lower", extent=extent, aspect="auto",
                       cmap=cmap, vmin=CLASS_I - 0.5, vmax=CLASS_IV + 0.5)
        ax.set_xlabel(r"Spot colatitude  $\theta$  (deg)")
        ax.set_ylabel(r"Observer inclination  $i$  (deg)")
        ax.set_title(f"{title}\n(antipodal pair, u = {compactness:g})", fontsize=10)
    cbar = fig.colorbar(im, ax=axes[:2], fraction=0.040, pad=0.015,
                        ticks=[CLASS_I, CLASS_II, CLASS_III, CLASS_IV])
    cbar.ax.set_yticklabels(["I", "II", "III", "IV"], fontsize=9)
    cbar.set_label("visibility class", fontsize=9)

    ax = axes[2]
    dphi, pf_iso, pf_real, dpf = v.zhao.T
    ax.plot(dphi, pf_iso, color="#7f8c8d", lw=2.0, label=r"PF$_{\rm iso}$")
    ax.plot(dphi, pf_real, color="#c0392b", lw=2.0, label=r"PF$_{\rm real}$")
    ax.plot(dphi, dpf, color="#2980b9", lw=2.4, label=r"$\Delta$PF")
    ax.axhline(0.0, color="k", lw=0.6, alpha=0.4)
    ax.set_xlabel(r"Azimuthal separation  $\Delta\phi$  (cycles)")
    ax.set_ylabel("pulsed fraction")
    ax.set_title("Zhao extension: the three layers they omit\n"
                 "(antipodal pair, PF + swept Δφ + tiling)", fontsize=10)
    ax.legend(fontsize=9, loc="center right")
    ax.set_xlim(0.0, 0.5)

    fig.suptitle("Validation: the engine reproduces the canonical antipodal visibility "
                 "map (Poutanen & Beloborodov 2006), and the Zhao geometry extended to "
                 "ΔPF vs. separation", fontsize=11)
    fig.text(0.5, 0.005,
             "Classes:  I = one spot visible all rotation   ·   II = one always visible, "
             "the other eclipses   ·   III = each eclipses for part, never both at once   "
             "·   IV = both visible all rotation",
             ha="center", fontsize=8.5)
    fig.savefig(path, dpi=130, bbox_inches="tight")
    print(f"✓ figure saved to {path}")


def print_summary(v: Validation):
    print("\n" + "=" * 74)
    print("Phase-diagram validation — Beloborodov single-spot + PB06 antipodal map")
    print("=" * 74)
    print(f"\nTier 1 (Beloborodov eclipse condition): {v.belo_mismatches} mismatches "
          f"across the (i, θ, u) grid  →  {'PASS' if v.belo_mismatches == 0 else 'FAIL'}")
    print(f"Tier 2 (PB06 Fig. 5 antipodal classes): {v.class_mismatches} cells where the "
          f"engine class ≠ analytic class  →  {'PASS' if v.class_mismatches == 0 else 'FAIL'}")
    present = sorted(set(v.analytic.flatten().tolist()))
    print("  classes present in the map: " + ", ".join(CLASS_NAMES[k] for k in present))
    dphi, pf_iso, pf_real, dpf = v.zhao.T
    print("\nZhao extension (antipodal pair, the three layers they omit):")
    print(f"  Δφ = 0.00 (overlap): PF_iso = {pf_iso[0]:.2f}, ΔPF = {dpf[0]:+.2f}")
    print(f"  Δφ = 0.50 (tiling):  PF_iso = {pf_iso[-1]:.2f}, ΔPF = {dpf[-1]:+.2f}")
    print("  → varying the azimuthal separation moves the systematic between hidden "
          "and live;\n    pulsed fraction is the axis Zhao's harmonic-amplitude analysis "
          "never uses.")


if __name__ == "__main__":
    v = run()
    plot_result(v)
    print_summary(v)
