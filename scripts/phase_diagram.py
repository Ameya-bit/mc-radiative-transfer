"""Geometry phase diagram: where the beaming systematic lands, as a general rule.

The two real-star anchors are single points in a much larger geometry space:
v0.9.1 (J0030) found the isotropic→realistic systematic *saturated* out of the
pulsed fraction (ΔPF ≈ 0, the effect hiding in waveform shape), while v0.9.2
(J0740) found it *live* in the pulsed fraction (ΔPF ≈ +0.16…+0.23). The headline
claim of the project is that the discriminator between those two outcomes is not
whether a single spot eclipses but whether the two spots **tile** the rotation.
This script turns that claim from two examples into a map.

The canonical model. Two equal-weight point spots at the **same** colatitude θ,
separated by azimuth **Δφ**, at a fixed spacetime (inclination i, compactness u).
Sweeping (Δφ, θ) and coloring by ΔPF = PF_real − PF_iso (at the profile-shape
optical depth τ = SHAPE_TAU) shows the systematic's observability across geometry
in one picture. Everything is assembled from the verified core via
:mod:`anchor_lib` (each spot is one :func:`mcrt.compute_profile`; a spot at
longitude φ₀ is that profile rolled in phase; light is additive) — **no change to
``mcrt.pulse``**.

The analytic boundary (the tiling criterion, made into a curve). Two identical
anti-phased spots each spend a fraction ``f_ecl(θ)`` of the cycle eclipsed behind
the star, their dark windows centered half a turn apart. The combined flux touches
zero — pinning PF = 1 for *both* beamings, so ΔPF collapses to ≈ 0 — exactly when
those two dark windows overlap, i.e. when

    Δφ  <  f_ecl(θ)          (saturated: systematic hides in waveform shape)
    Δφ  >  f_ecl(θ)          (tiling:   systematic is live in the pulsed fraction)

So the saturation boundary is simply ``Δφ_crit(θ) = f_ecl(θ)``, the single-spot
eclipse fraction — already computed by
:func:`anchor_lib.single_spot_eclipsed_fraction`. The numerical flux-floor → 0
contour is overlaid on top of it as an independent check that the two agree.

The real stars, placed on the map. J0030 (both fits) and J0740 (both fits) are
marked at their (Δφ, mean colatitude). Each marker's ΔPF is computed from that
fit's **full** geometry (its own spots, i, u), not the canonical equal-θ cell, so
the markers carry the true result; they should fall on the saturated side (J0030)
and the tiling side (J0740) of the boundary. Caveat: Riley's J0740 spots straddle
the equator (77°/108°), so its mean-θ ≈ 93° plotting position is only a position —
its true ΔPF comes from the real, non-equal-colatitude pair.

Two panels share the axes and color scale, one per representative spacetime
(J0740-like compact/edge-on; J0030-like moderate/intermediate-i), so each star sits
on a map drawn at its own kind of spacetime.

Scope — which field configurations the heatmap covers. The colatitude × separation
plane is the natural space for all three canonical hot-spot arrangements, but the
*background heatmap* fixes both spots at the SAME colatitude θ, so it is a slice:

    * a CENTERED DIPOLE has antipodal poles — colatitudes θ and 180°−θ, separation 0.5
      — so it is NOT a single point on the equal-θ map (it coincides only at the
      equatorial corner θ = 90°, Δφ = 0.5). The dipole/antipodal case is the validation
      limit (reproduce Poutanen & Beloborodov 2006), handled separately, not the headline;
    * OFFSET-DIPOLE and MULTIPOLAR same-hemisphere pairs (J0030-like) live directly on
      this slice — the regime where the systematic actually varies, and the novel result.

The underlying machinery (`multi_spot_flux`) has no such restriction — it takes each
spot's own colatitude (e.g. J0740 Riley's 77°/108°), so the real-star markers carry the
true unequal-colatitude ΔPF even where the equal-θ background is only schematic.

Run from the repository root:  python3 scripts/phase_diagram.py
"""

from typing import NamedTuple, Sequence

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe

from mcrt import ExactBending, beaming_lookup, compute_profile, pulsed_fraction
from anchor_lib import (
    N_PHASE,
    SHAPE_TAU,
    Anchor,
    load_library,
    shape_tau_index,
    single_spot_eclipsed_fraction,
    two_spot_flux,
)
from j0030_anchor import RILEY as J0030_RILEY, MILLER as J0030_MILLER
from j0740_anchor import RILEY as J0740_RILEY, MILLER as J0740_MILLER

FIGURE_PATH = "data/phase_diagram.png"
RESULTS_PATH = "data/phase_diagram.npz"

# Sweep grid. Δφ ∈ [0, 0.5] — anti-phasing maxes out at half a cycle (0.5);
# beyond that the pair is the mirror image. θ spans the full 0–180° spherical
# colatitude. Δφ resolution finer than 1/N_PHASE is wasted (azimuth lands on the
# phase grid via np.roll), so ~100 columns is plenty.
DPHI_AXIS = np.linspace(0.0, 0.5, 101)
THETA_AXIS_DEG = np.linspace(0.0, 180.0, 91)
THETA_AXIS = np.deg2rad(THETA_AXIS_DEG)


class Regime(NamedTuple):
    """One representative spacetime panel, with the real fits that live there."""

    label: str
    inclination: float          # i, radians (the panel's fixed spacetime)
    compactness: float          # u = 2GM/Rc²
    stars: Sequence[Anchor]     # real fits marked on this panel (carry their own i, u)


REGIMES = (
    Regime(
        label="J0740-like  (compact, near edge-on:  u = 0.49,  i = 87.6°)",
        inclination=J0740_RILEY.inclination,
        compactness=J0740_RILEY.compactness,
        stars=(J0740_RILEY, J0740_MILLER),
    ),
    Regime(
        label="J0030-like  (moderate, intermediate i:  u = 0.31,  i = 53.9°)",
        inclination=J0030_RILEY.inclination,
        compactness=J0030_RILEY.compactness,
        stars=(J0030_RILEY, J0030_MILLER),
    ),
)


def delta_pf_grid(inclination, compactness, beaming_shape,
                  dphi_axis=DPHI_AXIS, theta_axis=THETA_AXIS):
    """ΔPF and the isotropic flux floor over the (Δφ, θ) grid for one spacetime.

    For each colatitude the single-spot isotropic and realistic base profiles are
    built once (the expensive part); each Δφ column is then the cheap additive sum
    of the base with its phase-rolled copy (the equal-weight, equal-colatitude
    canonical pair). Returns two ``(n_theta, n_dphi)`` arrays: ΔPF = PF_real −
    PF_iso, and F_min/F_max of the isotropic pair (0 ⇔ the pair eclipses to zero,
    i.e. PF is saturated).
    """
    n_t, n_d = len(theta_axis), len(dphi_axis)
    dpf = np.zeros((n_t, n_d))
    floor = np.zeros((n_t, n_d))
    for it, theta in enumerate(theta_axis):
        iso_base = compute_profile(inclination, theta, compactness).flux
        real_base = compute_profile(inclination, theta, compactness,
                                    beaming=beaming_shape).flux
        for jd, dphi in enumerate(dphi_axis):
            shift = int(round(dphi * N_PHASE)) % N_PHASE
            iso = iso_base + np.roll(iso_base, shift)
            real = real_base + np.roll(real_base, shift)
            dpf[it, jd] = pulsed_fraction(real) - pulsed_fraction(iso)
            hi = iso.max()
            floor[it, jd] = float(iso.min() / hi) if hi > 0 else 0.0
    return dpf, floor


def saturation_boundary(inclination, compactness, theta_axis=THETA_AXIS):
    """Analytic tiling boundary Δφ_crit(θ) = single-spot eclipse fraction.

    Below this Δφ the two anti-phased dark windows overlap and the combined pulse
    hits zero (PF saturates, ΔPF → 0); above it the spots tile the rotation and
    the systematic stays live in the pulsed fraction.
    """
    return np.array([
        single_spot_eclipsed_fraction(inclination, theta, compactness)
        for theta in theta_axis
    ])


def star_marker(anchor, beaming_shape):
    """Where a published fit sits on the map, and its true (full-geometry) ΔPF.

    Δφ is the fit's azimuthal spot separation and the plotting colatitude is the
    mean of its two spots, but the ΔPF is computed from the fit's *actual* spots,
    inclination, and compactness — not the canonical equal-θ cell — so the marker
    reports the real result regardless of how well the fit matches the canonical
    pair (it does not, exactly, for equator-straddling Riley J0740).

    The marker's ΔPF uses **exact Schwarzschild bending** (Gate G1, v0.9.9.1): the
    per-star result is the quantitative number the paper cites, so it carries the
    bending correction. The broad heatmap behind it stays on the linear map for
    cost (stated tolerance: the exact−linear ΔPF gap is ≤ 0.007 at the u = 0.49
    compact panel, smaller elsewhere — far below the marker/boundary resolution).
    """
    bending = ExactBending(anchor.compactness)
    iso = two_spot_flux(anchor.inclination, anchor.compactness, anchor.spots, None,
                        bending=bending)
    real = two_spot_flux(anchor.inclination, anchor.compactness, anchor.spots,
                         beaming_shape, bending=bending)
    # Fold the separation into [0, 0.5]: a pair at Δφ and at 1−Δφ are mirror images
    # with identical PF, and the sweep axis only spans the non-redundant half. Without
    # this, Miller J0740 (Δφ = 0.56) would silently fall off the right edge of the map.
    dphi = abs(anchor.spots[1].azimuth - anchor.spots[0].azimuth)
    dphi = min(dphi, 1.0 - dphi)
    theta_mean_deg = float(np.rad2deg(np.mean([s.colatitude for s in anchor.spots])))
    return {
        "label": anchor.label,
        "dphi": float(dphi),
        "theta_mean_deg": theta_mean_deg,
        "delta_pf_true": float(pulsed_fraction(real) - pulsed_fraction(iso)),
    }


def build(beaming_shape):
    """Compute the grid, analytic boundary, and star markers for every regime."""
    out = {}
    for regime in REGIMES:
        dpf, floor = delta_pf_grid(regime.inclination, regime.compactness, beaming_shape)
        out[regime.label] = {
            "delta_pf": dpf,
            "flux_floor": floor,
            "boundary": saturation_boundary(regime.inclination, regime.compactness),
            "markers": [star_marker(s, beaming_shape) for s in regime.stars],
        }
    return out


def plot_result(results, path=FIGURE_PATH):
    """Two (Δφ, θ) panels: ΔPF heatmap, analytic tiling boundary, real-star markers."""
    vmax = max(results[r.label]["delta_pf"].max() for r in REGIMES)
    fig, axes = plt.subplots(1, len(REGIMES), figsize=(15, 6.2), sharey=True)
    mesh = None
    for ax, regime in zip(axes, REGIMES):
        r = results[regime.label]

        # ΔPF heatmap. pcolormesh wants cell edges; the axes hold cell centers.
        mesh = ax.pcolormesh(DPHI_AXIS, THETA_AXIS_DEG, r["delta_pf"],
                             cmap="magma", vmin=0.0, vmax=vmax, shading="auto")

        # Analytic tiling boundary Δφ_crit(θ) = f_ecl(θ) (solid), and the
        # independent numerical flux-floor→0 contour (dashed) — they should agree.
        ax.plot(r["boundary"], THETA_AXIS_DEG, color="#39d0ff", lw=2.4,
                label=r"analytic boundary $\Delta\phi_{\rm crit}=f_{\rm ecl}(\theta)$")
        ax.contour(DPHI_AXIS, THETA_AXIS_DEG, r["flux_floor"],
                   levels=[1e-3], colors="white", linewidths=1.2, linestyles="--")

        ax.axhline(90.0, color="#7f8c8d", lw=0.7, alpha=0.6)  # equator reference
        for i, m in enumerate(r["markers"]):
            ax.scatter(m["dphi"], m["theta_mean_deg"], s=170, marker="*",
                       facecolor=plt.cm.magma(m["delta_pf_true"] / vmax),
                       edgecolor="white", linewidth=1.4, zorder=5)
            # Labels flip to the left for right-side markers so they never clip the
            # panel edge, and stagger up/down by marker index so co-located fits (the
            # two J0740 teams land on the same cell) do not overprint. A dark stroke
            # keeps white text legible over bright cells.
            near_right = m["dphi"] > 0.30
            dy = 9 if i % 2 == 0 else -17
            ax.annotate(f"{m['label'].split(' (')[0]}  ΔPF={m['delta_pf_true']:+.2f}",
                        (m["dphi"], m["theta_mean_deg"]), textcoords="offset points",
                        xytext=(-10 if near_right else 10, dy),
                        ha="right" if near_right else "left", fontsize=7.5,
                        color="white", zorder=6,
                        path_effects=[pe.withStroke(linewidth=2.2, foreground="black")])
        ax.set_xlabel(r"Azimuthal spot separation  $\Delta\phi$  (cycles)")
        ax.set_title(regime.label, fontsize=9.5)
        ax.set_xlim(0.0, 0.5)
        ax.legend(loc="lower right", fontsize=7.5, framealpha=0.85)
    axes[0].set_ylabel(r"Spot colatitude  $\theta$  (deg)")

    cbar = fig.colorbar(mesh, ax=axes, fraction=0.046, pad=0.02)
    cbar.set_label(r"$\Delta\mathrm{PF} = \mathrm{PF}_{\rm real} - \mathrm{PF}_{\rm iso}$"
                   f"  (at τ = {SHAPE_TAU:g})")
    fig.suptitle("Where the beaming systematic lands: tiling (right of boundary) keeps it "
                 "live in PF;\nfailing to tile (left) saturates PF and the systematic hides "
                 "in waveform shape", fontsize=11)
    fig.savefig(path, dpi=130, bbox_inches="tight")
    print(f"✓ figure saved to {path}")


def save_results(results, path=RESULTS_PATH):
    payload = {"dphi_axis": DPHI_AXIS, "theta_axis_deg": THETA_AXIS_DEG,
               "shape_tau": SHAPE_TAU}
    for regime in REGIMES:
        key = "j0740" if "J0740" in regime.label else "j0030"
        r = results[regime.label]
        payload[f"{key}_delta_pf"] = r["delta_pf"]
        payload[f"{key}_flux_floor"] = r["flux_floor"]
        payload[f"{key}_boundary"] = r["boundary"]
        for m in r["markers"]:
            mkey = "riley" if "Riley" in m["label"] else "miller"
            payload[f"{key}_{mkey}_dphi"] = m["dphi"]
            payload[f"{key}_{mkey}_theta_deg"] = m["theta_mean_deg"]
            payload[f"{key}_{mkey}_delta_pf"] = m["delta_pf_true"]
    np.savez(path, **payload)
    print(f"✓ results saved to {path}")


def print_summary(results):
    print("\n" + "=" * 74)
    print("Geometry phase diagram — tiling vs. saturation in (Δφ, θ) space")
    print("=" * 74)
    for regime in REGIMES:
        r = results[regime.label]
        print(f"\n{regime.label}")
        for m in r["markers"]:
            # f_ecl at the marker's mean colatitude on THIS panel's spacetime, for
            # the side-of-boundary read (Δφ vs Δφ_crit).
            f_ecl = single_spot_eclipsed_fraction(
                regime.inclination, np.deg2rad(m["theta_mean_deg"]), regime.compactness)
            side = "TILING → PF-live" if m["dphi"] > f_ecl else "no-tile → PF-saturated"
            print(f"  {m['label'].split(' (')[0]:<26} "
                  f"Δφ = {m['dphi']:.2f}  vs  Δφ_crit(θ̄) = {f_ecl:.2f}   "
                  f"→ {side};  true ΔPF = {m['delta_pf_true']:+.2f}")
    print("\nReading: a star is PF-visible (live ΔPF) when its azimuthal separation Δφ\n"
          "exceeds the single-spot eclipse fraction f_ecl(θ) — i.e. when the two spots\n"
          "tile the rotation. J0740's wide-separation spots clear the boundary; J0030's\n"
          "same-hemisphere spots fall short, so PF saturates and the systematic retreats\n"
          "into waveform shape. Tiling — not single-spot eclipse — is the discriminator.")


def run():
    tau_values, mu_centers, intensity_by_tau = load_library()
    shape_idx = shape_tau_index(tau_values)  # raises if SHAPE_TAU left the grid
    beaming_shape = beaming_lookup(mu_centers, intensity_by_tau[shape_idx])
    return build(beaming_shape)


if __name__ == "__main__":
    results = run()
    save_results(results)
    plot_result(results)
    print_summary(results)
