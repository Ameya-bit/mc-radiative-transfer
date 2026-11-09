"""Journal-grade redesigns of the three load-bearing paper figures (APSS draft).

Design contract (see README.md in this folder for sources):
  * one message per figure, stated in the caption, not in an in-figure title;
  * no internal codenames (Track D, C2->C4), no literal underscores (F_max);
  * Okabe-Ito colorblind-safe palette, consistent encoding across figures
    (Riley = blue, Miller = vermillion, isotropic = neutral gray dashed);
  * panel letters (a)/(b)/(c); direct labeling preferred over legends;
  * sized for the Springer sn-jnl text block (~6.6 in): full-width 170 mm,
    fonts >= 7 pt at print size; vector PDF + 600 dpi PNG for review.

All curves are recomputed with the repo's own verified pipeline (mcrt +
scripts/anchor_lib); all scalars come from the committed npz evidence base.
Run from anywhere:  python make_updated_figs.py
"""

import sys
from pathlib import Path

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[4]
sys.path[:0] = [str(ROOT / "src"), str(ROOT / "scripts")]

from mcrt import compute_profile                      # noqa: E402
from mcrt.bending import ExactBending                 # noqa: E402
from mcrt.beaming import beaming_lookup               # noqa: E402
from anchor_lib import multi_spot_flux, shape_tau_index  # noqa: E402
from j0740_anchor import RILEY as J0740_RILEY, MILLER as J0740_MILLER  # noqa: E402
from j0030_anchor import RILEY as J0030_RILEY, MILLER as J0030_MILLER  # noqa: E402

DATA = ROOT / "data"
OUT = Path(__file__).resolve().parent

# ---- palette (Okabe-Ito + gold accent) ------------------------------------
BLUE = "#0072B2"        # Riley fits
VERMILLION = "#D55E00"  # Miller fits
GRAY = "#5A5A5A"        # isotropic reference
LIGHT = "#BBBBBB"
GOLD = "#FFC845"        # accent text on dark maps (mentor: gold on Fig. 5)

# Mentor feedback applied: labels sized to read at print without zooming
plt.rcParams.update({
    "font.size": 9,
    "axes.labelsize": 9.5,
    "axes.linewidth": 0.7,
    "xtick.labelsize": 8.5,
    "ytick.labelsize": 8.5,
    "legend.fontsize": 8,
    "legend.frameon": False,
    "xtick.direction": "in",
    "ytick.direction": "in",
    "xtick.top": True,
    "ytick.right": True,
    "savefig.dpi": 600,
    "pdf.fonttype": 42,
})

MM = 1 / 25.4  # mm -> inch


def panel_tag(ax, letter, x=0.02, y=0.96):
    ax.text(x, y, f"({letter})", transform=ax.transAxes, fontsize=9,
            fontweight="bold", va="top")


def pulsed_fraction(flux):
    lo, hi = float(flux.min()), float(flux.max())
    return (hi - lo) / (hi + lo) if (hi + lo) > 0 else 0.0


def load_library():
    d = np.load(DATA / "beaming_library.npz")
    return d["tau_values"], d["mu_centers"], d["intensity_by_tau"]


# ===========================================================================
# Figure 1 — J0740: realistic beaming sharpens the pulse, ΔPF is the bias
# ===========================================================================

def fig_j0740():
    tau_values, mu_centers, intensity_by_tau = load_library()
    ti = shape_tau_index(tau_values)  # tau = 10 production library
    beaming = beaming_lookup(mu_centers, intensity_by_tau[ti])

    anchors = [("Riley 2021", J0740_RILEY, BLUE),
               ("Miller 2021", J0740_MILLER, VERMILLION)]

    j = np.load(DATA / "j0740_anchor.npz")
    a3 = np.load(DATA / "a3_seed_errors.npz")

    fig, axes = plt.subplots(
        1, 3, figsize=(170 * MM, 52 * MM), constrained_layout=True)

    phase = np.linspace(0, 1, 1024, endpoint=False)
    for ax, (name, anchor, color) in zip(axes[:2], anchors):
        bending = ExactBending(anchor.compactness)
        iso = multi_spot_flux(anchor.inclination, anchor.compactness,
                              anchor.spots, None, bending=bending)
        real = multi_spot_flux(anchor.inclination, anchor.compactness,
                               anchor.spots, beaming, bending=bending)
        # peak-normalize each curve: overall constants cancel from PF and
        # pulse shapes (Eq. flux), so only the modulation is physical
        ax.plot(phase, iso / iso.max(), ls="--", lw=1.5, color=GRAY,
                label="isotropic")
        ax.plot(phase, real / real.max(), lw=1.7, color=color,
                label=r"realistic $I(\mu)$")
        ax.fill_between(phase, iso / iso.max(), real / real.max(),
                        color=color, alpha=0.12, lw=0)
        dpf = pulsed_fraction(real) - pulsed_fraction(iso)
        ax.text(0.97, 0.04,
                name + "\n" + rf"$\Delta\mathrm{{PF}} = {dpf:+.3f}$",
                transform=ax.transAxes, ha="right", va="bottom", fontsize=7.5,
                color="0.15",
                bbox=dict(fc="white", ec="none", alpha=0.75, pad=1.5))
        ax.set_xlim(0, 1)
        lo = min((iso / iso.max()).min(), (real / real.max()).min())
        # headroom so the legend row clears both the curves and the (a)/(b)
        # panel tags (the tags previously collided with "isotropic")
        ax.set_ylim(lo - 0.05, 1.28)
        ax.set_xlabel(r"Rotational phase $\varphi/2\pi$")
        ax.legend(loc="upper right", ncols=2, fontsize=7, handlelength=1.2,
                  columnspacing=0.7, handletextpad=0.4, borderaxespad=0.2)
    axes[0].set_ylabel(r"Flux $F(\varphi)\,/\,F_{\max}$")
    panel_tag(axes[0], "a")
    panel_tag(axes[1], "b")

    # (c) ΔPF vs optical depth, both fits, seed-scatter error bars
    ax = axes[2]
    for key, color, name, marker in (
            ("riley", BLUE, "Riley 2021", "o"),
            ("miller", VERMILLION, "Miller 2021", "s")):
        ax.errorbar(j["tau_values"], j[f"{key}_delta_pf"],
                    yerr=a3[f"j0740_{key}_dpf_std"],
                    color=color, marker=marker, ms=3.5, lw=1.1, capsize=2,
                    label=name)
    ax.set_xscale("log")
    ax.axhline(0, color=LIGHT, lw=0.7, zorder=0)
    ax.set_xlabel(r"Slab optical depth $\tau$")
    ax.set_ylabel(r"$\Delta\mathrm{PF} = \mathrm{PF}_{\rm real} - \mathrm{PF}_{\rm iso}$")
    ax.legend(loc="lower right")
    panel_tag(ax, "c")

    fig.savefig(OUT / "fig_j0740_updated.pdf")
    fig.savefig(OUT / "fig_j0740_updated.png")
    plt.close(fig)


# ===========================================================================
# Figure 2 — routing by spin: ΔPF collapses, waveform shape survives
# ===========================================================================

def fig_routing():
    c3 = np.load(DATA / "c3_band_doppler.npz")
    c4 = np.load(DATA / "c4_caveat_audit.npz")
    j = np.load(DATA / "j0740_anchor.npz")
    a3 = np.load(DATA / "a3_seed_errors.npz")

    def s(key):
        return float(np.atleast_1d(c3[key])[0])

    def s4(key):
        return float(np.atleast_1d(c4[key])[0])

    # Rung order mirrors Table 2 of the manuscript (top = static baseline).
    rungs = [
        "Static (0 Hz), bolometric",
        r"346.5 Hz, bolometric $\delta^4$",
        "346.5 Hz, NICER band (photons)",
        r"$+$ light-travel delay, $\delta^4$",
        r"$+$ light-travel delay, band",
    ]
    riley = [s("j0740_riley_0hz_m"), s("j0740_riley_d4_m"),
             s("j0740_riley_band_ph_m"),
             s("j0740_riley_d4_m") + s4("j0740_riley_delay_d4_shift"),
             s("j0740_riley_band_ph_m") + s4("j0740_riley_delay_band_ph_shift")]
    miller = [s("j0740_miller_0hz_m"), s("j0740_miller_d4_m"),
              s("j0740_miller_band_ph_m"),
              s("j0740_miller_d4_m") + s4("j0740_miller_delay_d4_shift"),
              s("j0740_miller_band_ph_m") + s4("j0740_miller_delay_band_ph_shift")]
    riley_err = [s("j0740_riley_0hz_s"), s("j0740_riley_d4_s"),
                 s("j0740_riley_band_ph_s"), s4("j0740_riley_delay_d4_sigma"),
                 s4("j0740_riley_delay_band_ph_sigma")]
    miller_err = [s("j0740_miller_0hz_s"), s("j0740_miller_d4_s"),
                  s("j0740_miller_band_ph_s"), s4("j0740_miller_delay_d4_sigma"),
                  s4("j0740_miller_delay_band_ph_sigma")]

    # Waveform-shape rms on the rungs where it was measured
    # (static from the anchor sweep at tau=10; spinning from Tracks C3/C4).
    ti = shape_tau_index(j["tau_values"])
    shape_riley = [float(j["riley_shape_rms"][ti]), np.nan,
                   s("j0740_riley_band_shape_rms"), np.nan,
                   float(np.atleast_1d(c4["j0740_riley_warp_shape_band_ph_rms"])[0])]
    shape_miller = [float(j["miller_shape_rms"][ti]), np.nan,
                    s("j0740_miller_band_shape_rms"), np.nan,
                    float(np.atleast_1d(c4["j0740_miller_warp_shape_band_ph_rms"])[0])]

    ypos = np.arange(len(rungs))[::-1]
    off = 0.16

    fig, (ax1, ax2) = plt.subplots(
        1, 2, figsize=(170 * MM, 62 * MM), sharey=True,
        gridspec_kw={"width_ratios": [1.45, 1.0]}, constrained_layout=True)

    for vals, errs, color, name, marker, dy in (
            (riley, riley_err, BLUE, "Riley 2021", "o", +off),
            (miller, miller_err, VERMILLION, "Miller 2021", "s", -off)):
        ax1.errorbar(vals, ypos + dy, xerr=errs, color=color, marker=marker,
                     ms=4, ls="none", capsize=2, label=name)
        for v, yv in zip(vals, ypos + dy):
            # labels sit beside markers (never above/below) so adjacent
            # rungs cannot collide; flip side near the right axis edge
            if v > 0.15:
                kw = dict(xytext=(-7, -2), ha="right")
            else:
                kw = dict(xytext=(7, -2), ha="left")
            ax1.annotate(f"{v:+.3f}", (v, yv), textcoords="offset points",
                         fontsize=6.3, color=color, **kw)
    ax1.axvline(0, color=LIGHT, lw=0.7, zorder=0)
    ax1.set_yticks(ypos, rungs)
    ax1.set_ylim(-0.8, len(rungs) - 1 + 0.8)
    ax1.set_xlim(-0.01, 0.225)
    ax1.set_xlabel(r"Pulsed-fraction bias $\Delta\mathrm{PF}$")
    ax1.legend(loc="center right")
    panel_tag(ax1, "a", x=0.03)

    for vals, color, name, marker, dy in (
            (shape_riley, BLUE, "Riley 2021", "o", +off),
            (shape_miller, VERMILLION, "Miller 2021", "s", -off)):
        m = ~np.isnan(vals)
        ax2.plot(np.asarray(vals)[m], (ypos + dy)[m], color=color,
                 marker=marker, ms=4, ls="none")
    ax2.axvline(0, color=LIGHT, lw=0.7, zorder=0)
    ax2.set_xlim(0, 0.16)
    ax2.set_xlabel("Waveform-shape difference (rms)")
    ax2.annotate("shape bias\npersists", xy=(0.105, ypos[-1] - 0.45),
                 fontsize=7, color="0.25", ha="center")
    panel_tag(ax2, "b", x=0.05)

    fig.savefig(OUT / "fig_routing_updated.pdf")
    fig.savefig(OUT / "fig_routing_updated.png")
    plt.close(fig)


# ===========================================================================
# Figure 3 — phase diagram: where the bias lives, saturated region marked
# ===========================================================================

def fig_phase_diagram():
    import matplotlib.patheffects as pe
    halo = [pe.withStroke(linewidth=1.6, foreground="k")]

    d = np.load(DATA / "phase_diagram.npz")
    dphi = d["dphi_axis"]
    theta = d["theta_axis_deg"]

    # anchor label offsets (points) chosen to avoid overlap: the two J0740
    # anchors nearly coincide, so one label goes above, one below
    panels = [
        ("a", "j0740", r"J0740-like  ($u=0.49$, $i=87.6^\circ$)",
         [("Riley 2021", d["j0740_riley_dphi"], d["j0740_riley_theta_deg"],
           (-8, 10), "right"),
          ("Miller 2021", d["j0740_miller_dphi"], d["j0740_miller_theta_deg"],
           (-8, -16), "right")]),
        ("b", "j0030", r"J0030-like  ($u=0.31$, $i=53.9^\circ$)",
         [("Riley 2019", d["j0030_riley_dphi"], d["j0030_riley_theta_deg"],
           (-8, 8), "right"),
          ("Miller 2019", d["j0030_miller_dphi"], d["j0030_miller_theta_deg"],
           (-8, -16), "right")]),
    ]

    vmax = max(d["j0740_delta_pf"].max(), d["j0030_delta_pf"].max())

    fig, axes = plt.subplots(
        1, 2, figsize=(170 * MM, 72 * MM), sharey=True, constrained_layout=True)

    for ax, (letter, key, label, anchors) in zip(axes, panels):
        dpf = d[f"{key}_delta_pf"]
        floor = d[f"{key}_flux_floor"]
        im = ax.pcolormesh(dphi, theta, dpf, cmap="magma", vmin=0, vmax=vmax,
                           shading="auto", rasterized=True)
        # saturated regime boundary: summed flux touches zero -> PF pinned at 1
        ax.contour(dphi, theta, floor, levels=[1e-9], colors="w",
                   linewidths=0.9, linestyles="--")
        # analytic tiling boundary under both bending maps (R4): the linear map
        # is the one the criterion is derived with; the exact map bounds its
        # placement error
        r1 = np.load(DATA / "r1_revision_checks.npz")
        for src, style, lw, lab in (
                (d[f"{key}_boundary"], "-", 1.6, "linear"),
                (r1[f"r4_{key}_boundary_exact"], (0, (4, 2)), 1.4, "exact")):
            b = np.asarray(src, dtype=float)
            inside = b <= dphi.max()  # sentinel 1.0 = no boundary at this theta
            ax.plot(np.where(inside, b, np.nan), theta, color="#56B4E9",
                    ls=style, lw=lw, label=f"tiling boundary ({lab} bending)")
        if letter == "a":
            ax.legend(loc="lower right", fontsize=6.8, labelcolor=GOLD,
                      framealpha=0)
        for name, adphi, atheta, offset, ha in anchors:
            ax.plot(adphi, atheta, marker="*", ms=12, mfc=GOLD, mec="k",
                    mew=0.7, ls="none")
            ax.annotate(name, (adphi, atheta), textcoords="offset points",
                        xytext=offset, ha=ha, fontsize=7.5, color=GOLD,
                        path_effects=halo)
        ax.set_xlim(dphi.min(), dphi.max())
        ax.set_xlabel(r"Azimuthal spot separation $\Delta\varphi$ (cycles)")
        ax.set_title(label, fontsize=9)
        ax.text(0.03, 0.975, f"({letter})", transform=ax.transAxes, fontsize=10,
                fontweight="bold", va="top", color=GOLD, path_effects=halo)
    axes[0].set_ylabel(r"Spot colatitude $\theta_s$ (deg)")

    # regime annotations (panel b has the clearest split); gold per mentor note
    axes[1].text(0.13, 95, "PF saturated\n(bias moves to\nwaveform shape)",
                 color=GOLD, fontsize=8, ha="center", path_effects=halo)
    axes[1].text(0.33, 25, "PF sensitive\n(bias visible in $\\Delta$PF)",
                 color=GOLD, fontsize=8, ha="center", path_effects=halo)

    cbar = fig.colorbar(im, ax=axes, shrink=0.92, pad=0.015)
    cbar.set_label(r"$\Delta\mathrm{PF}$ at $\tau = 10$")

    fig.savefig(OUT / "fig_phase_diagram_updated.pdf")
    fig.savefig(OUT / "fig_phase_diagram_updated.png")
    plt.close(fig)


# ===========================================================================
# Figure 4 — validation gallery: H(mu), SD1a static, SD1c Doppler
# ===========================================================================

def fig_validation():
    import os
    os.chdir(ROOT)  # code_comparison helpers use repo-relative data paths
    from mcrt.theory import chandrasekhar_h, eddington_limb_darkening
    from code_comparison import (load_reference, our_profile_at,
                                 COMPACTNESS as SD_U)
    from c1_doppler_validate import monochromatic_flux

    d = np.load(DATA / "d_isotropic_h.npz")
    mu = d["mu_centers"]

    fig = plt.figure(figsize=(170 * MM, 68 * MM), constrained_layout=True)
    gs = fig.add_gridspec(2, 3, height_ratios=[2.6, 1.1])
    tops = [fig.add_subplot(gs[0, k]) for k in range(3)]
    bots = [fig.add_subplot(gs[1, k], sharex=tops[k]) for k in range(3)]

    # ---- (a) Monte Carlo vs Chandrasekhar H(mu) ---------------------------
    ax, axr = tops[0], bots[0]
    mu_fine = np.linspace(1e-3, 1, 300)
    h_fine = chandrasekhar_h(mu_fine)
    ax.plot(mu_fine, h_fine / chandrasekhar_h(np.array([1.0]))[0], color="k",
            lw=1.5, label=r"$H(\mu)$ exact")
    ax.plot(mu_fine, eddington_limb_darkening(mu_fine)
            / eddington_limb_darkening(np.array([1.0]))[0], color=GRAY,
            ls=":", lw=1.3, label=r"Eddington $1+\frac{3}{2}\mu$")
    ax.errorbar(mu, d["isotropic_i_mean"], yerr=d["isotropic_sem"], ls="none",
                marker="o", ms=3.4, color=BLUE, lw=1,
                label="MC isotropic")
    ax.errorbar(mu, d["thomson_i_mean"], yerr=d["thomson_sem"], ls="none",
                marker="s", ms=3.2, color=VERMILLION, lw=1,
                label="MC Thomson")
    ax.set_ylabel(r"$I(\mu)\,/\,I(1)$")
    ax.legend(loc="lower right", fontsize=7, handlelength=1.4,
              labelspacing=0.25)
    panel_tag(ax, "a")
    chi2 = float(np.atleast_1d(d["isotropic_reduced_chi2"])[0])
    ax.text(0.04, 0.66, rf"$\chi^2/\mathrm{{dof}} = {chi2:.2f}$ (isotropic)",
            transform=ax.transAxes, fontsize=7.5, color="0.25")

    axr.axhspan(-2, 2, color="#EAF3EA", zorder=0)
    axr.axhline(0, color="k", lw=0.6)
    axr.plot(mu, d["isotropic_residual_flux"], "o", ms=3, color=BLUE)
    axr.plot(mu, d["thomson_residual_flux"], "s", ms=2.8, color=VERMILLION)
    axr.set_ylim(-3.2, 3.2)
    axr.set_xlabel(r"$\mu = \cos\theta$")
    axr.set_ylabel(r"$\Delta/\sigma$ (flux)")

    # ---- (b) static SD1a benchmark ---------------------------------------
    ax, axr = tops[1], bots[1]
    phase, f_ref = load_reference(str(DATA / "l26_reference/SD1a_test_IM.txt"))
    ours = our_profile_at(phase, bending=ExactBending(SD_U))
    ref_n, ours_n = f_ref / f_ref.max(), ours / ours.max()
    ax.plot(phase, ref_n, color="k", lw=1.5, label="reference (SD1a)")
    ax.plot(phase[::3], ours_n[::3], "o", ms=3.4, color=BLUE,
            label="this work")
    ax.set_ylabel(r"$F(\varphi)\,/\,F_{\max}$")
    ax.legend(loc="lower center", bbox_to_anchor=(0.5, 0.10), fontsize=7,
              handlelength=1.4, labelspacing=0.25)
    panel_tag(ax, "b")
    resid = 100 * (ours_n - ref_n)
    ax.text(0.5, 0.42, rf"max $|\Delta| = {np.abs(resid).max():.2f}\%$",
            transform=ax.transAxes, fontsize=7.5, color="0.25", ha="center")

    axr.axhline(0, color="k", lw=0.6)
    axr.plot(phase, resid, color=BLUE, lw=1.1)
    axr.set_xlabel(r"Rotational phase $\varphi/2\pi$")
    axr.set_ylabel(r"$\Delta$ (% peak)")

    # ---- (c) rotating SD1c benchmark -------------------------------------
    ax, axr = tops[2], bots[2]
    ref_c = np.loadtxt(DATA / "l26_reference/SD1c_test_IM.txt")
    phase_c, f_refc = ref_c[:, 0], ref_c[:, 1]
    ours_c = monochromatic_flux(phase_c, 200.0, with_time_delay=True)
    refc_n, oursc_n = f_refc / f_refc.max(), ours_c / ours_c.max()
    ax.plot(phase, ref_n, color=LIGHT, lw=1.2, label="static (SD1a)")
    ax.plot(phase_c, refc_n, color="k", lw=1.5, label="reference (SD1c)")
    ax.plot(phase_c[::3], oursc_n[::3], "o", ms=3.4, color=VERMILLION,
            label="this work, 200 Hz")
    ax.legend(loc="lower center", bbox_to_anchor=(0.5, 0.08), fontsize=7,
              handlelength=1.4, labelspacing=0.25)
    panel_tag(ax, "c")
    resid_c = 100 * (oursc_n - refc_n)
    ax.text(0.5, 0.50, rf"max $|\Delta| = {np.abs(resid_c).max():.2f}\%$",
            transform=ax.transAxes, fontsize=7.5, color="0.25", ha="center")

    axr.axhline(0, color="k", lw=0.6)
    axr.plot(phase_c, resid_c, color=VERMILLION, lw=1.1)
    axr.set_xlabel(r"Rotational phase $\varphi/2\pi$")

    fig.savefig(OUT / "fig_validation_updated.pdf")
    fig.savefig(OUT / "fig_validation_updated.png")
    plt.close(fig)


# ===========================================================================
# Figure 5 — J0030: PF saturates, the bias moves into the waveform shape
# ===========================================================================

def fig_j0030():
    tau_values, mu_centers, intensity_by_tau = load_library()
    ti = shape_tau_index(tau_values)
    beaming = beaming_lookup(mu_centers, intensity_by_tau[ti])

    # linear bending map (None), matching the published J0030 anchor numbers
    anchors = [("Riley 2019", J0030_RILEY, BLUE),
               ("Miller 2019", J0030_MILLER, VERMILLION)]

    j = np.load(DATA / "j0030_anchor.npz")
    a3 = np.load(DATA / "a3_seed_errors.npz")

    fig, axes = plt.subplots(
        1, 3, figsize=(170 * MM, 55 * MM), constrained_layout=True)

    phase = np.linspace(0, 1, 1024, endpoint=False)
    for ax, (name, anchor, color) in zip(axes[:2], anchors):
        iso = multi_spot_flux(anchor.inclination, anchor.compactness,
                              anchor.spots, None)
        real = multi_spot_flux(anchor.inclination, anchor.compactness,
                               anchor.spots, beaming)
        rms = float(np.sqrt(np.mean(
            (real / real.max() - iso / iso.max()) ** 2)))
        ax.plot(phase, iso / iso.max(), ls="--", lw=1.5, color=GRAY,
                label="isotropic")
        ax.plot(phase, real / real.max(), lw=1.7, color=color,
                label=r"realistic $I(\mu)$")
        ax.fill_between(phase, iso / iso.max(), real / real.max(),
                        color=color, alpha=0.12, lw=0)
        ax.axhline(0, color=LIGHT, lw=0.8, zorder=0)
        ax.text(0.97, 0.68,
                name + "\n" + r"$\mathrm{PF} = 1$ (both)"
                + "\n" + rf"shape rms $= {rms:.3f}$",
                transform=ax.transAxes, ha="right", va="top", fontsize=7.5,
                color="0.15",
                bbox=dict(fc="white", ec="none", alpha=0.75, pad=1.5))
        ax.set_xlim(0, 1)
        ax.set_ylim(-0.06, 1.2)
        ax.set_xlabel(r"Rotational phase $\varphi/2\pi$")
        ax.legend(loc="upper left", handlelength=1.4, fontsize=7,
                  labelspacing=0.25)
    axes[0].set_ylabel(r"Flux $F(\varphi)\,/\,F_{\max}$")
    panel_tag(axes[0], "a")
    panel_tag(axes[1], "b")

    ax = axes[2]
    for key, color, name, marker in (
            ("riley", BLUE, "Riley 2019", "o"),
            ("miller", VERMILLION, "Miller 2019", "s")):
        ax.errorbar(j["tau_values"], j[f"{key}_shape_rms"],
                    yerr=a3[f"j0030_{key}_shape_std"],
                    color=color, marker=marker, ms=4, lw=1.3, capsize=2,
                    label=name)
    ax.set_xscale("log")
    ax.axhline(0, color=LIGHT, lw=0.7, zorder=0)
    ax.set_ylim(0, 0.08)
    ax.set_xlabel(r"Slab optical depth $\tau$")
    ax.set_ylabel("Waveform-shape rms")
    ax.legend(loc="lower right", fontsize=7)
    panel_tag(ax, "c")

    fig.savefig(OUT / "fig_j0030_updated.pdf")
    fig.savefig(OUT / "fig_j0030_updated.png")
    plt.close(fig)


if __name__ == "__main__":
    fig_j0740()
    fig_routing()
    fig_phase_diagram()
    fig_validation()
    fig_j0030()
    print("wrote", *sorted(p.name for p in OUT.glob("fig_*_updated.*")))
