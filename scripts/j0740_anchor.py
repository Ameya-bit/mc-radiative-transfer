"""Second real-star anchor: the beaming swap at PSR J0740+6620, where it does NOT eclipse.

v0.9.1 planted the isotropic→realistic beaming swap at PSR J0030+0451 and found the
pulsed-fraction systematic *saturated* to ΔPF ≈ 0: J0030's two spots hug the same pole,
dive behind the star, and pin PF = 1 for every beaming, so the systematic could only
show up in the waveform *shape*. The obvious question that leaves open is whether a real
star ever lets the systematic land back in the pulsed fraction. This script answers it
at the other end of the NICER mass–radius sample.

    Riley et al. (2021, ApJL 918 L27)   — X-PSI / Amsterdam, ST-U model
    Miller et al. (2021, ApJL 918 L28)  — Illinois–Maryland, two-circle

The result — and why J0740 is the complement of J0030. Riley's fit places J0740's two
hot spots in *opposite* hemispheres (colatitudes 77° and 108°) viewed almost edge-on
(i ≈ 87.6°). With its extreme compactness (u ≈ 0.494) light bending lets us see far
enough around the back that **neither spot ever sets** — each only grazes the limb
(μ_min ≈ 0.005). And because the spots are ~half a cycle apart in azimuth (Δφ ≈ 0.44),
they *tile* the rotation: when one is at the limb the other faces us, so the summed
flux never falls near zero (floor F_min/F_max ≈ 0.70). The isotropic pulsed fraction is
therefore low and far from saturation (PF ≈ 0.18), leaving ample headroom. Swapping in
the realistic beaming raises it directly: **ΔPF = +0.137 ± 0.003 at τ ≈ 10** (production
library, exact Schwarzschild bending), the same sign, size, and τ-dependence v0.9.0 found
on invented geometries — now on a published star and *visible in the pulsed fraction itself*.

Both teams agree (the cross-check anchor). Miller's independent fit agrees with Riley on
M, R, inclination, anti-phasing, and the two-circular-spot model, but places the spots
essentially *on the equator* (colatitude ≈ 92°). There each spot's *center* does dip
behind the star for ~21% of the cycle — yet because the two spots are anti-phased, the
eclipses fall at opposite phases and the *combined* pulse still never reaches zero
(floor ≈ 0.63). So PF stays unsaturated for Miller too, and the swap lands an even larger
**ΔPF = +0.195 ± 0.005 at τ ≈ 10**. The headline is therefore robust to the choice of team:
at J0740 the beaming systematic is a +0.137…+0.195 pulsed-fraction effect either way.

The discriminator, made precise. The thing that flips J0030's verdict is *not* whether a
single spot eclipses — Miller's J0740 spots eclipse 21% each and still leave a live ΔPF.
It is whether the two spots **tile** the rotation. J0030's spots hug the *same* far
hemisphere and each is hidden ~45% of the cycle, so even anti-phased their visible
windows leave a dark gap: the summed flux touches zero, PF pins at 1, ΔPF ≈ 0, and the
systematic retreats into waveform shape. J0740's spots sit in opposite regions and are
shallow enough that anti-phasing fully covers the cycle: the pulse stays off zero, PF
keeps its headroom, and the systematic shows up where NICER reads it — in the pulsed
fraction. PF-visible at J0740, shape-only at J0030: tiling is the knob.

Built with no change to the verified ``mcrt.pulse`` core: each spot is a
:func:`mcrt.compute_profile` call, a spot at a different longitude is that profile
rolled in phase, and light is additive — all via :mod:`anchor_lib`, shared with the
J0030 script so the two cannot drift apart.

Run from the repository root:  python scripts/j0740_anchor.py
"""

import numpy as np
import matplotlib.pyplot as plt

from mcrt import ExactBending, beaming_lookup, pulsed_fraction
from anchor_lib import (
    BACKGROUND_FRACS,
    N_PHASE,
    SHAPE_TAU,
    Anchor,
    Spot,
    delta_pf_vs_background,
    load_library,
    shape_tau_index,
    single_spot_eclipsed_fraction,
    single_spot_min_visible_mu,
    two_spot_flux,
    waveform_shape_change,
)

FIGURE_PATH = "data/pulse_profile_j0740.png"
RESULTS_PATH = "data/j0740_anchor.npz"


# --- Published geometries (verified against the papers' tables) ----------------
#
# Riley 2021 ST-U, headline NICER×XMM table (medians): M = 2.072 M⊙, R_eq = 12.39 km;
# compactness GM/Rc² = 0.247 ⇒ u = 0.494. cos i = 0.0424 ⇒ i = 1.5284 rad (87.57°).
# ST-U = two independent, single-temperature CIRCULAR caps (no crescent/masking), so
# area = π sin²ζ is exact and the point reduction is clean. The two caps:
#   spot p: Θp = 1.35 rad (77.3°), ζp = 0.147 rad, log₁₀T = 5.988
#   spot s: Θs = 1.89 rad (108.3°), ζs = 0.146 rad, log₁₀T = 5.992
# Phases are reported in offset Earth / Earth-antipode frames and flagged bimodal; the
# ML-point separation is Δφ ≈ 0.442 cyc (anti-phased, opposite hemispheres). Bolometric
# weight ∝ area(sin²ζ) × T⁴; the two caps are near-identical in size and temperature, so
# the weighting barely matters (robustness-checked against equal weights).
_R_WP = np.sin(0.147) ** 2 * (10 ** 5.988) ** 4
_R_WS = np.sin(0.146) ** 2 * (10 ** 5.992) ** 4
RILEY = Anchor(
    label="Riley 2021 (X-PSI ST-U)",
    inclination=1.5284,
    compactness=0.494,
    spots=(
        Spot(colatitude=1.35, azimuth=0.000, weight=_R_WP),  # northern cap
        Spot(colatitude=1.89, azimuth=0.442, weight=_R_WS),  # southern cap, anti-phased
    ),
    note="weight ∝ area(sin²ζ) × T⁴; off-equator spots, no eclipse",
)

# Miller 2021 two-circle, headline NICER+XMM table (Table 7 medians): M = 2.062 M⊙,
# R_e = 13.713 km; GM/Rc² = 0.222 ⇒ u = 0.444. θobs = 1.527 rad (87.49°). Two uniform
# circular spots, essentially ON the equator — and that is the divergence from Riley:
#   spot 1: θc1 = 1.600 rad (91.67°), Δθ1 = 0.098 rad, kT1 = 0.094 keV
#   spot 2: θc2 = 1.612 rad (92.36°), Δθ2 = 0.096 rad, kT2 = 0.094 keV (equal temp)
# Δφ2 = 0.558 cyc (anti-phased). Equatorial spots + equatorial observer ⇒ each spot's
# CENTER dips behind the star ~21% of the cycle. But anti-phasing puts those eclipses at
# opposite phases, so the combined pulse still never reaches zero (floor ≈ 0.63): PF stays
# unsaturated and ΔPF is large (+0.195) — the same regime as Riley, not J0030. A useful
# second-team confirmation that J0740's systematic is PF-visible regardless of the fit.
_MI_W1 = np.sin(0.098) ** 2 * 0.094 ** 4
_MI_W2 = np.sin(0.096) ** 2 * 0.094 ** 4
MILLER = Anchor(
    label="Miller 2021 (Illinois–Maryland)",
    inclination=1.527,
    compactness=0.444,
    spots=(
        Spot(colatitude=1.600, azimuth=0.000, weight=_MI_W1),
        Spot(colatitude=1.612, azimuth=0.558, weight=_MI_W2),
    ),
    note="weight ∝ area(sin²ζ) × T⁴; equatorial spots, centers eclipse at median",
)

ANCHORS = (RILEY, MILLER)


def sweep_anchor(anchor, tau_values, mu_centers, intensity_by_tau, bending=None):
    """ΔPF(τ) (the headline) plus the visibility numbers that explain its size.

    Both J0740 anchors keep PF unsaturated — Riley because no spot sets, Miller because
    the anti-phased spots tile the cycle even though each center dips behind ~21% of the
    time — so ΔPF is a live pulsed-fraction signal for both. The eclipse fraction and
    flux floor returned alongside it are what make that statement checkable.

    ``bending`` selects the light-bending map: ``None`` is Beloborodov's linear default
    (used by the ``a3_seed_errors`` linear cross-check against the v0.9.7 convergence
    study); the standalone :func:`run` passes an :class:`mcrt.bending.ExactBending` per
    anchor, because Gate G1 (v0.9.9.1) switched the J0740 anchors to exact Schwarzschild
    bending — the −0.0066 (2.1σ) linear bias at u = 0.494 exceeds the seed error bar.
    """
    iso = two_spot_flux(anchor.inclination, anchor.compactness, anchor.spots, None,
                        bending=bending)
    pf_iso = pulsed_fraction(iso)

    n_tau = len(tau_values)
    shape_idx = shape_tau_index(tau_values)  # raises if SHAPE_TAU left the grid
    pf_real = np.zeros(n_tau)
    delta_pf = np.zeros(n_tau)
    shape_rms = np.zeros(n_tau)
    real_at_shape_tau = None
    for ti in range(n_tau):
        beaming = beaming_lookup(mu_centers, intensity_by_tau[ti])
        real = two_spot_flux(anchor.inclination, anchor.compactness, anchor.spots,
                             beaming, bending=bending)
        pf_real[ti] = pulsed_fraction(real)
        delta_pf[ti] = pf_real[ti] - pf_iso
        shape_rms[ti], _ = waveform_shape_change(iso, real)
        if ti == shape_idx:
            real_at_shape_tau = real

    # Per-spot visibility: eclipse fraction (0 ⇒ never sets) and the grazing μ_min.
    eclipse = [single_spot_eclipsed_fraction(anchor.inclination, s.colatitude,
                                             anchor.compactness, bending=bending)
               for s in anchor.spots]
    min_mu = [single_spot_min_visible_mu(anchor.inclination, s.colatitude,
                                         anchor.compactness, bending=bending)
              for s in anchor.spots]

    return {
        "pf_iso": pf_iso,
        "pf_real": pf_real,
        "delta_pf": delta_pf,
        "shape_rms": shape_rms,
        "flux_floor": float(iso.min() / iso.max()),  # F_min/F_max; >0 ⇒ no eclipse
        "eclipse_fraction": eclipse,
        "min_visible_mu": min_mu,
        "delta_pf_bg": delta_pf_vs_background(iso, real_at_shape_tau),
        "iso_flux": iso,
        "real_flux": real_at_shape_tau,
    }


def plot_result(tau_values, results, path=FIGURE_PATH):
    """A: ΔPF(τ) — the live PF systematic (headline). B: the non-eclipsing pulse.
    C: ΔPF vs background — robustness of the headline."""
    fig, (ax_pf, ax_wave, ax_bg) = plt.subplots(1, 3, figsize=(16, 4.6))
    colors = {"Riley 2021 (X-PSI ST-U)": "#2c7fb8",
              "Miller 2021 (Illinois–Maryland)": "#c0392b"}

    # A — the headline: realistic beaming raises the pulsed fraction directly, because
    # J0740's pulse is unsaturated. Both teams' fits give a live, positive ΔPF (Riley
    # +0.137, Miller +0.195 at τ≈10) — the systematic is PF-visible regardless of the fit.
    for anchor in ANCHORS:
        r, c = results[anchor.label], colors[anchor.label]
        ax_pf.semilogx(tau_values, r["delta_pf"], "o-", color=c,
                       label=f"{anchor.label}  (PF_iso={r['pf_iso']:.2f})")
    ax_pf.axhline(0.0, color="#7f8c8d", lw=0.8)
    ax_pf.set_xlabel("Optical depth τ")
    ax_pf.set_ylabel(r"$\Delta\mathrm{PF} = \mathrm{PF}_{\rm real} - \mathrm{PF}_{\rm iso}$")
    ax_pf.set_title("Beaming systematic lands in the pulsed fraction\n"
                    "(J0740's anti-phased spots tile the cycle — PF unsaturated)")
    ax_pf.legend(fontsize=8)
    ax_pf.grid(True, alpha=0.3)

    # B — why PF has headroom: the Riley two-spot pulse is double-peaked and never near
    # zero (F_min/F_max ≈ 0.7); realistic beaming sharpens it (deeper dips, higher peaks).
    anchor = RILEY
    r = results[anchor.label]
    cyc = np.linspace(0.0, 1.0, N_PHASE, endpoint=False)
    ax_wave.plot(cyc, r["iso_flux"] / r["iso_flux"].max(), "-", color="#7f8c8d", lw=2,
                 label="two-spot isotropic")
    ax_wave.plot(cyc, r["real_flux"] / r["real_flux"].max(), "-", color="#2c7fb8", lw=2,
                 label=fr"two-spot realistic $I(\mu;\tau={SHAPE_TAU:g})$")
    ax_wave.axhline(0.0, color="#bdc3c7", lw=0.8)
    ax_wave.set_xlabel("Rotational phase (cycles)")
    ax_wave.set_ylabel(r"Normalized flux $F(\phi)/F_{\max}$")
    ax_wave.set_title(f"Non-eclipsing two-spot pulse — {anchor.label}\n"
                      f"(floor $F_{{\\min}}/F_{{\\max}}={r['flux_floor']:.2f}>0$)")
    ax_wave.legend(fontsize=8)
    ax_wave.grid(True, alpha=0.3)

    # C — robustness: a real unpulsed background only dilutes the (positive) Riley ΔPF;
    # it does not depend on a knife-edge zero the way J0030's caveat panel did.
    for anchor in ANCHORS:
        r, c = results[anchor.label], colors[anchor.label]
        ax_bg.plot(BACKGROUND_FRACS, r["delta_pf_bg"], "o-", color=c, label=anchor.label)
    ax_bg.axhline(0.0, color="#7f8c8d", lw=0.8)
    ax_bg.set_xlabel("Assumed unpulsed background (fraction of peak)")
    ax_bg.set_ylabel(r"$\Delta\mathrm{PF}$")
    ax_bg.set_title(f"Robustness: ΔPF vs background (τ={SHAPE_TAU:g})\n"
                    "(both teams stay positive; not a saturation artifact)")
    ax_bg.legend(fontsize=8)
    ax_bg.grid(True, alpha=0.3)

    fig.tight_layout()
    fig.savefig(path, dpi=130)
    print(f"✓ figure saved to {path}")


def save_results(tau_values, results, path=RESULTS_PATH):
    payload = {"tau_values": tau_values, "background_fracs": BACKGROUND_FRACS}
    for anchor in ANCHORS:
        key = "riley" if "Riley" in anchor.label else "miller"
        r = results[anchor.label]
        payload[f"{key}_pf_iso"] = r["pf_iso"]
        payload[f"{key}_pf_real"] = r["pf_real"]
        payload[f"{key}_delta_pf"] = r["delta_pf"]
        payload[f"{key}_shape_rms"] = r["shape_rms"]
        payload[f"{key}_flux_floor"] = r["flux_floor"]
        payload[f"{key}_eclipse_fraction"] = np.array(r["eclipse_fraction"])
        payload[f"{key}_min_visible_mu"] = np.array(r["min_visible_mu"])
        payload[f"{key}_delta_pf_bg"] = r["delta_pf_bg"]
    np.savez(path, **payload)
    print(f"✓ results saved to {path}")


def print_summary(tau_values, results):
    print("\n" + "=" * 74)
    print("Rung D′ — beaming swap at PSR J0740+6620 (the non-eclipsing complement)")
    print("=" * 74)
    for anchor in ANCHORS:
        r = results[anchor.label]
        peak = int(np.argmax(r["delta_pf"]))
        ecl = ", ".join(f"{e:.0%}" for e in r["eclipse_fraction"])
        mu = ", ".join(f"{m:.3f}" for m in r["min_visible_mu"])
        print(f"\n{anchor.label}   (u = {anchor.compactness:.3f}, "
              f"i = {np.rad2deg(anchor.inclination):.1f}°)")
        print(f"  weighting: {anchor.note}")
        print(f"  per-spot eclipse fraction:      {ecl}   (0% ⇒ never sets)")
        print(f"  per-spot grazing μ_min:         {mu}")
        print(f"  two-spot flux floor F_min/F_max: {r['flux_floor']:.3f}  "
              f"({'combined pulse off zero — PF unsaturated' if r['flux_floor'] > 1e-3 else 'combined pulse hits zero — PF saturates'})")
        print(f"  isotropic pulsed fraction:      PF = {r['pf_iso']:.3f}")
        print(f"  ΔPF peak:                       {r['delta_pf'][peak]:+.3f} at τ≈{tau_values[peak]:g}  "
              f"(PF_real = {r['pf_real'][peak]:.3f})")
        print(f"  ΔPF with background "
              f"({BACKGROUND_FRACS[1]:.0%}–{BACKGROUND_FRACS[-1]:.0%}): "
              f"{r['delta_pf_bg'][1]:+.3f} … {r['delta_pf_bg'][-1]:+.3f}")
    print("\nReading: J0740's two hot spots are anti-phased, so they tile the rotation and\n"
          "the combined pulse never reaches zero (floor 0.6–0.7). PF is therefore low and\n"
          "unsaturated, and the beaming systematic lands directly in the pulsed fraction —\n"
          "ΔPF = +0.137 (Riley) to +0.195 (Miller), the same effect v0.9.0 found on invented\n"
          "geometries, now on a real star and confirmed by both teams. (Miller's spots even\n"
          "dip behind 21% each, but anti-phasing keeps the pulse off zero.) Contrast J0030:\n"
          "same-hemisphere spots, deep 45% eclipses that don't tile → PF saturates, ΔPF ≈ 0,\n"
          "systematic in shape. Tiling — not single-spot eclipse — decides where it shows up.")


def run():
    """Standalone run — **exact Schwarzschild bending** (Gate G1, v0.9.9.1).

    Each anchor gets its own :class:`mcrt.bending.ExactBending(u)`; the linear
    map is retained only in ``a3_seed_errors`` as the convergence cross-check.
    """
    tau_values, mu_centers, intensity_by_tau = load_library()
    results = {}
    for a in ANCHORS:
        bending = ExactBending(a.compactness)
        results[a.label] = sweep_anchor(a, tau_values, mu_centers,
                                        intensity_by_tau, bending=bending)
    return tau_values, results


if __name__ == "__main__":
    tau_values, results = run()
    save_results(tau_values, results)
    plot_result(tau_values, results)
    print_summary(tau_values, results)
