"""Real-star anchor: the isotropic-vs-realistic beaming swap at PSR J0030+0451.

Earlier work measured the beaming systematic ΔPF = PF_real − PF_iso on three
*invented* always-visible geometries (``beaming_pulse_sweep.py``), where it reached
the ≈ +0.15 scale. This script anchors the same swap at PSR J0030+0451's published geometry —
the canonical NICER mass–radius pulsar — using two independent fits:

    Riley et al. (2019, ApJL 887 L21)   — X-PSI / Amsterdam, ST+PST model
    Miller et al. (2019, ApJL 887 L24)  — Illinois–Maryland, three-oval-spot model

It is a *differential* statement (``to_finish.md`` §6): a relative change at a
realistic operating point, not a fit to J0030's data.

The result — and why it is interesting. Both fits place J0030's hot spots in the
*same* hemisphere, tilted away from us (colatitudes 128–167°), viewed at i ≈ 50–54°.
So the spots dive behind the star: a point spot is eclipsed for ~45% of the rotation,
and even summing the two published spots the flux still touches zero (they hug the
same pole and do not tile the cycle). A zero minimum pins the pulsed fraction at
PF = 1 for *both* beamings, so the **pulsed-fraction** systematic saturates to ΔPF ≈ 0.

This is not a null finding — it is the same physics Rung C already showed: ΔPF is
large at intermediate geometries and vanishes at near-saturated ones (the Rung C
"high-contrast" case gave only +0.018). J0030 sits squarely in that saturated corner.
The systematic does not disappear; it moves out of the pulsed fraction and into the
**waveform shape** — assuming isotropy distorts the visible pulse by several percent
RMS. So the systematic's *observability* is geometry-dependent: PF-visible for
intermediate geometries, shape-visible at J0030.

Idealization note. Real spots are extended and sit on an unpulsed background, so the
true minimum is non-zero rather than a hard eclipse. We do not model a background
(it would be a tuned knob); panel C shows that *were* one present, a PF systematic
would re-emerge — and with the *opposite* sign to Rung C, because J0030's spots only
ever reach modest emission angles (μ ≲ 0.45), so limb darkening dims their peak
instead of sharpening the contrast.

Built here in the script with no change to the verified ``mcrt.pulse`` core: each
spot is a :func:`mcrt.pulse.compute_profile` call, a spot at a different longitude is
that profile rolled in phase (``np.roll``; separations land on the 1024-point grid),
and light is additive so the star's flux is the weighted sum.

Run from the repository root:  python scripts/j0030_anchor.py
"""

import numpy as np
import matplotlib.pyplot as plt

from mcrt import beaming_lookup, compute_profile, pulsed_fraction
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
    two_spot_flux,
    waveform_shape_change,
)

FIGURE_PATH = "data/pulse_profile_j0030.png"
RESULTS_PATH = "data/j0030_anchor.npz"


# --- Published geometries (verified against the papers' tables) ----------------
#
# Riley 2019 ST+PST, Table 2: M = 1.34 M⊙, R_eq = 12.71 km; compactness GM/Rc² =
# 0.156 ⇒ u = 0.312. i = 0.94 rad. Spot 1 (ST, small circular): Θp = 2.23 rad,
# log₁₀T = 6.11. Spot 2 (PST crescent): Θs = 2.91 rad, log₁₀T = 6.11 (equal temp).
# Azimuthal separation ≈ 0.45 cyc (φp = 0.46 vs φs = −0.59 "from Earth antipode").
# The crescent is not cleanly point-reducible to an area, so Riley uses EQUAL weights;
# the equal-vs-physical check (Miller) shows the swap barely depends on weighting.
RILEY = Anchor(
    label="Riley 2019 (X-PSI ST+PST)",
    inclination=0.94,
    compactness=0.312,
    spots=(
        Spot(colatitude=2.23, azimuth=0.00, weight=1.0),  # ST small circular
        Spot(colatitude=2.91, azimuth=0.45, weight=1.0),  # PST crescent
    ),
    note="equal weight (crescent area not cleanly point-reducible; temps equal)",
)

# Miller 2019 three-oval-spot, Table 8 (medians): M = 1.443 M⊙, R_e = 13.019 km;
# compactness GM/Rc² = 0.163 ⇒ u = 0.326. θobs = 0.878 rad. The two MAIN spots (the
# third is tiny — "very small contribution" — and is dropped):
#   Spot 1: θc1 = 2.270 rad, Δθ1 = 0.036, oval ratio f1 = 5.352, kT1 = 0.117 keV
#   Spot 2: θc2 = 2.417 rad, Δθ2 = 0.033, oval ratio f2 = 15.769, kT2 = 0.115 keV,
#           Δφ2 = 0.460 cyc (explicit longitude separation)
# Oval area ∝ Δθ²·f, bolometric weight ∝ area × T⁴ — a clean physical weighting.
_M_W1 = (0.036**2 * 5.352) * 0.117**4
_M_W2 = (0.033**2 * 15.769) * 0.115**4
MILLER = Anchor(
    label="Miller 2019 (Illinois–Maryland)",
    inclination=0.878,
    compactness=0.326,
    spots=(
        Spot(colatitude=2.270, azimuth=0.000, weight=_M_W1),
        Spot(colatitude=2.417, azimuth=0.460, weight=_M_W2),
    ),
    note="weight ∝ area(Δθ²·f) × T⁴ from Table 8",
)

ANCHORS = (RILEY, MILLER)


def sweep_anchor(anchor, tau_values, mu_centers, intensity_by_tau):
    """Per-τ waveform-shape change (the headline) plus the saturation/caveat numbers."""
    iso = two_spot_flux(anchor.inclination, anchor.compactness, anchor.spots, None)
    shape_idx = shape_tau_index(tau_values)  # raises if SHAPE_TAU left the grid
    shape_rms = np.zeros(len(tau_values))
    shape_max = np.zeros(len(tau_values))
    real_at_shape_tau = None
    for ti in range(len(tau_values)):
        beaming = beaming_lookup(mu_centers, intensity_by_tau[ti])
        real = two_spot_flux(anchor.inclination, anchor.compactness, anchor.spots, beaming)
        shape_rms[ti], shape_max[ti] = waveform_shape_change(iso, real)
        if ti == shape_idx:
            real_at_shape_tau = real

    return {
        "pf_iso": pulsed_fraction(iso),
        "eclipse_fraction": single_spot_eclipsed_fraction(
            anchor.inclination, anchor.spots[0].colatitude, anchor.compactness),
        "shape_rms": shape_rms,
        "shape_max": shape_max,
        "delta_pf_bg": delta_pf_vs_background(iso, real_at_shape_tau),
        "iso_flux": iso,
        "real_flux": real_at_shape_tau,
    }


def plot_result(tau_values, results, path=FIGURE_PATH):
    """A: shape change vs τ (headline). B: the eclipsing two-spot pulse. C: ΔPF caveat."""
    fig, (ax_shape, ax_wave, ax_bg) = plt.subplots(1, 3, figsize=(16, 4.6))
    colors = {"Riley 2019 (X-PSI ST+PST)": "#2c7fb8",
              "Miller 2019 (Illinois–Maryland)": "#c0392b"}

    # A — the systematic that survives the eclipse: waveform-shape change, rising with τ.
    for anchor in ANCHORS:
        r, c = results[anchor.label], colors[anchor.label]
        ax_shape.semilogx(tau_values, r["shape_rms"], "o-", color=c, label=f"{anchor.label}: RMS")
        ax_shape.semilogx(tau_values, r["shape_max"], "s--", color=c, alpha=0.5,
                          markersize=4, label=f"{anchor.label}: max-local")
    ax_shape.set_xlabel("Optical depth τ")
    ax_shape.set_ylabel("Waveform shape change (normalized flux)")
    ax_shape.set_title("Beaming reshapes the J0030 pulse\n(PF is saturated; the shape is not)")
    ax_shape.legend(fontsize=7)
    ax_shape.grid(True, alpha=0.3)

    # B — why PF is blind: the two-spot pulse eclipses to zero; a single spot too.
    anchor = MILLER
    r = results[anchor.label]
    cyc = np.linspace(0.0, 1.0, N_PHASE, endpoint=False)
    single = compute_profile(anchor.inclination, anchor.spots[0].colatitude,
                             anchor.compactness).flux
    ax_wave.plot(cyc, single / single.max(), ":", color="#95a5a6", lw=1.5,
                 label="one spot only")
    ax_wave.plot(cyc, r["iso_flux"] / r["iso_flux"].max(), "-", color="#7f8c8d", lw=2,
                 label="two-spot isotropic")
    ax_wave.plot(cyc, r["real_flux"] / r["real_flux"].max(), "-", color="#c0392b", lw=2,
                 label=fr"two-spot realistic $I(\mu;\tau={SHAPE_TAU:g})$")
    ax_wave.set_xlabel("Rotational phase (cycles)")
    ax_wave.set_ylabel(r"Normalized flux $F(\phi)/F_{\max}$")
    ax_wave.set_title(f"Eclipsing two-spot pulse — {anchor.label}")
    ax_wave.legend(fontsize=8)
    ax_wave.grid(True, alpha=0.3)

    # C — caveat: a background would un-saturate PF and expose a sign-flipped systematic.
    for anchor in ANCHORS:
        r, c = results[anchor.label], colors[anchor.label]
        ax_bg.plot(BACKGROUND_FRACS, r["delta_pf_bg"], "o-", color=c, label=anchor.label)
    ax_bg.axhline(0.0, color="#7f8c8d", lw=0.8)
    ax_bg.set_xlabel("Assumed unpulsed background (fraction of peak)")
    ax_bg.set_ylabel(r"$\Delta\mathrm{PF}$ if de-saturated")
    ax_bg.set_title(f"Caveat: PF systematic vs. background\n(τ={SHAPE_TAU:g}; "
                    "negative — opposite to Rung C)")
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
        payload[f"{key}_eclipse_fraction"] = r["eclipse_fraction"]
        payload[f"{key}_shape_rms"] = r["shape_rms"]
        payload[f"{key}_shape_max"] = r["shape_max"]
        payload[f"{key}_delta_pf_bg"] = r["delta_pf_bg"]
    np.savez(path, **payload)
    print(f"✓ results saved to {path}")


def print_summary(tau_values, results):
    print("\n" + "=" * 74)
    print("Rung D — beaming swap at PSR J0030+0451's published geometry")
    print("=" * 74)
    for anchor in ANCHORS:
        r = results[anchor.label]
        peak = int(np.argmax(r["shape_rms"]))
        print(f"\n{anchor.label}   (u = {anchor.compactness:.3f}, "
              f"i = {np.rad2deg(anchor.inclination):.1f}°)")
        print(f"  weighting: {anchor.note}")
        print(f"  single-spot eclipse fraction:   {r['eclipse_fraction']:.0%} of the cycle")
        print(f"  two-spot pulsed fraction:       PF = {r['pf_iso']:.3f}  "
              f"(saturated → ΔPF ≈ 0 in PF)")
        print(f"  waveform shape change at τ≈{tau_values[peak]:g}:  "
              f"RMS {r['shape_rms'][peak]:.3f}, max-local {r['shape_max'][peak]:.3f}")
        print(f"  caveat — ΔPF if de-saturated by background "
              f"({BACKGROUND_FRACS[1]:.0%}–{BACKGROUND_FRACS[-1]:.0%}): "
              f"{r['delta_pf_bg'][1]:+.3f} … {r['delta_pf_bg'][-1]:+.3f}")
    print("\nReading: J0030's same-hemisphere spots eclipse, so the pulsed fraction is "
          "saturated\nand blind to the beaming systematic — exactly the near-saturated "
          "corner where Rung C\nalso found ΔPF→0. The systematic instead reshapes the "
          "visible waveform (a few % RMS),\nso its observability is geometry-dependent: "
          "PF-visible mid-geometry, shape-visible here.")


def run():
    tau_values, mu_centers, intensity_by_tau = load_library()
    results = {a.label: sweep_anchor(a, tau_values, mu_centers, intensity_by_tau)
               for a in ANCHORS}
    return tau_values, results


if __name__ == "__main__":
    tau_values, results = run()
    save_results(tau_values, results)
    plot_result(tau_values, results)
    print_summary(tau_values, results)
