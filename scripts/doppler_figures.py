"""Figures for the v0.9.10 Doppler deep-dive (Track C, presentation only).

Three deterministic graphics, all recomputed from the production library + the
committed Track-C code paths (no numbers are stored here — if a driver changes,
rerunning this script regenerates consistent figures):

1. ``doppler_dpf_vs_spin.png`` — ΔPF(ν) for both J0740 anchors (δ⁴, exact
   bending, pooled curve): the smooth, monotonic collapse of the frozen headline
   with spin, with the real 346.5 Hz marked. The §5 ν-table as a curve.
2. ``doppler_band_amplification.png`` — why the NICER band amplifies Doppler:
   (a) the observed blackbody photon spectrum with the band shaded — the hump
   sits *below* the band, so NICER watches the Wien tail; a ±β temperature
   wobble swings the in-band counts hugely; (b) the flux weight vs δ — the
   in-band photon weight is steeper than both bolometric conventions
   (effective exponent ≈ 6 vs 3/4) across the actual δ range at J0740.
3. ``doppler_dilution_ladder.png`` — ΔPF under the full treatment ladder
   (frozen → δ³ → δ⁴ → band energy → band photon → band photon + time-delay
   warp) for both anchors, with per-seed error bars: the C2 → C3 → C4 story in
   one picture, ending at the quotable all-geometry residuals.

Run from the repository root:  PYTHONPATH=src:scripts python3 scripts/doppler_figures.py
"""

import matplotlib.pyplot as plt
import numpy as np

from mcrt import ExactBending, Rotation, band_boost, beaming_lookup, pulsed_fraction, \
    spot_speed
from anchor_lib import multi_spot_flux, shape_tau_index
from c2_doppler_coupling import J0740_SPIN_HZ, RADIUS_KM
from c3_band_doppler import BAND_KEV, KT_KEV, RILEY, anchor_band
from c4_caveat_audit import delayed_multi_spot_flux, per_seed_dpf
import j0740_anchor as j0740

LIBRARY_PATH = "data/beaming_library.npz"
SPIN_FIG = "data/doppler_dpf_vs_spin.png"
BAND_FIG = "data/doppler_band_amplification.png"
LADDER_FIG = "data/doppler_dilution_ladder.png"

# Fixed two-hue assignment used across the repo's figures (CVD-distinct pair).
COLOR = {"riley": "#1f77b4", "miller": "#d62728"}
GRAY = "#7f7f7f"


def anchor_key(anchor) -> str:
    return "riley" if "Riley" in anchor.label else "miller"


def dpf_vs_spin_figure(pooled_curve, mu_centers):
    """ΔPF(ν) curves, δ⁴ / exact bending / pooled library curve."""
    beaming = beaming_lookup(mu_centers, pooled_curve)
    spins = np.linspace(0.0, 600.0, 41)
    fig, ax = plt.subplots(figsize=(8.0, 4.6), layout="constrained")
    for anchor in j0740.ANCHORS:
        bend = ExactBending(anchor.compactness)
        r_km = RADIUS_KM[anchor.label]
        key = anchor_key(anchor)
        dpf = []
        for nu in spins:
            rot = None if nu == 0.0 else Rotation(nu, r_km)
            iso = multi_spot_flux(anchor.inclination, anchor.compactness, anchor.spots,
                                  None, bending=bend, rotation=rot)
            real = multi_spot_flux(anchor.inclination, anchor.compactness, anchor.spots,
                                   beaming, bending=bend, rotation=rot)
            dpf.append(pulsed_fraction(real) - pulsed_fraction(iso))
        ax.plot(spins, dpf, color=COLOR[key], lw=2.0)
        ax.annotate(anchor.label.split(" (")[0], xy=(spins[-1], dpf[-1]),
                    xytext=(-4, 8), textcoords="offset points", ha="right",
                    color=COLOR[key], fontsize=9)
        nu_real = float(np.interp(J0740_SPIN_HZ, spins, dpf))
        ax.plot([J0740_SPIN_HZ], [nu_real], "o", ms=7, color=COLOR[key],
                mec="white", mew=1.5)
    ax.axvline(J0740_SPIN_HZ, color=GRAY, lw=0.9, ls=":")
    ax.annotate("real spin 346.5 Hz", xy=(J0740_SPIN_HZ, 0.185), xytext=(6, 0),
                textcoords="offset points", color=GRAY, fontsize=9)
    ax.set_xlabel("spin frequency ν  [Hz]")
    ax.set_ylabel("ΔPF = PF_real − PF_iso   (δ⁴, exact bending)")
    ax.set_title("The Doppler boost drains the pulsed-fraction systematic smoothly with spin")
    ax.grid(alpha=0.25, lw=0.6)
    ax.set_ylim(bottom=0.0)
    ax.margins(x=0.02)
    fig.savefig(SPIN_FIG, dpi=130)
    print(f"saved {SPIN_FIG}")


def band_amplification_figure():
    """(a) spectrum + shaded band; (b) flux weight vs δ — the Wien-tail amplifier."""
    u = next(a.compactness for a in j0740.ANCHORS if "Riley" in a.label)
    kt, (e1, e2) = KT_KEV[RILEY], BAND_KEV[RILEY]
    beta = spot_speed(J0740_SPIN_HZ, RADIUS_KM[RILEY], np.pi / 2.0, u)
    gamma = 1.0 / np.sqrt(1.0 - beta**2)
    delta_hi, delta_lo = 1.0 / (gamma * (1.0 - beta)), 1.0 / (gamma * (1.0 + beta))
    fig, (ax_s, ax_w) = plt.subplots(1, 2, figsize=(10.0, 4.4), layout="constrained")

    energy = np.linspace(0.01, 1.6, 400)
    t0 = np.sqrt(1.0 - u) * kt
    static = energy**2 / np.expm1(energy / t0)
    norm = static.max()
    for delta, label, color in ((delta_hi, f"approaching (δ = {delta_hi:.2f})",
                                 COLOR["riley"]),
                                (1.0, "static (δ = 1)", "#333333"),
                                (delta_lo, f"receding (δ = {delta_lo:.2f})",
                                 COLOR["miller"])):
        spec = energy**2 / np.expm1(energy / (delta * t0))
        ax_s.plot(energy, spec / norm, color=color, lw=1.8, label=label)
    ax_s.axvspan(e1, e2, color=GRAY, alpha=0.15, lw=0)
    ax_s.annotate("NICER band", xy=(0.5 * (e1 + e2), 0.5), ha="center",
                  color="#555555", fontsize=9)
    ax_s.set_yscale("log")
    ax_s.set_ylim(1e-9, 3.0)
    ax_s.set_xlabel("observed photon energy E  [keV]")
    ax_s.set_ylabel(r"photon spectrum  $E^2/(e^{E/T_\mathrm{obs}} - 1)$  (normalized)")
    ax_s.set_title(f"The band sits on the Wien tail (kT = {kt} keV, u = {u:.3f})",
                   fontsize=10)
    ax_s.legend(fontsize=8, loc="lower left")
    ax_s.grid(alpha=0.25, lw=0.6)

    deltas = np.linspace(delta_lo - 0.01, delta_hi + 0.01, 200)
    ax_w.plot(deltas, deltas**3, color=GRAY, lw=1.6, ls=":", label="δ³ (bolometric photon)")
    ax_w.plot(deltas, deltas**4, color=GRAY, lw=1.6, ls="--", label="δ⁴ (bolometric energy)")
    ax_w.plot(deltas, band_boost(deltas, u, anchor_band(RILEY)), color=COLOR["riley"],
              lw=2.2, label="in-band photon weight (n_eff ≈ 6)")
    ax_w.axvspan(delta_lo, delta_hi, color=GRAY, alpha=0.12, lw=0)
    ax_w.annotate("δ range at β = 0.127", xy=(1.0, 1.72), ha="center",
                  color="#555555", fontsize=9)
    ax_w.axhline(1.0, color="#bbbbbb", lw=0.8)
    ax_w.axvline(1.0, color="#bbbbbb", lw=0.8)
    ax_w.set_xlabel("Doppler factor δ")
    ax_w.set_ylabel("flux weight  (normalized to δ = 1)")
    ax_w.set_title("In-band response is steeper than any bolometric δⁿ", fontsize=10)
    ax_w.legend(fontsize=8, loc="upper left")
    ax_w.grid(alpha=0.25, lw=0.6)

    fig.suptitle("Why the NICER band amplifies the Doppler dilution (Riley J0740 parameters)",
                 fontsize=11)
    fig.savefig(BAND_FIG, dpi=130)
    print(f"saved {BAND_FIG}")


def dilution_ladder_figure(seed_curves, mu_centers):
    """ΔPF under the C2→C4 treatment ladder, per-seed error bars, both anchors."""
    treatments = [
        ("frozen (0 Hz)", lambda a, r: None),
        ("δ³ bolometric photon", lambda a, r: Rotation(J0740_SPIN_HZ, r, photon_flux=True)),
        ("δ⁴ bolometric energy", lambda a, r: Rotation(J0740_SPIN_HZ, r)),
        ("band-limited energy", lambda a, r: Rotation(J0740_SPIN_HZ, r,
                                                      band=anchor_band(a.label))),
        ("band-limited photon", lambda a, r: Rotation(J0740_SPIN_HZ, r, photon_flux=True,
                                                      band=anchor_band(a.label))),
        ("band photon + time-delay warp", "warped"),
    ]
    fig, ax = plt.subplots(figsize=(8.0, 4.8), layout="constrained")
    y_pos = np.arange(len(treatments))[::-1]
    for offset, anchor in zip((+0.13, -0.13), j0740.ANCHORS):
        bend = ExactBending(anchor.compactness)
        r_km = RADIUS_KM[anchor.label]
        key = anchor_key(anchor)
        means, sigmas = [], []
        for _, rot_factory in treatments:
            if rot_factory == "warped":
                rot = Rotation(J0740_SPIN_HZ, r_km, photon_flux=True,
                               band=anchor_band(anchor.label))

                def flux_fn(a, b, _bend=bend, _rot=rot, _r=r_km):
                    return delayed_multi_spot_flux(a, b, _bend, _rot, _r, J0740_SPIN_HZ)
            else:
                rot = rot_factory(anchor, r_km)

                def flux_fn(a, b, _bend=bend, _rot=rot):
                    return multi_spot_flux(a.inclination, a.compactness, a.spots, b,
                                           bending=_bend, rotation=_rot)
            dpf = per_seed_dpf(anchor, mu_centers, seed_curves, flux_fn)
            means.append(dpf.mean())
            sigmas.append(dpf.std(ddof=1))
        ax.errorbar(means, y_pos + offset, xerr=sigmas, fmt="o", ms=6.5,
                    color=COLOR[key], ecolor=COLOR[key], elinewidth=1.4, capsize=3,
                    label=anchor.label.split(" (")[0])
        for x, y in zip(means, y_pos + offset):
            ax.annotate(f"{x:+.3f}", xy=(x, y), xytext=(6, -3),
                        textcoords="offset points", fontsize=8, color="#444444")
    ax.set_yticks(y_pos)
    ax.set_yticklabels([name for name, _ in treatments], fontsize=9)
    ax.axvline(0.0, color="#bbbbbb", lw=0.8)
    ax.set_xlabel("ΔPF = PF_real − PF_iso  (τ = 10, exact bending; error bars = seed σ)")
    ax.set_title("The dilution ladder: each physical refinement, C2 → C3 → C4")
    ax.legend(fontsize=9, loc="lower right")
    ax.grid(axis="x", alpha=0.25, lw=0.6)
    ax.margins(x=0.08)
    fig.savefig(LADDER_FIG, dpi=130)
    print(f"saved {LADDER_FIG}")


def main():
    d = np.load(LIBRARY_PATH)
    mu = d["mu_centers"]
    ti10 = shape_tau_index(d["tau_values"])
    dpf_vs_spin_figure(d["intensity_by_tau"][ti10], mu)
    band_amplification_figure()
    dilution_ladder_figure(d["intensity_by_tau_seed"][ti10], mu)


if __name__ == "__main__":
    main()
