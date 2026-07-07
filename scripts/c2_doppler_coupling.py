"""Track C2 — the measured Doppler coupling to the headline ΔPF (Gate G2).

C1 (v0.9.10) built the Doppler + aberration layer ``mcrt.rotating`` and validated it
against the 200 Hz SD1c code-comparison waveform (residual ~1%, the linear-bending
floor). C2 answers the question the track exists for: **does the rotational velocity
of the real stars change our pulsed-fraction systematic ΔPF**, and by how much?

The stars are not slow — J0740 spins at 346.5 Hz (β ≈ 0.12 c at the caps), J0030 at
205.5 Hz. Our headline was computed in the frozen slow-rotation limit. The naïve
expectation was that the Doppler *boost* (δ⁴) cancels in ΔPF = PF_real − PF_iso
(it hits both terms identically) leaving only **aberration** (μ' = δ cos α) to
survive. **The measurement shows that expectation is wrong**: PF is a nonlinear
(max/min) statistic, and the boost — strongest exactly at the limb phases where
limb-darkening acts — dominates the coupling (deep-dive §6). C2 measures the
*difference of differences*,

    coupling = ΔPF(ν) − ΔPF(0),   ΔPF(ν) = PF_real(ν) − PF_iso(ν),

with iso and real sharing the *same* Doppler treatment, on the production 5-seed
library and the exact Schwarzschild bending the anchors already use (Gate G1). σ is
the seed-to-seed Monte-Carlo error. The coupling is a **systematic** (a fixed
modeling bias; more photons shrink σ, not the coupling), so Gate G2 keys on
|coupling| vs σ:

    |coupling| ≤ σ  → one quantified caveat sentence ("Doppler shifts ΔPF by < X").
    |coupling| >  σ → report ΔPF(ν) as a quantified correction alongside the headline.

The δ⁴ **bolometric-energy** boost is used (the production I(μ) is a grey energy
intensity); the δ³ photon-flux boost is a documented alternative that changes the
boost — not the aberration coupling — by O(β).

Deliberately excluded (geometry-only, cancel differentially in ΔPF — quantified
caveats, not modelled): photon light-travel-time delay (~2% of phase) and rotational
oblateness (~3% at 346 Hz, AlGendy & Morsink 2014).

**Shape routing (the companion measurement).** The ΔPF collapse does *not* mean the
beaming signal vanishes: the iso-vs-real *waveform shape* difference (peak-normalized
RMS) and the fundamental-harmonic difference Δ(A1/A0) are nearly spin-invariant. What
spin does is pump a large common second harmonic (the δ-boost is a once-per-cycle
modulation multiplying a once-per-cycle pulse) into both models, saturating the
min/max PF statistic as a probe of beaming while the signal survives in shape —
the same routing logic as the J0030 eclipse story, with spin as a second router.
:func:`shape_routing` measures exactly that, and the four-pulse figure shows it.

Pure interpolation + geometry (no Monte Carlo), seconds to run. The ν = 0 branch
reproduces the exact-bending headline (A3/B2) bit-for-bit; only the rotating branch
is new.

Run from the repository root:  PYTHONPATH=src:scripts python3 scripts/c2_doppler_coupling.py
"""

import matplotlib.pyplot as plt
import numpy as np

from mcrt import ExactBending, Rotation, beaming_lookup, pulsed_fraction, spot_speed
from anchor_lib import SHAPE_TAU, multi_spot_flux, shape_tau_index, waveform_shape_change
import j0030_anchor as j0030
import j0740_anchor as j0740

LIBRARY_PATH = "data/beaming_library.npz"
RESULTS_PATH = "data/c2_doppler_coupling.npz"
FIGURE_PATH = "data/pulse_profile_doppler_routing.png"

# Measured spin frequencies (as seen at infinity) and published radii per fit.
J0740_SPIN_HZ = 346.53
J0030_SPIN_HZ = 205.53
RADIUS_KM = {
    "Riley 2021 (X-PSI ST-U)": 12.39,        # Riley 2021 J0740, R_eq
    "Miller 2021 (Illinois–Maryland)": 13.713,  # Miller 2021 J0740, R_e
    "Riley 2019 (X-PSI ST+PST)": 12.71,      # Riley 2019 J0030, R_eq
    "Miller 2019 (Illinois–Maryland)": 13.019,  # Miller 2019 J0030, R_e
}


def per_seed_delta_pf(anchor, mu_centers, seed_curves_at_tau, bending, rotation):
    """(mean, std, pf_iso) of ΔPF = PF_real − PF_iso across seeds, one (bend, rot) pair.

    Both PF_iso and each PF_real use the *same* ``bending`` and ``rotation``, so the
    returned spread is the seed MC error and the ν→ν' difference is the pure Doppler
    coupling (nothing else changes between the two calls).
    """
    iso = multi_spot_flux(anchor.inclination, anchor.compactness, anchor.spots,
                          None, bending=bending, rotation=rotation)
    pf_iso = pulsed_fraction(iso)

    n_seeds = seed_curves_at_tau.shape[0]
    dpf = np.empty(n_seeds)
    for si in range(n_seeds):
        beaming = beaming_lookup(mu_centers, seed_curves_at_tau[si])
        real = multi_spot_flux(anchor.inclination, anchor.compactness, anchor.spots,
                               beaming, bending=bending, rotation=rotation)
        dpf[si] = pulsed_fraction(real) - pf_iso
    return dpf.mean(), dpf.std(ddof=1), pf_iso


def j0740_coupling(seed_curves, mu_centers, ti10):
    """Per-anchor ΔPF(0) vs ΔPF(346.5 Hz) ± σ and the |coupling|/σ ratio."""
    rows = []
    curves = seed_curves[ti10]  # (n_seeds, n_bins)
    print("=== J0740 (non-eclipsing, PF live): ΔPF(τ=10) = PF_real − PF_iso ===")
    print(f"{'anchor':<34}{'β_eq':>7}{'ΔPF(0 Hz)':>16}{'ΔPF(346.5 Hz)':>18}"
          f"{'coupling':>11}{'|Δ|/σ':>8}")
    for anchor in j0740.ANCHORS:
        bend = ExactBending(anchor.compactness)
        rot = Rotation(spin_hz=J0740_SPIN_HZ, radius_km=RADIUS_KM[anchor.label])
        beta_eq = spot_speed(J0740_SPIN_HZ, RADIUS_KM[anchor.label],
                             np.pi / 2.0, anchor.compactness)

        s0_m, s0_s, _ = per_seed_delta_pf(anchor, mu_centers, curves, bend, None)
        sr_m, sr_s, _ = per_seed_delta_pf(anchor, mu_centers, curves, bend, rot)
        coupling = sr_m - s0_m
        sigma = max(s0_s, sr_s)
        ratio = abs(coupling) / sigma if sigma > 0 else np.inf
        key = "riley" if "Riley" in anchor.label else "miller"
        print(f"{anchor.label:<34}{beta_eq:>7.3f}"
              + f"{s0_m:+.4f}±{s0_s:.4f}".rjust(16)
              + f"{sr_m:+.4f}±{sr_s:.4f}".rjust(18)
              + f"{coupling:+.4f}".rjust(11)
              + f"{ratio:.1f}σ".rjust(8))
        rows.append({"key": key, "u": anchor.compactness, "beta_eq": beta_eq,
                     "dpf0_m": s0_m, "dpf0_s": s0_s, "dpfr_m": sr_m, "dpfr_s": sr_s,
                     "coupling": coupling, "sigma": sigma, "ratio": ratio})
    return rows


def j0030_unmoved(pooled, mu_centers, ti10):
    """Confirm the eclipsing anchor stays saturated (ΔPF ≈ 0) under rotation.

    Rotation cannot un-eclipse a spot — visibility (cos α ≥ 0) is a static-frame
    geometric condition, unchanged by the velocity — so PF stays pinned at 1 and the
    systematic remains shape-only. Uses the pooled library curve (qualitative point).
    """
    print("\n=== J0030 (eclipsing, PF saturated): rotation must not move ΔPF ===")
    print(f"{'anchor':<34}{'β_eq':>7}{'ΔPF(0 Hz)':>12}{'ΔPF(205.5 Hz)':>15}"
          f"{'floor(0)':>10}{'floor(ν)':>10}")
    beaming = beaming_lookup(mu_centers, pooled[ti10])
    rows = []
    for anchor in j0030.ANCHORS:
        bend = ExactBending(anchor.compactness)
        rot = Rotation(spin_hz=J0030_SPIN_HZ, radius_km=RADIUS_KM[anchor.label])
        beta_eq = spot_speed(J0030_SPIN_HZ, RADIUS_KM[anchor.label],
                             np.pi / 2.0, anchor.compactness)
        out = {"key": "riley" if "Riley" in anchor.label else "miller",
               "u": anchor.compactness, "beta_eq": beta_eq}
        for tag, r in (("0", None), ("rot", rot)):
            iso = multi_spot_flux(anchor.inclination, anchor.compactness, anchor.spots,
                                  None, bending=bend, rotation=r)
            real = multi_spot_flux(anchor.inclination, anchor.compactness, anchor.spots,
                                   beaming, bending=bend, rotation=r)
            out[f"dpf_{tag}"] = pulsed_fraction(real) - pulsed_fraction(iso)
            out[f"floor_{tag}"] = float(iso.min() / iso.max())
        print(f"{anchor.label:<34}{beta_eq:>7.3f}"
              + f"{out['dpf_0']:+.4f}".rjust(12)
              + f"{out['dpf_rot']:+.4f}".rjust(15)
              + f"{out['floor_0']:.4f}".rjust(10)
              + f"{out['floor_rot']:.4f}".rjust(10))
        rows.append(out)
    return rows


def harmonics(flux, n_max: int = 2):
    """Relative harmonic amplitudes [A1/A0, …, A_n/A0] of a phase-periodic curve.

    A_k is the amplitude of the k-cycles-per-rotation Fourier component; A1 is the
    **fundamental** (the basic once-per-rotation wave a waveform fitter leans on),
    A2 the second harmonic (where the Doppler boost injects its power).
    """
    flux = np.asarray(flux, dtype=float)
    coeffs = np.fft.rfft(flux) / flux.size
    a0 = coeffs[0].real
    return [2.0 * float(abs(coeffs[k])) / a0 for k in range(1, n_max + 1)]


def shape_routing(pooled_curve, mu_centers, ti10):
    """Where the beaming systematic lives at each spin: PF vs waveform shape.

    For each J0740 anchor at ν = 0 and 346.5 Hz (pooled library curve): PF_iso,
    PF_real, the peak-normalized iso-vs-real shape difference (RMS, max), and the
    harmonic content A1/A0, A2/A0 of both models. The paper-level read: ΔPF
    collapses ~70% while shape RMS and Δ(A1/A0) barely move — the systematic is
    re-routed out of the min/max statistic, not removed.
    """
    print("\n=== Shape routing (pooled curve): the systematic re-routes, not vanishes ===")
    print(f"{'anchor':<34}{'ν':>10}{'PF_iso':>8}{'PF_real':>8}{'ΔPF':>8}"
          f"{'shapeRMS':>10}{'Δ(A1/A0)':>10}{'Δ(A2/A0)':>10}")
    beaming = beaming_lookup(mu_centers, pooled_curve)
    rows = []
    for anchor in j0740.ANCHORS:
        bend = ExactBending(anchor.compactness)
        key = "riley" if "Riley" in anchor.label else "miller"
        for tag, rot in (("0", None),
                         ("rot", Rotation(spin_hz=J0740_SPIN_HZ,
                                          radius_km=RADIUS_KM[anchor.label]))):
            iso = multi_spot_flux(anchor.inclination, anchor.compactness, anchor.spots,
                                  None, bending=bend, rotation=rot)
            real = multi_spot_flux(anchor.inclination, anchor.compactness, anchor.spots,
                                   beaming, bending=bend, rotation=rot)
            rms, mx = waveform_shape_change(iso, real)
            (a1_i, a2_i), (a1_r, a2_r) = harmonics(iso), harmonics(real)
            row = {"key": key, "tag": tag,
                   "pf_iso": pulsed_fraction(iso), "pf_real": pulsed_fraction(real),
                   "shape_rms": rms, "shape_max": mx,
                   "a1_iso": a1_i, "a1_real": a1_r, "a2_iso": a2_i, "a2_real": a2_r}
            nu_label = "0 Hz" if tag == "0" else f"{J0740_SPIN_HZ:g} Hz"
            print(f"{anchor.label:<34}{nu_label:>10}"
                  f"{row['pf_iso']:>8.4f}{row['pf_real']:>8.4f}"
                  f"{row['pf_real'] - row['pf_iso']:>+8.4f}{rms:>10.4f}"
                  f"{a1_r - a1_i:>+10.4f}{a2_r - a2_i:>+10.4f}")
            rows.append(row)
    return rows


def routing_figure(pooled_curve, mu_centers):
    """Four-pulse figure (Riley J0740): iso/real × 0/346.5 Hz, peak-normalized.

    The visual of the routing claim: between the panels the two curves converge in
    *amplitude* (ΔPF collapses) while the gap between their *shapes* persists.
    """
    anchor = next(a for a in j0740.ANCHORS if "Riley" in a.label)
    bend = ExactBending(anchor.compactness)
    beaming = beaming_lookup(mu_centers, pooled_curve)
    fig, axes = plt.subplots(1, 2, figsize=(10.0, 4.4), sharey=True, layout="constrained")
    for ax, (title, rot) in zip(axes, [
            ("frozen (0 Hz)", None),
            (f"{J0740_SPIN_HZ:g} Hz, δ⁴", Rotation(spin_hz=J0740_SPIN_HZ,
                                                   radius_km=RADIUS_KM[anchor.label]))]):
        iso = multi_spot_flux(anchor.inclination, anchor.compactness, anchor.spots,
                              None, bending=bend, rotation=rot)
        real = multi_spot_flux(anchor.inclination, anchor.compactness, anchor.spots,
                               beaming, bending=bend, rotation=rot)
        phase = np.arange(iso.size) / iso.size
        ax.plot(phase, iso / iso.max(), color="#7f7f7f", lw=1.8, label="isotropic")
        ax.plot(phase, real / real.max(), color="#d62728", lw=1.8,
                label=f"realistic (τ = {SHAPE_TAU:g})")
        dpf = pulsed_fraction(real) - pulsed_fraction(iso)
        rms, _ = waveform_shape_change(iso, real)
        ax.set_title(f"{title}:  ΔPF = {dpf:+.3f},  shape RMS = {rms:.3f}", fontsize=10)
        ax.set_xlabel("rotational phase  φ / 2π  [cycles]")
        ax.margins(x=0.01)
    axes[0].set_ylabel("normalized flux  F(φ) / F_max")
    axes[0].legend(loc="lower left", fontsize=9)
    fig.suptitle("J0740 (Riley): spin drains ΔPF but not the waveform-shape difference "
                 "— the beaming systematic re-routes", fontsize=11)
    fig.savefig(FIGURE_PATH, dpi=130)
    print(f"\n✓ four-pulse routing figure saved to {FIGURE_PATH}")


def main():
    d = np.load(LIBRARY_PATH)
    mu = d["mu_centers"]
    tau = d["tau_values"]
    seed_curves = d["intensity_by_tau_seed"]  # (n_tau, n_seeds, n_bins)
    pooled = d["intensity_by_tau"]
    n_seeds = seed_curves.shape[1]
    ti10 = shape_tau_index(tau)

    print(f"C2 — Doppler coupling to ΔPF(τ={SHAPE_TAU:g}), {n_seeds} production seeds, "
          f"exact bending, δ⁴ energy boost\n")

    j0740_rows = j0740_coupling(seed_curves, mu, ti10)
    j0030_rows = j0030_unmoved(pooled, mu, ti10)
    routing_rows = shape_routing(pooled[ti10], mu, ti10)
    routing_figure(pooled[ti10], mu)

    # --- Gate G2 verdict on the leading (most compact) anchor ------------------
    riley = next(r for r in j0740_rows if r["key"] == "riley")
    correct = riley["ratio"] > 1.0
    print("\n" + "=" * 78)
    print(f"Gate G2 (Riley, u = {riley['u']:.3f}, β_eq = {riley['beta_eq']:.3f} at "
          f"{J0740_SPIN_HZ:g} Hz):")
    print(f"  ΔPF(0) = {riley['dpf0_m']:+.4f} ± {riley['dpf0_s']:.4f}   "
          f"ΔPF(ν) = {riley['dpfr_m']:+.4f} ± {riley['dpfr_s']:.4f}   "
          f"coupling = {riley['coupling']:+.4f} ({riley['ratio']:.1f}σ)")
    if correct:
        print("  |coupling| > σ  →  report ΔPF(ν) as a quantified CORRECTION alongside "
              "the headline\n  (the Doppler-corrected number, with oblateness/time-delay "
              "as stated caveats).")
    else:
        print("  |coupling| ≤ σ  →  one quantified CAVEAT sentence: Doppler shifts ΔPF by "
              f"< {riley['sigma']:.3f}\n  (below the seed error bar); the frozen headline "
              "stands.")
    j0030_ok = all(abs(r["dpf_rot"]) < 1e-3 for r in j0030_rows)
    print(f"  J0030 saturation preserved under rotation: "
          f"{'YES (ΔPF stays ≈ 0)' if j0030_ok else 'NO — investigate'}")

    payload = {"tau_shape": np.float64(SHAPE_TAU), "n_seeds": np.int64(n_seeds),
               "j0740_spin_hz": np.float64(J0740_SPIN_HZ),
               "j0030_spin_hz": np.float64(J0030_SPIN_HZ),
               "gate_correct": np.bool_(correct)}
    for r in j0740_rows:
        k = r["key"]
        for f in ("u", "beta_eq", "dpf0_m", "dpf0_s", "dpfr_m", "dpfr_s",
                  "coupling", "sigma", "ratio"):
            payload[f"j0740_{k}_{f}"] = np.float64(r[f])
    for r in j0030_rows:
        k = r["key"]
        payload[f"j0030_{k}_dpf0"] = np.float64(r["dpf_0"])
        payload[f"j0030_{k}_dpf_rot"] = np.float64(r["dpf_rot"])
    for r in routing_rows:
        stem = f"j0740_{r['key']}"
        for f in ("pf_iso", "pf_real", "shape_rms", "shape_max",
                  "a1_iso", "a1_real", "a2_iso", "a2_real"):
            payload[f"{stem}_{f}_{r['tag']}"] = np.float64(r[f])
    np.savez(RESULTS_PATH, **payload)
    print(f"\n✓ C2 results saved to {RESULTS_PATH}")


if __name__ == "__main__":
    main()
