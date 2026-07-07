"""Track A3: per-seed downstream error bars — the ΔPF ± σ the paper quotes.

The production library (v0.9.7.1) ships 5 independent per-seed I(μ; τ) curves in
``intensity_by_tau_seed``. Pushing each seed's curve through the *same* verified
anchor pipeline that produces the pooled headline yields the seed-to-seed spread
of the observable — i.e. the Monte Carlo error bar on ΔPF, in the currency of the
paper's claim. Pure interpolation + geometry (no Monte Carlo here), so it costs
seconds.

Reuses ``j0740_anchor.sweep_anchor`` / ``j0030_anchor.sweep_anchor`` unchanged: a
single seed's ``intensity_by_tau_seed[:, si, :]`` slice already has the
(n_tau, n_bins) shape those functions expect. The pooled-curve result is quoted
alongside as the central value.

Run from the repository root:  PYTHONPATH=src:scripts python3 scripts/a3_seed_errors.py
"""

import numpy as np

from anchor_lib import shape_tau_index
import j0030_anchor as j0030
import j0740_anchor as j0740

LIBRARY_PATH = "data/beaming_library.npz"
RESULTS_PATH = "data/a3_seed_errors.npz"

# v0.9.7 convergence-study converged reference (Riley J0740 ΔPF at τ = 10), the
# gate this A3 must reproduce within ~2σ (docs/deep-dives/v0.9.7-convergence-redo.md §3.3).
V097_DPF_RILEY_TAU10 = (0.1486, 0.0026)


def per_seed_metric(module, anchor, tau, mu, seed_curves, metric_key):
    """(mean, std) across seeds of a sweep_anchor metric, per τ.

    ``seed_curves`` is intensity_by_tau_seed with shape (n_tau, n_seeds, n_bins);
    seed ``si`` is passed as the (n_tau, n_bins) slice the anchor already expects.
    """
    n_seeds = seed_curves.shape[1]
    vals = np.array([
        module.sweep_anchor(anchor, tau, mu, seed_curves[:, si, :])[metric_key]
        for si in range(n_seeds)
    ])  # (n_seeds, n_tau)
    return vals.mean(axis=0), vals.std(axis=0, ddof=1)


def main():
    d = np.load(LIBRARY_PATH)
    tau, mu = d["tau_values"], d["mu_centers"]
    seed_curves = d["intensity_by_tau_seed"]
    pooled = d["intensity_by_tau"]
    n_seeds = seed_curves.shape[1]
    ti10 = shape_tau_index(tau)  # τ = 10, the headline optical depth

    payload = {"tau_values": tau, "n_seeds": np.int64(n_seeds)}
    print(f"A3 — downstream ΔPF ± σ from {n_seeds} per-seed library curves\n")

    # --- J0740: the live pulsed-fraction systematic (the headline) ---
    print("=== J0740 (non-eclipsing): ΔPF(τ) = PF_real − PF_iso ===")
    print(f"{'anchor':<34}{'ΔPF(τ=10) ± σ':>22}{'ΔPF(τ) peak':>18}")
    for anchor in j0740.ANCHORS:
        m, s = per_seed_metric(j0740, anchor, tau, mu, seed_curves, "delta_pf")
        pooled_dpf = j0740.sweep_anchor(anchor, tau, mu, pooled)["delta_pf"]
        pk = int(np.argmax(m))
        print(f"{anchor.label:<34}{m[ti10]:+.4f} ± {s[ti10]:.4f}"
              f"      {m[pk]:+.4f} @ τ={tau[pk]:g}")
        key = "riley" if "Riley" in anchor.label else "miller"
        payload[f"j0740_{key}_dpf_mean"] = m
        payload[f"j0740_{key}_dpf_std"] = s
        payload[f"j0740_{key}_dpf_pooled"] = pooled_dpf

    # Gate: Riley ΔPF(τ=10) must land within ~2σ of the v0.9.7 converged value.
    rm, rs = payload["j0740_riley_dpf_mean"], payload["j0740_riley_dpf_std"]
    ref, ref_s = V097_DPF_RILEY_TAU10
    z = abs(rm[ti10] - ref) / np.hypot(rs[ti10], ref_s)
    print(f"\nGate: Riley ΔPF(τ=10) = {rm[ti10]:+.4f} ± {rs[ti10]:.4f}  vs "
          f"v0.9.7 {ref:+.4f} ± {ref_s:.4f}  → {z:.1f}σ  "
          f"({'PASS' if z <= 2 else 'CHECK'})")

    # --- J0030: ΔPF saturates (PF pinned), so the systematic is in waveform shape ---
    print("\n=== J0030 (eclipsing): waveform-shape RMS(τ) (ΔPF ≈ 0 → shape) ===")
    for anchor in j0030.ANCHORS:
        m, s = per_seed_metric(j0030, anchor, tau, mu, seed_curves, "shape_rms")
        pk = int(np.argmax(m))
        print(f"{anchor.label:<34} peak RMS {m[pk]:.4f} ± {s[pk]:.4f} @ τ={tau[pk]:g}")
        key = "riley" if "Riley" in anchor.label else "miller"
        payload[f"j0030_{key}_shape_mean"] = m
        payload[f"j0030_{key}_shape_std"] = s

    np.savez(RESULTS_PATH, **payload)
    print(f"\n✓ per-seed error bars saved to {RESULTS_PATH}")


if __name__ == "__main__":
    main()
