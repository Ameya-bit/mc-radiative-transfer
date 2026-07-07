"""Track B2 — the measured linear→exact bending shift in the headline ΔPF.

B1 (v0.9.9) built the exact Schwarzschild map ``ExactBending(u)`` and validated it
(SD1a code-comparison residual 0.80% → 0.11%). B2 answers the question the whole
track exists for: **how much was Beloborodov's linear approximation biasing the
production ΔPF at J0740**, where the faint phase lives at grazing emission and the
compactness (u = 0.494 Riley, 0.444 Miller) sits at the edge of the linear map's
stated validity?

This driver re-runs both J0740 anchors through the *same* per-seed pipeline that
produced the linear headline (A3), once with the linear map and once with exact
bending, on the production 5-seed library. For each it reports

    ΔPF_linear ± σ,  ΔPF_exact ± σ,  shift = ΔPF_exact − ΔPF_linear,  |shift| / σ,

with σ the seed-to-seed (Monte Carlo) error bar. The shift is a **systematic**
(a fixed modeling bias; more photons shrink σ but not the shift), so Gate G1 keys
on |shift| vs σ, not on sign:

    |shift| ≤ σ  → keep Beloborodov everywhere; quote the check in one sentence.
    |shift| >  σ → the linear map biases the headline beyond its own error bar;
                   report the exact number for the anchors + phase-diagram markers.

It also confirms the eclipsing anchor J0030 (u = 0.312/0.326, PF saturated) is
**unmoved** by the bending swap — PF stays pinned regardless of the map, so the
systematic there is still shape-only, and no anchor decision hangs on it.

Pure interpolation + geometry (no Monte Carlo), so it costs seconds. The linear
branch reproduces A3 bit-for-bit (``bending=None``); only the exact branch is new.

Run from the repository root:  PYTHONPATH=src:scripts python3 scripts/b2_exact_bending.py
"""

import numpy as np

from mcrt import ExactBending, beaming_lookup, pulsed_fraction
from anchor_lib import SHAPE_TAU, multi_spot_flux, shape_tau_index
import j0030_anchor as j0030
import j0740_anchor as j0740

LIBRARY_PATH = "data/beaming_library.npz"
RESULTS_PATH = "data/b2_exact_bending.npz"


def per_seed_delta_pf(anchor, mu_centers, seed_curves_at_tau, bending):
    """(mean, std) of ΔPF = PF_real − PF_iso across seeds, for one bending map.

    ``seed_curves_at_tau`` is ``intensity_by_tau_seed[ti]`` with shape
    (n_seeds, n_bins) — the per-seed I(μ) at the headline τ. Both PF_iso and each
    PF_real use the *same* ``bending`` map, so the returned spread is the seed MC
    error and the exact−linear difference is the pure systematic.
    """
    iso = multi_spot_flux(anchor.inclination, anchor.compactness, anchor.spots,
                          None, bending=bending)
    pf_iso = pulsed_fraction(iso)

    n_seeds = seed_curves_at_tau.shape[0]
    dpf = np.empty(n_seeds)
    for si in range(n_seeds):
        beaming = beaming_lookup(mu_centers, seed_curves_at_tau[si])
        real = multi_spot_flux(anchor.inclination, anchor.compactness,
                               anchor.spots, beaming, bending=bending)
        dpf[si] = pulsed_fraction(real) - pf_iso
    return dpf.mean(), dpf.std(ddof=1), pf_iso


def j0740_shift(seed_curves, mu_centers, ti10):
    """Per-anchor linear vs exact ΔPF(τ=10) ± σ and the gate verdict."""
    rows = []
    print("=== J0740 (non-eclipsing): ΔPF(τ=10) = PF_real − PF_iso ===")
    print(f"{'anchor':<34}{'ΔPF_linear':>14}{'ΔPF_exact':>14}"
          f"{'shift':>10}{'|shift|/σ':>11}")
    for anchor in j0740.ANCHORS:
        bend = ExactBending(anchor.compactness)
        curves = seed_curves[ti10]  # (n_seeds, n_bins)
        lin_m, lin_s, lin_pfiso = per_seed_delta_pf(anchor, mu_centers, curves, None)
        exa_m, exa_s, exa_pfiso = per_seed_delta_pf(anchor, mu_centers, curves, bend)
        shift = exa_m - lin_m
        sigma = max(lin_s, exa_s)  # judge against the (comparable) seed error
        ratio = abs(shift) / sigma if sigma > 0 else np.inf
        key = "riley" if "Riley" in anchor.label else "miller"
        print(f"{anchor.label:<34}"
              + f"{lin_m:+.4f}±{lin_s:.4f}".rjust(14)
              + f"{exa_m:+.4f}±{exa_s:.4f}".rjust(14)
              + f"{shift:+.4f}".rjust(10)
              + f"{ratio:.1f}σ".rjust(11))
        rows.append({
            "key": key, "u": anchor.compactness,
            "lin_m": lin_m, "lin_s": lin_s, "exa_m": exa_m, "exa_s": exa_s,
            "shift": shift, "sigma": sigma, "ratio": ratio,
            "pf_iso_lin": lin_pfiso, "pf_iso_exa": exa_pfiso,
        })
    return rows


def j0030_unmoved(pooled, mu_centers, ti10):
    """Confirm the eclipsing anchor's ΔPF stays ≈ 0 (saturated) under exact bending.

    Uses the pooled library curve — the point is qualitative (PF stays pinned, so
    the shape-only verdict is bending-independent), not a per-seed error bar.
    """
    print("\n=== J0030 (eclipsing): PF saturated → bending must not move ΔPF ===")
    print(f"{'anchor':<34}{'ΔPF_linear':>12}{'ΔPF_exact':>12}"
          f"{'floor_lin':>11}{'floor_exa':>11}")
    beaming = beaming_lookup(mu_centers, pooled[ti10])
    rows = []
    for anchor in j0030.ANCHORS:
        bend = ExactBending(anchor.compactness)
        out = {}
        for tag, b in (("lin", None), ("exa", bend)):
            iso = multi_spot_flux(anchor.inclination, anchor.compactness,
                                  anchor.spots, None, bending=b)
            real = multi_spot_flux(anchor.inclination, anchor.compactness,
                                   anchor.spots, beaming, bending=b)
            out[f"dpf_{tag}"] = pulsed_fraction(real) - pulsed_fraction(iso)
            out[f"floor_{tag}"] = float(iso.min() / iso.max())
        key = "riley" if "Riley" in anchor.label else "miller"
        print(f"{anchor.label:<34}"
              + f"{out['dpf_lin']:+.4f}".rjust(12)
              + f"{out['dpf_exa']:+.4f}".rjust(12)
              + f"{out['floor_lin']:.4f}".rjust(11)
              + f"{out['floor_exa']:.4f}".rjust(11))
        rows.append({"key": key, "u": anchor.compactness, **out})
    return rows


def main():
    d = np.load(LIBRARY_PATH)
    mu = d["mu_centers"]
    tau = d["tau_values"]
    seed_curves = d["intensity_by_tau_seed"]  # (n_tau, n_seeds, n_bins)
    pooled = d["intensity_by_tau"]
    n_seeds = seed_curves.shape[1]
    ti10 = shape_tau_index(tau)

    print(f"B2 — linear→exact bending shift in ΔPF(τ={SHAPE_TAU:g}), "
          f"{n_seeds} production seeds\n")

    j0740_rows = j0740_shift(seed_curves, mu, ti10)
    j0030_rows = j0030_unmoved(pooled, mu, ti10)

    # --- Gate G1 verdict on the leading (most compact) anchor ------------------
    riley = next(r for r in j0740_rows if r["key"] == "riley")
    switch = riley["ratio"] > 1.0
    print("\n" + "=" * 70)
    print(f"Gate G1 (Riley, u = {riley['u']:.3f}): "
          f"|shift| = {abs(riley['shift']):.4f}, σ = {riley['sigma']:.4f} → "
          f"{riley['ratio']:.1f}σ")
    if switch:
        print("  |shift| > σ  →  SWITCH anchors + phase-diagram markers to exact "
              "bending;\n  the broad phase-diagram sweep stays linear with a stated "
              "tolerance.")
    else:
        print("  |shift| ≤ σ  →  KEEP Beloborodov everywhere; quote the check in "
              "one sentence.")
    j0030_ok = all(abs(r["dpf_exa"]) < 1e-3 for r in j0030_rows)
    print(f"  J0030 saturation preserved under exact bending: "
          f"{'YES (ΔPF stays ≈ 0)' if j0030_ok else 'NO — investigate'}")

    payload = {"tau_shape": np.float64(SHAPE_TAU), "n_seeds": np.int64(n_seeds),
               "gate_switch": np.bool_(switch)}
    for r in j0740_rows:
        k = r["key"]
        payload[f"j0740_{k}_u"] = np.float64(r["u"])
        payload[f"j0740_{k}_dpf_linear_mean"] = np.float64(r["lin_m"])
        payload[f"j0740_{k}_dpf_linear_std"] = np.float64(r["lin_s"])
        payload[f"j0740_{k}_dpf_exact_mean"] = np.float64(r["exa_m"])
        payload[f"j0740_{k}_dpf_exact_std"] = np.float64(r["exa_s"])
        payload[f"j0740_{k}_shift"] = np.float64(r["shift"])
        payload[f"j0740_{k}_shift_over_sigma"] = np.float64(r["ratio"])
    for r in j0030_rows:
        k = r["key"]
        payload[f"j0030_{k}_dpf_linear"] = np.float64(r["dpf_lin"])
        payload[f"j0030_{k}_dpf_exact"] = np.float64(r["dpf_exa"])
    np.savez(RESULTS_PATH, **payload)
    print(f"\n✓ B2 results saved to {RESULTS_PATH}")


if __name__ == "__main__":
    main()
