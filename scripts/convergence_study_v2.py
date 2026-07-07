"""Per-τ photon-count convergence study (v2) — sizing the beaming library.

The v0.7.0 study (scripts/convergence_study.py) measured error-vs-N at a single
optical depth (τ = 10) and recommended 200k injected photons — a number the
library (scripts/tau_sweep.py) then applied uniformly across all six τ. That
silently ignores that the histogram is built only from *escaped* photons, and
the escape fraction falls from ~92% (τ = 0.1) to ~4% (τ = 30): the thickest,
most physically interesting rows carry ~22× worse statistics than the row the
study validated. This study fixes the three v1 limitations documented in
docs/next-steps.md §5:

  1. **Every library τ is swept**, not just τ = 10.
  2. **Convergence is measured against escaped photons** — the real statistical
     currency — so one calibration converts to an injected-N budget per τ.
  3. **ΔPF joins the tracked observables.** The paper's claim is the pulsed-
     fraction shift ΔPF ≈ +0.15 (J0740, Riley 2021 geometry); each (τ, N, seed)
     run's I(μ) is pushed through the anchor pipeline so we measure directly
     how much of ΔPF's spread is Monte Carlo noise. The mean should hold steady
     with N while the seed-to-seed spread collapses as N^(-1/2).

Design notes:
  * 10 seeds for N ≤ 1e5, 5 above — a 5-sample std is itself ~35% uncertain
    (1/√(2(n−1))), which is why the v1 error column wobbled non-monotonically.
  * Runs are independent across (τ, N, seed), so the sweep parallelizes over
    processes; RNG streams come from SeedSequence spawn keys derived from the
    (τ, N, seed) coordinates, so results are reproducible run-to-run and
    insensitive to grid edits or task scheduling order.
  * τ = 10 and 30 get an extra N = 3e6 point so the production readoff for the
    expensive thick rows rests on measurement, not extrapolation.

Production targets (5-seed pooled library, pooled error = per-seed / √5):
    pooled σ(b)        ≤ 0.02  → per-seed target 0.0447
    pooled tail-bin σ  ≤ 0.02  → per-seed target 0.0447 (absolute, I normalized to μ≈1)
    pooled σ(ΔPF)      ≤ 0.01  → per-seed target 0.0224

Usage (repository root):
    PYTHONPATH=src python3 scripts/convergence_study_v2.py            # full sweep
    PYTHONPATH=src python3 scripts/convergence_study_v2.py --quick    # plumbing check
    PYTHONPATH=src python3 scripts/convergence_study_v2.py --summarize-only
"""

import argparse
import os
import time
from concurrent.futures import ProcessPoolExecutor, as_completed

import numpy as np

from mcrt import (
    Simulation,
    beaming_lookup,
    extract_intensity,
    fit_limb_darkening_slope,
    loglog_slope,
    n_for_target_error,
    pulsed_fraction,
)

TAU_VALUES = [0.1, 0.3, 1.0, 3.0, 10.0, 30.0]
BASE_N_VALUES = [1_000, 3_000, 10_000, 30_000, 100_000, 300_000, 1_000_000]
# Thick rows are where the production budget will land — anchor them with a
# measured point instead of extrapolating the -1/2 line past the swept range.
EXTRA_N = {10.0: [3_000_000], 30.0: [3_000_000]}
SEEDS_SMALL_N = 10   # N ≤ SEED_BOOST_MAX_N
SEEDS_LARGE_N = 5
SEED_BOOST_MAX_N = 100_000

BASE_SEED = 20260706
N_BINS = 20
RESULTS_PATH = "data/convergence_v2_results.npz"
ERROR_PLOT_PATH = "data/convergence_v2_error_vs_escaped.png"
DPF_PLOT_PATH = "data/convergence_v2_dpf.png"

# Per-seed error targets (see module docstring for the pooled-library logic).
POOL_SEEDS = 5
TARGET_SLOPE_B = 0.02 * np.sqrt(POOL_SEEDS)
TARGET_INTENSITY = 0.02 * np.sqrt(POOL_SEEDS)
TARGET_DPF = 0.01 * np.sqrt(POOL_SEEDS)

# Rough s-per-1e6-injected by τ (measured 2026-07-06, scalar engine, M-series
# laptop) — used only to schedule long tasks first, never for physics.
COST_PER_M = {0.1: 4.2, 0.3: 5.8, 1.0: 12.2, 3.0: 28.3, 10.0: 85.3, 30.0: 252.0}

QUICK_TAU_VALUES = [0.1, 30.0]
QUICK_N_VALUES = [1_000, 10_000, 100_000]
QUICK_SEEDS = 3


def n_seeds_for(n):
    return SEEDS_SMALL_N if n <= SEED_BOOST_MAX_N else SEEDS_LARGE_N


def build_tasks(tau_values=TAU_VALUES, base_n=BASE_N_VALUES, extra_n=EXTRA_N,
                seeds_fn=n_seeds_for):
    """Flat (τ, N, seed) task list, longest-cost first for load balancing."""
    tasks = []
    for tau in tau_values:
        for n in list(base_n) + list(extra_n.get(tau, [])):
            for seed in range(seeds_fn(n)):
                tasks.append((tau, n, seed))
    return sorted(tasks, key=lambda t: t[1] * COST_PER_M.get(t[0], 100.0),
                  reverse=True)


def run_one(task):
    """One (τ, N, seed) simulation → observables. Runs in a worker process.

    RNG stream: SeedSequence(BASE_SEED, spawn_key=(τ×10, N, seed)) — keyed on
    the task coordinates, so a run is reproducible regardless of which worker
    executes it or how the grid is later extended.
    """
    tau, n, seed = task
    ss = np.random.SeedSequence(BASE_SEED, spawn_key=(int(round(tau * 10)), int(n), seed))
    rng = np.random.default_rng(ss)

    t0 = time.perf_counter()
    sim = Simulation(tau_total=tau, num_photons=n, rng=rng)
    sim.run()
    runtime = time.perf_counter() - t0

    escaped_mu = sim.results["escaped_mu"]
    escaped = len(escaped_mu)

    # Undersampled corner (few escapers / empty normalization bin): record NaNs
    # rather than a garbage curve; the analysis masks them.
    with np.errstate(divide="ignore", invalid="ignore"):
        mu_centers, intensity = extract_intensity(escaped_mu, n_bins=N_BINS)
    if not np.all(np.isfinite(intensity)):
        intensity = np.full(N_BINS, np.nan)
        slope = np.nan
    else:
        slope = fit_limb_darkening_slope(mu_centers, intensity)

    return {
        "tau": tau, "n": n, "seed": seed,
        "escaped": escaped,
        "slope": slope,
        "intensity": intensity,
        "mu_centers": mu_centers,
        "runtime": runtime,
    }


def compute_delta_pf(mu_centers, intensity_rows):
    """ΔPF (realistic − isotropic PF) for each I(μ) row, Riley J0740 geometry.

    This is the paper's headline observable: the anchor pipeline is pure
    interpolation + geometry, so pushing every (τ, N, seed) beaming curve
    through it costs nothing and measures the Monte Carlo noise *in the
    currency of the claim*. NaN curves (undersampled runs) → NaN ΔPF.
    """
    from j0740_anchor import RILEY          # published geometry, single home
    from anchor_lib import multi_spot_flux  # verified two-spot mechanics

    iso_flux = multi_spot_flux(RILEY.inclination, RILEY.compactness, RILEY.spots,
                               beaming=None)
    pf_iso = pulsed_fraction(iso_flux)

    dpf = np.full(len(intensity_rows), np.nan)
    for i, intensity in enumerate(intensity_rows):
        if not np.all(np.isfinite(intensity)):
            continue
        beaming = beaming_lookup(mu_centers, intensity)
        real_flux = multi_spot_flux(RILEY.inclination, RILEY.compactness,
                                    RILEY.spots, beaming=beaming)
        dpf[i] = pulsed_fraction(real_flux) - pf_iso
    return dpf


def run_sweep(tasks, max_workers):
    """Execute all tasks in parallel; return flat record arrays."""
    records = []
    total = len(tasks)
    done = 0
    with ProcessPoolExecutor(max_workers=max_workers) as pool:
        futures = [pool.submit(run_one, t) for t in tasks]
        for fut in as_completed(futures):
            records.append(fut.result())
            done += 1
            if done % 25 == 0 or done == total:
                print(f"  {done}/{total} runs complete", flush=True)

    mu_centers = next(r["mu_centers"] for r in records)
    flat = {
        "tau": np.array([r["tau"] for r in records]),
        "n": np.array([r["n"] for r in records], dtype=float),
        "seed": np.array([r["seed"] for r in records]),
        "escaped": np.array([r["escaped"] for r in records], dtype=float),
        "slope": np.array([r["slope"] for r in records]),
        "intensity": np.array([r["intensity"] for r in records]),
        "runtime": np.array([r["runtime"] for r in records]),
        "mu_centers": mu_centers,
    }
    flat["delta_pf"] = compute_delta_pf(mu_centers, flat["intensity"])
    return flat


def summarize(flat):
    """Group by (τ, N): mean escaped + across-seed error per observable."""
    mu = flat["mu_centers"]
    tail_idx = 0
    bulk_idx = int(np.argmin(np.abs(mu - 0.7)))

    out = {}
    for tau in sorted(set(flat["tau"])):
        sel_tau = flat["tau"] == tau
        n_values = np.array(sorted(set(flat["n"][sel_tau])))
        rows = {"n": n_values,
                "escaped_mean": np.zeros(len(n_values)),
                "slope_mean": np.zeros(len(n_values)),
                "slope_err": np.zeros(len(n_values)),
                "tail_err": np.zeros(len(n_values)),
                "bulk_err": np.zeros(len(n_values)),
                "dpf_mean": np.zeros(len(n_values)),
                "dpf_err": np.zeros(len(n_values))}
        for i, n in enumerate(n_values):
            sel = sel_tau & (flat["n"] == n)
            rows["escaped_mean"][i] = flat["escaped"][sel].mean()
            slopes = flat["slope"][sel]
            dpfs = flat["delta_pf"][sel]
            tail = flat["intensity"][sel][:, tail_idx]
            bulk = flat["intensity"][sel][:, bulk_idx]
            # NaN-aware: undersampled seeds are masked, not averaged in.
            rows["slope_mean"][i] = np.nanmean(slopes) if np.isfinite(slopes).any() else np.nan
            rows["slope_err"][i] = _nanstd(slopes)
            rows["tail_err"][i] = _nanstd(tail)
            rows["bulk_err"][i] = _nanstd(bulk)
            rows["dpf_mean"][i] = np.nanmean(dpfs) if np.isfinite(dpfs).any() else np.nan
            rows["dpf_err"][i] = _nanstd(dpfs)
        # Escape fraction from the largest-N runs (tightest estimate).
        rows["escape_frac"] = rows["escaped_mean"][-1] / n_values[-1]
        out[tau] = rows
    return out, {"tail_idx": tail_idx, "bulk_idx": bulk_idx}


def _nanstd(values):
    finite = np.asarray(values)[np.isfinite(values)]
    if len(finite) < 2:
        return np.nan
    return float(np.std(finite, ddof=1))


def print_report(summary):
    """Per-τ fitted slopes, escaped-N targets, and the production budget."""
    print("\n=== Convergence in escaped-photon currency ===")
    header = (f"{'τ':>5} {'f_esc':>7} | {'slope b':>8} {'tail':>8} {'ΔPF':>8} "
              f"| {'N_esc(b)':>10} {'N_esc(tail)':>12} {'N_esc(ΔPF)':>11} "
              f"| {'binding N_esc':>13} {'N_inject/seed':>14}")
    print(header)
    print("-" * len(header))

    budget = {}
    for tau, rows in summary.items():
        esc = rows["escaped_mean"]
        targets = [("slope_err", TARGET_SLOPE_B), ("tail_err", TARGET_INTENSITY),
                   ("dpf_err", TARGET_DPF)]
        slopes, needs = [], []
        for key, target in targets:
            slopes.append(loglog_slope(esc, rows[key]))
            rec = n_for_target_error(esc, rows[key], target)
            needs.append(rec.n if rec is not None else np.nan)
        binding = np.nanmax(needs)
        inject = binding / rows["escape_frac"]
        budget[tau] = {"binding_escaped": binding, "inject_per_seed": inject,
                       "escape_frac": rows["escape_frac"]}
        print(f"{tau:>5} {rows['escape_frac']:>7.3f} "
              f"| {slopes[0]:>8.2f} {slopes[1]:>8.2f} {slopes[2]:>8.2f} "
              f"| {needs[0]:>10.2e} {needs[1]:>12.2e} {needs[2]:>11.2e} "
              f"| {binding:>13.2e} {inject:>14.2e}")

    print(f"\nTargets (per-seed, {POOL_SEEDS}-seed pooled library): "
          f"σ(b) ≤ {TARGET_SLOPE_B:.3f}, tail-bin σ ≤ {TARGET_INTENSITY:.3f}, "
          f"σ(ΔPF) ≤ {TARGET_DPF:.3f}")
    print("N_esc columns: escaped photons per seed needed to meet each target "
          "(fitted -1/2 line); N_inject = binding / escape fraction.")
    return budget


def print_dpf_table(summary):
    print("\n=== ΔPF (Riley J0740): mean holds, spread collapses ===")
    for tau, rows in summary.items():
        cells = "  ".join(
            f"{m:+.3f}±{e:.3f}" if np.isfinite(m) and np.isfinite(e) else "   n/a    "
            for m, e in zip(rows["dpf_mean"], rows["dpf_err"]))
        print(f"  τ={tau:>4}: {cells}")
    any_tau = next(iter(summary.values()))
    print(f"  (columns are N = {[f'{v:.0e}' for v in any_tau['n']]}"
          " — thick rows have one extra 3e6 column)")


def plot_error_vs_escaped(summary, path=ERROR_PLOT_PATH):
    """2×2 panels: per-observable error vs escaped photons, one curve per τ.

    If curves from different τ collapse onto one line per panel, escaped count
    is the universal statistical currency and the per-τ budget conversion holds.
    """
    import matplotlib.pyplot as plt
    from matplotlib import cm

    panels = [("slope_err", "limb-slope b error", TARGET_SLOPE_B),
              ("tail_err", "tail-bin (μ≈0.03) error", TARGET_INTENSITY),
              ("bulk_err", "bulk-bin (μ≈0.68) error", TARGET_INTENSITY),
              ("dpf_err", "σ(ΔPF), Riley J0740", TARGET_DPF)]
    taus = list(summary.keys())
    colors = cm.viridis(np.linspace(0.0, 0.9, len(taus)))

    fig, axes = plt.subplots(2, 2, figsize=(12, 9))
    for ax, (key, label, target) in zip(axes.flat, panels):
        for tau, color in zip(taus, colors):
            rows = summary[tau]
            esc, err = rows["escaped_mean"], rows[key]
            mask = np.isfinite(err) & (err > 0)
            ax.loglog(esc[mask], err[mask], "o-", color=color, ms=4,
                      label=f"τ = {tau}")
        ref_esc = np.array([1e2, 3e6])
        anchor_rows = summary[taus[-1]]
        fin = np.isfinite(anchor_rows[key]) & (anchor_rows[key] > 0)
        if fin.any():
            e0, r0 = anchor_rows["escaped_mean"][fin][-1], anchor_rows[key][fin][-1]
            ax.loglog(ref_esc, r0 * (ref_esc / e0) ** -0.5, "k:", alpha=0.6,
                      label=r"$N_{\rm esc}^{-1/2}$")
        ax.axhline(target, color="#c0392b", ls="--", lw=1,
                   label=f"per-seed target {target:.3f}")
        ax.set_xlabel("Escaped photons")
        ax.set_ylabel("Error (std across seeds)")
        ax.set_title(label)
        ax.grid(True, which="both", alpha=0.3)
        ax.legend(fontsize=7)
    fig.suptitle("Convergence v2: all τ collapse in escaped-photon currency")
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)
    print(f"✓ Error-vs-escaped plot saved to {path}")


def plot_dpf(summary, path=DPF_PLOT_PATH):
    """Mean ΔPF ± seed spread vs escaped photons: signal steady, noise collapsing."""
    import matplotlib.pyplot as plt
    from matplotlib import cm

    taus = list(summary.keys())
    colors = cm.viridis(np.linspace(0.0, 0.9, len(taus)))

    plt.figure(figsize=(9, 6))
    for tau, color in zip(taus, colors):
        rows = summary[tau]
        mask = np.isfinite(rows["dpf_mean"]) & np.isfinite(rows["dpf_err"])
        plt.errorbar(rows["escaped_mean"][mask], rows["dpf_mean"][mask],
                     yerr=rows["dpf_err"][mask], fmt="o-", ms=4, capsize=3,
                     color=color, label=f"τ = {tau}")
    plt.xscale("log")
    plt.xlabel("Escaped photons")
    plt.ylabel("ΔPF = PF(realistic) − PF(isotropic)")
    plt.title("ΔPF (Riley J0740 geometry): mean vs Monte Carlo noise")
    plt.grid(True, alpha=0.3)
    plt.legend(fontsize=8)
    plt.savefig(path)
    plt.close()
    print(f"✓ ΔPF convergence plot saved to {path}")


def save_flat(flat, path=RESULTS_PATH):
    np.savez(path, **flat)
    print(f"✓ Raw sweep records saved to {path}")


def load_flat(path=RESULTS_PATH):
    data = np.load(path)
    return {k: data[k] for k in data.files}


def main():
    parser = argparse.ArgumentParser(description="Per-τ convergence study (v2).")
    parser.add_argument("--quick", action="store_true",
                        help="Tiny sweep (2 τ, 3 N, 3 seeds) for plumbing checks.")
    parser.add_argument("--summarize-only", action="store_true",
                        help=f"Skip the sweep; re-report/plot from {RESULTS_PATH}.")
    parser.add_argument("--workers", type=int,
                        default=max(1, (os.cpu_count() or 4) - 2))
    args = parser.parse_args()

    if args.summarize_only:
        print(f"Loading cached sweep from {RESULTS_PATH}...")
        flat = load_flat()
    else:
        if args.quick:
            tasks = build_tasks(QUICK_TAU_VALUES, QUICK_N_VALUES, {},
                                seeds_fn=lambda n: QUICK_SEEDS)
            print(f"QUICK sweep: {len(tasks)} runs on {args.workers} workers...")
        else:
            tasks = build_tasks()
            total_photons = sum(t[1] for t in tasks)
            est = sum(t[1] * COST_PER_M.get(t[0], 100.0) / 1e6 for t in tasks)
            print(f"FULL sweep: {len(tasks)} runs, {total_photons/1e6:.1f}M photons, "
                  f"~{est/60:.0f} CPU-min on {args.workers} workers "
                  f"(~{est/60/args.workers:.0f} min wall if balanced)...")
        t0 = time.perf_counter()
        flat = run_sweep(tasks, max_workers=args.workers)
        print(f"Sweep finished in {(time.perf_counter() - t0)/60:.1f} min wall.")
        save_flat(flat)

    summary, _ = summarize(flat)
    print_dpf_table(summary)
    print_report(summary)
    plot_error_vs_escaped(summary)
    plot_dpf(summary)


if __name__ == "__main__":
    main()
