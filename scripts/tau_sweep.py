"""Build the production beaming-function library I(μ; τ) across optical depths.

The emergent specific intensity I(μ) carries the limb darkening set by how many
times a photon scatters before escaping:

    τ ≪ 1 (thin)  → little scattering → nearly isotropic emission, slope b → 0
    τ ≫ 1 (thick) → many scatterings  → approaches Eddington / Chandrasekhar, b → ~1.5

Downstream pulse-profile code interpolates I(μ) for a chosen τ from the saved
lookup table (data/beaming_library.npz) instead of re-running the Monte Carlo.

This is the **v0.9.8 production rework** (next-steps.md Track A1), replacing the
uniform-200k-injected build with the escape-matched, multi-seed design the v0.9.7
convergence redo prescribed (docs/deep-dives/v0.9.7-convergence-redo.md §3.4):

  * **Escape-matched injected N per (τ, seed).** The I(μ) histogram is built only
    from *escaped* photons, and the escape fraction falls from ~92% (τ=0.1) to
    ~4% (τ=30). A fixed injected count is therefore 22× worse statistics at τ=30
    than at τ=0.1. We instead target a uniform 4×10⁵ *escaped* photons per row,
    injecting more where the slab is thicker (≈ 80M photons total).
  * **5 independent seeds per τ**, each on its own reproducible RNG stream
    (SeedSequence spawn keys), so the library ships per-bin seed error bars and
    b(τ) ± σ — the referee's first question about the ΔPF headline.
  * **Fix 8a (normalization):** each I(μ) curve is normalized by the *fitted*
    intercept a of I = a(1 + bμ), not a single noisy edge bin, so the stored curve
    does not inherit one bin's Monte Carlo noise.
  * **Fix 8b (bin-center bias):** counts are divided by the *mean escaped μ* in
    each bin (pooled over all runs), not the geometric bin center — grazing
    escapes pile toward a bin's high edge, so the center biases the low-μ tail.

Run from the repository root:
    PYTHONPATH=src python3 scripts/tau_sweep.py            # full production run (~40 min)
    PYTHONPATH=src python3 scripts/tau_sweep.py --quick    # plumbing check (~1 min)
"""

import argparse
import os
import threading
import time
from concurrent.futures import ProcessPoolExecutor, as_completed

import numpy as np
import matplotlib

matplotlib.use("Agg")  # headless: this runs unattended in the background
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib import cm  # noqa: E402

from mcrt import (  # noqa: E402
    Simulation,
    chandrasekhar_h,
    eddington_limb_darkening,
    fit_limb_darkening,
)

# --- Production budget (v0.9.7 §3.4): 4×10⁵ escaped photons per (τ, seed). ------
TAU_VALUES = [0.1, 0.3, 1.0, 3.0, 10.0, 30.0]
TARGET_ESCAPED = 400_000
INJECTED_PER_SEED = {           # escaped-matched injected counts (v0.9.7 §3.4 table)
    0.1: 440_000,
    0.3: 500_000,
    1.0: 720_000,
    3.0: 1_300_000,
    10.0: 3_400_000,
    30.0: 9_500_000,
}
N_SEEDS = 5
N_BINS = 20
MU_FLOOR = 0.1                  # μ→0 bins are noisy grazing escapes; excluded from the fit
BASE_SEED = 20260708           # per-(τ, seed) streams spawn off this; distinct from the study base

# Rough s-per-1e6-injected by τ (v0.9.7 measurement) — schedules long tasks first,
# never used for physics. Same numbers as convergence_study_v2.
COST_PER_M = {0.1: 4.2, 0.3: 5.8, 1.0: 12.2, 3.0: 28.3, 10.0: 85.3, 30.0: 252.0}
ESCAPE_FRAC_APPROX = {0.1: 0.916, 0.3: 0.793, 1.0: 0.553, 3.0: 0.302, 10.0: 0.117, 30.0: 0.042}

HEARTBEAT_SEC = 120            # keep the log ticking during the long τ=30 tasks

LIBRARY_PATH = "data/beaming_library.npz"
CURVES_PATH = "data/beaming_tau_curves.png"
SLOPE_PATH = "data/beaming_slope_vs_tau.png"

# Quick-mode budget: ~3000 escaped per row, 2 seeds — a seconds-scale plumbing check.
QUICK_TARGET_ESCAPED = 3_000
QUICK_SEEDS = 2


# --- The parallel sweep ---------------------------------------------------------
def build_tasks(injected_map, n_seeds):
    """Flat (τ, injected_N, seed) task list, longest-cost first for load balancing."""
    tasks = [(tau, injected_map[tau], seed)
             for tau in TAU_VALUES for seed in range(n_seeds)]
    return sorted(tasks, key=lambda t: t[1] * COST_PER_M.get(t[0], 100.0), reverse=True)


def run_one(task):
    """One (τ, N, seed) simulation → per-bin escape histogram. Runs in a worker.

    RNG stream: SeedSequence(BASE_SEED, spawn_key=(τ×10, N, seed)) — keyed on the
    task coordinates, so each run is reproducible regardless of which worker takes
    it or how the grid is later edited (the v0.9.7 seeding convention).

    Returns raw per-bin ``counts`` and ``sum_mu`` (Σμ of escapes in each bin) so the
    pooled mean-escaped-μ grid (fix 8b) can be assembled exactly after all runs land.
    ``b_quick`` is a cheap bin-center slope for the progress log only.
    """
    tau, injected, seed = task
    ss = np.random.SeedSequence(BASE_SEED, spawn_key=(int(round(tau * 10)), int(injected), seed))
    rng = np.random.default_rng(ss)

    t0 = time.perf_counter()
    sim = Simulation(tau_total=tau, num_photons=injected, rng=rng)
    sim.run()
    runtime = time.perf_counter() - t0

    escaped_mu = np.asarray(sim.results["escaped_mu"], dtype=float)
    counts, edges = np.histogram(escaped_mu, bins=N_BINS, range=(0.0, 1.0))
    sum_mu, _ = np.histogram(escaped_mu, bins=N_BINS, range=(0.0, 1.0), weights=escaped_mu)
    centers = 0.5 * (edges[:-1] + edges[1:])

    b_quick = np.nan
    if (counts > 0).all():
        with np.errstate(divide="ignore", invalid="ignore"):
            raw = counts / centers
        b_quick = fit_limb_darkening(centers, raw / raw[-1], mu_floor=MU_FLOOR)[1]

    return {
        "tau": tau, "injected": int(injected), "seed": seed,
        "escaped": int(escaped_mu.size),
        "counts": counts.astype(float), "sum_mu": sum_mu,
        "edges": edges, "runtime": runtime, "b_quick": b_quick,
    }


def run_sweep(tasks, max_workers):
    """Execute all tasks in parallel with a live progress log + heartbeat."""
    records = []
    total = len(tasks)
    t0 = time.perf_counter()
    state = {"done": 0}
    stop = threading.Event()

    def heartbeat():
        while not stop.wait(HEARTBEAT_SEC):
            el = (time.perf_counter() - t0) / 60
            print(f"  [heartbeat] {el:5.1f} min elapsed — "
                  f"{state['done']}/{total} tasks complete", flush=True)

    hb = threading.Thread(target=heartbeat, daemon=True)
    hb.start()
    with ProcessPoolExecutor(max_workers=max_workers) as pool:
        futures = [pool.submit(run_one, t) for t in tasks]
        for fut in as_completed(futures):
            r = fut.result()
            records.append(r)
            state["done"] += 1
            frac = r["escaped"] / r["injected"] * 100
            print(f"  [{state['done']:2d}/{total}] τ={r['tau']:>4} seed={r['seed']} "
                  f"inj={r['injected']:>9,} esc={r['escaped']:>8,} ({frac:4.1f}%) "
                  f"b≈{r['b_quick']:5.3f}  {r['runtime'] / 60:5.1f} min", flush=True)
    stop.set()
    print(f"  sweep wall time: {(time.perf_counter() - t0) / 60:.1f} min", flush=True)
    return records


# --- Pooling: fixes 8a + 8b -----------------------------------------------------
def pooled_mu_grid(records, n_bins):
    """Shared μ grid = mean escaped μ per bin, pooled over every run (fix 8b).

    Photons within a bin are not uniformly distributed (grazing escapes crowd the
    high edge), so the geometric bin center over-states the μ the count represents.
    The pooled mean escaped μ is the low-variance, essentially τ-independent
    estimator of that representative μ; it serves as both the divisor (counts/μ̄)
    and the stored x-location for every τ row. Empty pooled bins → bin center.
    """
    total_counts = np.zeros(n_bins)
    total_sum_mu = np.zeros(n_bins)
    for r in records:
        total_counts += r["counts"]
        total_sum_mu += r["sum_mu"]
    edges = records[0]["edges"]
    centers = 0.5 * (edges[:-1] + edges[1:])
    with np.errstate(divide="ignore", invalid="ignore"):
        mu = np.where(total_counts > 0, total_sum_mu / np.where(total_counts > 0, total_counts, 1.0),
                      centers)
    return mu


def curve_from_counts(counts, mu_grid, mu_floor):
    """One escape histogram → normalized I(μ) curve, plus its fit (a, b).

    Divides by the pooled mean-μ grid (8b), fits I = a(1 + bμ) over μ > mu_floor,
    and normalizes the whole curve by the fitted intercept a (8a). Empty bins →
    NaN (does not occur at production N; kept as a guard).
    """
    with np.errstate(divide="ignore", invalid="ignore"):
        raw = np.where(counts > 0, counts / mu_grid, np.nan)
    fit_mask = np.isfinite(raw) & (mu_grid > mu_floor)
    if fit_mask.sum() < 3:
        return np.full_like(raw, np.nan), np.nan, np.nan
    a, b = fit_limb_darkening(mu_grid[fit_mask], raw[fit_mask], mu_floor=0.0)
    return raw / a, a, b


def _nanstd(values, ddof=1):
    finite = np.asarray(values)[np.isfinite(values)]
    return float(np.std(finite, ddof=ddof)) if len(finite) >= 2 else np.nan


def pool_by_tau(records, mu_grid, n_seeds):
    """Group runs by τ, build per-seed curves, and reduce to pooled ± seed error."""
    n_tau, n_bins = len(TAU_VALUES), N_BINS
    out = {
        "intensity_by_tau": np.full((n_tau, n_bins), np.nan),
        "intensity_std_by_tau": np.full((n_tau, n_bins), np.nan),
        "intensity_by_tau_seed": np.full((n_tau, n_seeds, n_bins), np.nan),
        "b_of_tau": np.full(n_tau, np.nan),
        "b_std_of_tau": np.full(n_tau, np.nan),
        "b_by_tau_seed": np.full((n_tau, n_seeds), np.nan),
        "escaped_by_tau_seed": np.full((n_tau, n_seeds), np.nan),
        "escape_frac_by_tau": np.full(n_tau, np.nan),
    }
    injected = np.full(n_tau, np.nan)  # sourced from the records, so --quick reports honestly

    for ti, tau in enumerate(TAU_VALUES):
        recs = sorted((r for r in records if r["tau"] == tau), key=lambda r: r["seed"])
        injected[ti] = recs[0]["injected"]  # same for every seed of this τ
        for si, r in enumerate(recs):
            inten, _a, b = curve_from_counts(r["counts"], mu_grid, MU_FLOOR)
            out["intensity_by_tau_seed"][ti, si] = inten
            out["b_by_tau_seed"][ti, si] = b
            out["escaped_by_tau_seed"][ti, si] = r["escaped"]
        seed_curves = out["intensity_by_tau_seed"][ti]
        out["intensity_by_tau"][ti] = np.nanmean(seed_curves, axis=0)
        out["intensity_std_by_tau"][ti] = [_nanstd(seed_curves[:, j]) for j in range(n_bins)]
        out["b_of_tau"][ti] = np.nanmean(out["b_by_tau_seed"][ti])
        out["b_std_of_tau"][ti] = _nanstd(out["b_by_tau_seed"][ti])
        out["escape_frac_by_tau"][ti] = np.nanmean(out["escaped_by_tau_seed"][ti]) / injected[ti]
    out["injected_per_seed"] = injected
    return out


def save_library(mu_grid, pooled, n_seeds, path=LIBRARY_PATH):
    provenance = (
        "v0.9.8 production beaming library. Escape-matched injected N targeting "
        f"{TARGET_ESCAPED:,} escaped photons per (tau, seed), {n_seeds} seeds. "
        "Fix 8a: I(mu) normalized by fitted intercept a of I=a(1+b*mu). "
        "Fix 8b: counts divided by pooled mean escaped mu per bin (mu_centers), "
        "not the geometric bin center. Per-(tau,seed) SeedSequence(BASE_SEED, "
        "spawn_key=(tau*10, injected, seed)) streams. See "
        "docs/deep-dives/v0.9.7-convergence-redo.md and docs/next-steps.md Track A."
    )
    np.savez(
        path,
        # --- backward-compatible keys (downstream loaders read these) ---
        tau_values=np.asarray(TAU_VALUES, dtype=float),
        mu_centers=mu_grid,
        intensity_by_tau=pooled["intensity_by_tau"],
        b_of_tau=pooled["b_of_tau"],
        # --- new error / per-seed / provenance keys ---
        intensity_std_by_tau=pooled["intensity_std_by_tau"],
        b_std_of_tau=pooled["b_std_of_tau"],
        escape_frac_by_tau=pooled["escape_frac_by_tau"],
        intensity_by_tau_seed=pooled["intensity_by_tau_seed"],
        b_by_tau_seed=pooled["b_by_tau_seed"],
        escaped_by_tau_seed=pooled["escaped_by_tau_seed"],
        injected_per_seed=pooled["injected_per_seed"],
        n_seeds=n_seeds,
        n_bins=N_BINS,
        mu_floor=MU_FLOOR,
        base_seed=BASE_SEED,
        target_escaped=TARGET_ESCAPED,
        generated_utc=time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        provenance=provenance,
    )
    print(f"✓ Production beaming library saved to {path}")


def print_summary(pooled):
    print("\n=== Production beaming library summary ===")
    print(f"{'τ':>5} {'injected':>10} {'esc mean':>10} {'f_esc':>7} {'b ± σ':>16}")
    print("-" * 52)
    for ti, tau in enumerate(TAU_VALUES):
        esc = np.nanmean(pooled["escaped_by_tau_seed"][ti])
        print(f"{tau:>5} {pooled['injected_per_seed'][ti]:>10,.0f} {esc:>10,.0f} "
              f"{pooled['escape_frac_by_tau'][ti]:>7.3f} "
              f"{pooled['b_of_tau'][ti]:>7.3f} ± {pooled['b_std_of_tau'][ti]:.3f}")
    print("\nSanity (vs v0.9.7 converged): b rises monotonically to a ~1.76 plateau; "
          "τ=10 and τ=30 statistically identical.")


# --- Figures --------------------------------------------------------------------
def plot_curves(mu_grid, pooled, path=CURVES_PATH):
    mu_th = np.linspace(1e-3, 1.0, 200)
    edd = eddington_limb_darkening(mu_th) / eddington_limb_darkening(1.0)
    H = chandrasekhar_h(mu_th)
    H = H / chandrasekhar_h(1.0)[0]

    colors = cm.viridis(np.linspace(0.0, 0.9, len(TAU_VALUES)))
    plt.figure(figsize=(8, 6))
    for tau, intensity, std, color in zip(TAU_VALUES, pooled["intensity_by_tau"],
                                          pooled["intensity_std_by_tau"], colors):
        plt.errorbar(mu_grid, intensity, yerr=std, fmt="o-", color=color, ms=4,
                     capsize=2, lw=1.2, label=f"τ = {tau}")
    plt.plot(mu_th, edd, "r--", lw=1.5, label="Eddington (1 + 1.5μ)")
    plt.plot(mu_th, H, "k-", lw=1.5, label="Chandrasekhar H(μ)")
    plt.xlabel(r"$\mu = \cos(\theta)$")
    plt.ylabel(r"Specific intensity $I(\mu)$ (normalized by fitted intercept $a$)")
    plt.title("Production beaming library: I(μ; τ) with per-bin seed error bars")
    plt.legend(fontsize=9)
    plt.grid(True, alpha=0.3)
    plt.savefig(path)
    plt.close()
    print(f"✓ I(μ; τ) curve family saved to {path}")


def plot_slope_vs_tau(pooled, path=SLOPE_PATH):
    plt.figure(figsize=(8, 5))
    plt.errorbar(TAU_VALUES, pooled["b_of_tau"], yerr=pooled["b_std_of_tau"],
                 fmt="o-", color="#2c7fb8", capsize=3, label="MC best-fit slope b(τ) ± σ")
    plt.axhline(1.5, color="#c0392b", ls="--", label="Eddington (b = 1.5)")
    plt.axhline(0.0, color="#7f8c8d", ls=":", label="Isotropic (b = 0)")
    plt.xscale("log")
    plt.xlabel("Optical depth τ")
    plt.ylabel("Best-fit limb-darkening slope b")
    plt.title("Limb darkening vs. optical depth (monotone rise to a plateau)")
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.savefig(path)
    plt.close()
    print(f"✓ slope-vs-τ plot saved to {path}")


def main():
    parser = argparse.ArgumentParser(description="Production beaming-library sweep (v0.9.8).")
    parser.add_argument("--quick", action="store_true",
                        help="Tiny sweep (~3k escaped/row, 2 seeds) for a plumbing check.")
    parser.add_argument("--workers", type=int, default=max(1, (os.cpu_count() or 4) - 2))
    parser.add_argument("--out", default=LIBRARY_PATH,
                        help="Override the library output path (used by --quick smoke tests).")
    args = parser.parse_args()

    if args.quick:
        injected_map = {t: max(1_500, round(QUICK_TARGET_ESCAPED / ESCAPE_FRAC_APPROX[t]))
                        for t in TAU_VALUES}
        n_seeds = QUICK_SEEDS
        print("QUICK plumbing sweep (not for physics).")
    else:
        injected_map = INJECTED_PER_SEED
        n_seeds = N_SEEDS

    tasks = build_tasks(injected_map, n_seeds)
    total_inj = sum(t[1] for t in tasks)
    cpu_min = sum(t[1] * COST_PER_M.get(t[0], 100.0) / 1e6 for t in tasks) / 60
    print(f"Sweep: {len(tasks)} tasks ({len(TAU_VALUES)} τ × {n_seeds} seeds), "
          f"{total_inj / 1e6:.1f}M photons injected, ~{cpu_min:.0f} CPU-min on "
          f"{args.workers} workers.")
    print(f"Wall time is set by the {n_seeds} longest (τ=30) tasks running in parallel "
          f"(~{injected_map[30.0] * COST_PER_M[30.0] / 1e6 / 60:.0f} min each).", flush=True)

    records = run_sweep(tasks, max_workers=args.workers)
    mu_grid = pooled_mu_grid(records, N_BINS)
    pooled = pool_by_tau(records, mu_grid, n_seeds)

    save_library(mu_grid, pooled, n_seeds, path=args.out)
    print_summary(pooled)
    plot_curves(mu_grid, pooled)
    plot_slope_vs_tau(pooled)
    print("\n✓ Track A2 complete — library + figures written.")


if __name__ == "__main__":
    main()
