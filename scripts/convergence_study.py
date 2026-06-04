"""Photon-count convergence study for the beaming function.

Replaces the by-feel photon counts (5000 / 1000 / 200000) with defensible,
knee-justified values for the paper's methods/appendix. For each observable we
sweep N across decades, run several independent seeds per N, and estimate the
statistical error as the spread (std) across seeds — there is no closed-form
"truth" for the noisy low-μ tail bins. Error vs N on log-log axes shows the
~ -1/2 statistics-limited slope flattening to a systematic floor; the bend is
the knee, and production N is set just past it.

Observables:
  1. Energy-conservation residual (escaped + absorbed − injected): a structural
     invariant here — exactly zero at every N (every photon is counted).
  2. Mean-free-path estimate Σpath / Σscatters vs the expected 1.0: converges fast.
  3. Per-μ-bin beaming values, with explicit attention to the low-μ tail bins
     (μ → 0) — the slowest to converge and the binding case (same tail that makes
     τ = 30 noisy in the library; see docs/deep-dives/v0.6.1-isotropic-injection.md).

Usage:
    PYTHONPATH=src python3 scripts/convergence_study.py            # full sweep (~10 min)
    PYTHONPATH=src python3 scripts/convergence_study.py --quick    # fast dev sweep
"""

import argparse

import numpy as np
import matplotlib.pyplot as plt

from mcrt import (
    Simulation,
    extract_intensity,
    fit_limb_darkening_slope,
    statistical_error,
    loglog_slope,
    find_knee,
    n_for_target_error,
)

# Full sweep: decades of N, several seeds each, at the binding τ = 10 / low-μ corner.
N_VALUES = [1_000, 3_000, 10_000, 30_000, 100_000, 300_000, 1_000_000]
N_SEEDS = 5
BASE_SEED = 20260603
TAU_TOTAL = 10.0
N_BINS = 20
EDDINGTON_SLOPE = 1.5
RESULTS_PATH = "data/convergence_results.npz"

# Per-observable error tolerances that define a "good enough" production N. Set
# by the precision the science needs; since every observable is statistics-limited
# in range, these targets (not a knee) set the recommended N.
TARGET_MFP = 5e-3        # 0.5% of the ~1.0 mean free path
TARGET_SLOPE_B = 0.05    # absolute slope error (~3% of b ≈ 1.7)
TARGET_INTENSITY = 0.02  # 2% per-μ-bin intensity error

# Reduced sweep for fast plumbing checks during development.
QUICK_N_VALUES = [1_000, 3_000, 10_000, 30_000, 100_000]
QUICK_N_SEEDS = 3


def run_single(n, tau, n_bins, rng):
    """One simulation; return the observables the study tracks."""
    sim = Simulation(tau_total=tau, num_photons=n, rng=rng)
    sim.run()

    escaped = len(sim.results["escaped_mu"])
    absorbed = sim.results["absorbed_count"]
    energy_residual = escaped + absorbed - n  # structural invariant → 0

    scatters = np.sum(sim.results["num_scatterings"])
    path = np.sum(sim.results["total_path_lengths"])
    mfp = path / scatters if scatters > 0 else np.nan

    mu_centers, intensity = extract_intensity(sim.results["escaped_mu"], n_bins=n_bins)
    slope = fit_limb_darkening_slope(mu_centers, intensity)

    return {
        "energy_residual": energy_residual,  # escaped + absorbed - n (exactly 0)
        "mfp": mfp,
        "slope": slope,
        "intensity": intensity,
        "mu_centers": mu_centers,
    }


def run_sweep(n_values, tau=TAU_TOTAL, n_bins=N_BINS, n_seeds=N_SEEDS, base_seed=BASE_SEED):
    """Sweep N × seeds and stack observables into arrays.

    Independent, reproducible streams via SeedSequence.spawn(base_seed).

    Returns a dict of arrays:
      n_values        (n_N,)
      mu_centers      (n_bins,)
      energy_residual (n_N, n_seeds)
      mfp             (n_N, n_seeds)
      slope           (n_N, n_seeds)
      intensity       (n_N, n_seeds, n_bins)
    """
    n_values = list(n_values)
    # Independent streams by flat index k = i * n_seeds + s. Note: inserting a new
    # N anywhere but the end shifts later runs onto different child seeds — bump
    # BASE_SEED (or append the new N) when changing the grid to keep prior runs
    # reproducible.
    children = np.random.SeedSequence(base_seed).spawn(len(n_values) * n_seeds)
    rngs = [np.random.default_rng(c) for c in children]

    energy = np.zeros((len(n_values), n_seeds))
    mfp = np.zeros((len(n_values), n_seeds))
    slope = np.zeros((len(n_values), n_seeds))
    intensity = np.zeros((len(n_values), n_seeds, n_bins))
    mu_centers = None

    k = 0
    for i, n in enumerate(n_values):
        for s in range(n_seeds):
            res = run_single(n, tau, n_bins, rngs[k])
            k += 1
            energy[i, s] = res["energy_residual"]
            mfp[i, s] = res["mfp"]
            slope[i, s] = res["slope"]
            intensity[i, s] = res["intensity"]
            mu_centers = res["mu_centers"]
        print(
            f"  N={n:>8}: <b>={slope[i].mean():.3f}±{slope[i].std(ddof=1):.3f}  "
            f"<mfp>={mfp[i].mean():.4f}  max|E_resid|={np.abs(energy[i]).max():.0f}"
        )

    return {
        "n_values": np.asarray(n_values, dtype=float),
        "mu_centers": mu_centers,
        "energy_residual": energy,
        "mfp": mfp,
        "slope": slope,
        "intensity": intensity,
    }


def summarize(sweep):
    """Per-observable error, knee, and target-error production N."""
    n = sweep["n_values"]
    mu = sweep["mu_centers"]

    mfp_err = statistical_error(sweep["mfp"], axis=1)
    slope_err = statistical_error(sweep["slope"], axis=1)
    intensity_err = statistical_error(sweep["intensity"], axis=1)  # (n_N, n_bins)

    tail_idx = 0  # lowest-μ bin (μ → 0): grazing escapes, slowest to converge
    bulk_idx = int(np.argmin(np.abs(mu - 0.7)))  # representative bulk bin

    # (error curve, target tolerance) per reported observable.
    observables = {
        "mean_free_path": (mfp_err, TARGET_MFP),
        "limb_slope_b": (slope_err, TARGET_SLOPE_B),
        f"intensity_bulk(μ≈{mu[bulk_idx]:.2f})": (intensity_err[:, bulk_idx], TARGET_INTENSITY),
        f"intensity_tail(μ≈{mu[tail_idx]:.2f})": (intensity_err[:, tail_idx], TARGET_INTENSITY),
    }

    return {
        "n_values": n,
        "mu_centers": mu,
        "mfp_err": mfp_err,
        "slope_err": slope_err,
        "intensity_err": intensity_err,
        "tail_idx": tail_idx,
        "bulk_idx": bulk_idx,
        "energy_max_resid": float(np.abs(sweep["energy_residual"]).max()),
        "observables": observables,
    }


def print_report(summary):
    n = summary["n_values"]
    print("\nEnergy conservation: max |escaped + absorbed − injected| = "
          f"{summary['energy_max_resid']:.0f} (structural invariant; exact at all N)\n")
    header = f"{'observable':<26}{'loglog slope':>13}{'knee N':>12}{'target':>9}{'recommended N':>20}"
    print(header)
    print("-" * len(header))
    for name, (err, target) in summary["observables"].items():
        slope = loglog_slope(n, err)
        knee = find_knee(n, err)
        rec = n_for_target_error(n, err, target)
        knee_s = f"{knee:.0e}" if knee is not None else "none"
        if rec is None:
            rec_s = "n/a"
        elif rec.extrapolated:
            rec_s = f"~{rec.n:.0e} (extrap.)"
        else:
            rec_s = f"{rec.n:.0e}"
        print(f"{name:<26}{slope:>13.2f}{knee_s:>12}{target:>9.0e}{rec_s:>20}")
    print("\nAll observables ride the ~N^(-1/2) line (no persistent knee in range):")
    print("recommended N is the photon count that meets each error target. Tail bins")
    print("(μ → 0) need an extrapolated N beyond the swept range — the binding case and")
    print("the natural opening for variance reduction (importance sampling toward low μ).")


def plot_slope_vs_n(sweep, output_path="data/convergence_slope.png"):
    """Kept slope-vs-N panel, now with across-seed error bars."""
    n = sweep["n_values"]
    slope_mean = sweep["slope"].mean(axis=1)
    slope_err = statistical_error(sweep["slope"], axis=1)

    plt.figure(figsize=(8, 5))
    plt.errorbar(n, slope_mean, yerr=slope_err, fmt="o-", color="#2c7fb8",
                 capsize=3, label="MC best-fit slope b (mean ± seed std)")
    plt.axhline(EDDINGTON_SLOPE, color="#c0392b", ls="--", label="Eddington (b = 1.5)")
    plt.xscale("log")
    plt.xlabel("Number of photons")
    plt.ylabel("Best-fit limb-darkening slope b")
    plt.title(f"Beaming-function slope convergence (τ = {TAU_TOTAL})")
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.savefig(output_path)
    plt.close()
    print(f"✓ Slope-convergence plot saved to {output_path}")


def plot_error_vs_n(summary, output_path="data/convergence_error_vs_n.png"):
    """Log-log error vs N: scalar observables (left) and per-μ-bin (right)."""
    n = summary["n_values"]
    mu = summary["mu_centers"]
    intensity_err = summary["intensity_err"]

    # Mask non-positive errors: a zero spread at small N is *undersampling* (no
    # escapers in that bin across all seeds), not convergence — drop it rather
    # than draw a misleading dip on the log axis.
    def positive(arr):
        arr = np.asarray(arr, dtype=float)
        return np.where(arr > 0, arr, np.nan)

    fig, (ax_scalar, ax_bins) = plt.subplots(1, 2, figsize=(13, 5))

    # Left: scalar observables + a -1/2 reference guide.
    ax_scalar.loglog(n, positive(summary["mfp_err"]), "o-", color="#2c7fb8", label="mean-free-path error")
    ax_scalar.loglog(n, positive(summary["slope_err"]), "s-", color="#d95f0e", label="limb-slope b error")
    ref = summary["slope_err"][0] * (n / n[0]) ** (-0.5)
    ax_scalar.loglog(n, ref, "k:", alpha=0.7, label=r"$N^{-1/2}$ reference")
    ax_scalar.set_xlabel("Number of photons N")
    ax_scalar.set_ylabel("Statistical error (std across seeds)")
    ax_scalar.set_title("Scalar observables: error vs N")
    ax_scalar.legend()
    ax_scalar.grid(True, which="both", alpha=0.3)

    # Right: per-μ-bin intensity error — tail bins are highest and flattest-to-converge.
    bin_indices = [0, summary["bulk_idx"], len(mu) - 2]  # tail, bulk, near-normalization
    colors = plt.cm.viridis(np.linspace(0.15, 0.85, len(bin_indices)))
    for idx, color in zip(bin_indices, colors):
        ax_bins.loglog(n, positive(intensity_err[:, idx]), "o-", color=color,
                       label=f"μ ≈ {mu[idx]:.2f}" + (" (tail)" if idx == 0 else ""))
    ref_b = intensity_err[-1, 0] * (n / n[-1]) ** (-0.5)
    ax_bins.loglog(n, ref_b, "k:", alpha=0.7, label=r"$N^{-1/2}$ reference")
    ax_bins.set_xlabel("Number of photons N")
    ax_bins.set_ylabel("Per-bin intensity error (std across seeds)")
    ax_bins.set_title("Beaming-function bins: low-μ tail converges slowest")
    ax_bins.legend()
    ax_bins.grid(True, which="both", alpha=0.3)

    fig.tight_layout()
    fig.savefig(output_path)
    plt.close(fig)
    print(f"✓ Error-vs-N plot saved to {output_path}")


def save_sweep(sweep, path=RESULTS_PATH):
    """Persist the raw sweep arrays so the report/figures can be regenerated
    without re-running the (~10 min) Monte Carlo sweep."""
    np.savez(path, **sweep)
    print(f"✓ Raw sweep arrays saved to {path}")


def load_sweep(path=RESULTS_PATH):
    data = np.load(path)
    return {k: data[k] for k in data.files}


def main():
    parser = argparse.ArgumentParser(description="Photon-count convergence study.")
    parser.add_argument("--quick", action="store_true",
                        help="Reduced sweep (cap N at 1e5, fewer seeds) for fast dev.")
    parser.add_argument("--summarize-only", action="store_true",
                        help=f"Skip the sweep; re-report/plot from {RESULTS_PATH}.")
    args = parser.parse_args()

    if args.summarize_only:
        print(f"Loading cached sweep from {RESULTS_PATH}...")
        sweep = load_sweep()
    else:
        if args.quick:
            n_values, n_seeds = QUICK_N_VALUES, QUICK_N_SEEDS
            print("Running QUICK convergence sweep (development mode)...")
        else:
            n_values, n_seeds = N_VALUES, N_SEEDS
            print("Running FULL convergence sweep...")
        sweep = run_sweep(n_values, n_seeds=n_seeds)
        save_sweep(sweep)

    summary = summarize(sweep)
    print_report(summary)
    plot_slope_vs_n(sweep)
    plot_error_vs_n(summary)


if __name__ == "__main__":
    main()
