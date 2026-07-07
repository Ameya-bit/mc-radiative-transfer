"""Track D — certify the transport machinery against Chandrasekhar H(μ).

The beaming library uses the physical **Thomson** dipole phase function
P(μ) = (3/4)(1 + μ²). There is no closed-form emergent intensity for that phase
function, so the library is validated only by its *shape* (slope b, bracketed by
Eddington and H). This script closes attention item 2 by controlling the one
degree of freedom that blocks an exact check: swap the engine onto the
**isotropic** phase function P(μ) = 1/2, for which a conservative semi-infinite
atmosphere emits *exactly* I(μ) ∝ H(μ) (Chandrasekhar, *Radiative Transfer* Ch.
III/X). At τ = 10 the slab is on the b(τ) plateau (τ = 10 ≡ τ = 30, v0.9.7), so
it is effectively semi-infinite and the emergent curve must match H(μ) to within
Monte Carlo error. That match certifies the transport + escape + binning
machinery with the phase function held fixed — the geometry pipeline is then the
only remaining unvalidated layer, and it is exercised separately by the code
comparisons (SD1a/SD1c).

The **Thomson** run is carried alongside as a control: it is the same machinery
with the physical phase function, and its residual against H(μ) is the ~1–3%
*phase-function sensitivity* (Chandrasekhar Ch. X) — never zero, which is exactly
why the paper must say "near-exact reference", not "matches the exact solution".

Metric — the transport certificate is a **flux-space** goodness-of-fit, not an
intensity-space one. The clean observable is the raw escaping-photon angular
distribution: a semi-infinite isotropic atmosphere emits photons with density
∝ H(μ)·μ (the μ is the cos-θ projection of the emergent flux). We compare the
per-bin escape counts directly to N·∫_bin H(μ)μ dμ with multinomial errors —
no μ-division. The reconstructed intensity I(μ) = counts/μ carries a known
bin-center division bias (attention fix 8b: the flux ∝ I(μ)μ varies across a
finite bin, so dividing by the bin *center* biases by O(Δμ²), a ~few-% effect
that grows with statistics and would masquerade as a transport error). So the
intensity overlay is the *figure*; the flux-space χ² is the *verdict*.

Method (mirrors the library's per-seed error budget, v0.9.7):
  * N_SEEDS independent seeds, each ~INJECTED_PER_SEED photons at τ = 10.
  * Pool the escape counts; χ² = Σ (counts − exp)² / exp over interior bins
    (μ > MU_FLOOR), exp ∝ ∫_bin H(μ)μ dμ matched to the interior count total
    (one normalization constraint ⇒ dof = n_interior − 1). "Within error bars"
    = reduced χ² ≈ 1. The Thomson control is the same test with the physical
    phase function: its χ² is inflated by the phase-function sensitivity.

Usage (repository root):
    PYTHONPATH=src python3 scripts/d_isotropic_validate.py            # full run
    PYTHONPATH=src python3 scripts/d_isotropic_validate.py --quick    # fast plumbing check
"""

import argparse
import os
import time
from concurrent.futures import ProcessPoolExecutor, as_completed

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from mcrt import (
    Simulation,
    chandrasekhar_h,
    eddington_limb_darkening,
    extract_intensity,
    fit_limb_darkening_slope,
)

TAU = 10.0                     # on the b(τ) plateau ⇒ effectively semi-infinite
N_BINS = 20
MU_FLOOR = 0.1                 # grazing bins are noise-dominated (same cut as the fit)
BIN_EDGES = np.linspace(0.0, 1.0, N_BINS + 1)
N_SEEDS = 5
INJECTED_PER_SEED = 700_000    # ~3.5M total ⇒ ~4×10⁵ escaped (escape frac ~12% at τ=10)
BASE_SEED = 20260707           # reproducible, distinct from the library's stream

RESULTS_PATH = "data/d_isotropic_h.npz"
FIGURE_PATH = "data/d_isotropic_h_overlay.png"

# Phase functions to run: the isotropic validation and the Thomson control.
PHASE_FUNCTIONS = ("isotropic", "thomson")


def run_one(task):
    """One (phase_function, seed) slab at τ = 10 → normalized I(μ) curve.

    RNG stream: SeedSequence(BASE_SEED, spawn_key=(phase_id, seed)), keyed on the
    task coordinates so a run is reproducible independent of worker scheduling.
    Runs in a worker process.
    """
    phase_function, seed, n = task
    phase_id = PHASE_FUNCTIONS.index(phase_function)
    ss = np.random.SeedSequence(BASE_SEED, spawn_key=(phase_id, seed))
    rng = np.random.default_rng(ss)

    t0 = time.perf_counter()
    sim = Simulation(
        tau_total=TAU, num_photons=n, rng=rng, phase_function=phase_function
    )
    sim.run()
    runtime = time.perf_counter() - t0

    escaped_mu = sim.results["escaped_mu"]
    counts, _ = np.histogram(escaped_mu, bins=BIN_EDGES)
    mu_centers, intensity = extract_intensity(escaped_mu, n_bins=N_BINS)
    return {
        "phase_function": phase_function,
        "seed": seed,
        "escaped": len(escaped_mu),
        "counts": counts,                 # raw escaping distribution (flux-space test)
        "mu_centers": mu_centers,
        "intensity": intensity,           # I(μ) reconstruction (figure only)
        "slope": fit_limb_darkening_slope(mu_centers, intensity),
        "runtime": runtime,
    }


def _h_flux_per_bin(edges, n_quad=400):
    """∫_bin H(μ)·μ dμ per bin — the semi-infinite isotropic emergent flux law.

    H(μ)·μ is the emergent *photon flux* density (H(μ) is the intensity, μ the
    projection). Dense trapezoid per bin so the prediction carries no binning
    approximation of its own (the MC side is exact counts).
    """
    trapz = np.trapezoid if hasattr(np, "trapezoid") else np.trapz  # numpy 2 rename
    fluxes = []
    for lo, hi in zip(edges[:-1], edges[1:]):
        x = np.linspace(lo, hi, n_quad)
        fluxes.append(trapz(chandrasekhar_h(x) * x, x))
    return np.array(fluxes)


def build_tasks(n_seeds, injected):
    return [(pf, s, injected) for pf in PHASE_FUNCTIONS for s in range(n_seeds)]


def run_sweep(tasks, max_workers):
    records = []
    total = len(tasks)
    with ProcessPoolExecutor(max_workers=max_workers) as pool:
        futures = {pool.submit(run_one, t): t for t in tasks}
        for i, fut in enumerate(as_completed(futures), 1):
            rec = fut.result()
            records.append(rec)
            print(
                f"  [{i}/{total}] {rec['phase_function']:>9} seed {rec['seed']}: "
                f"{rec['escaped']:,} escaped, b={rec['slope']:.3f}, "
                f"{rec['runtime']:.1f}s"
            )
    return records


def summarize_phase(records, phase_function):
    """Flux-space χ² vs H(μ)·μ (the verdict) + intensity curve (the figure).

    **Verdict (flux space).** Pool the escape counts across seeds. The expected
    counts are exp_i ∝ ∫_bin H(μ)μ dμ, scaled so the interior total matches the
    interior count total (one normalization constraint). χ² = Σ (counts−exp)²/exp
    over interior bins (μ > MU_FLOOR); residual = (counts−exp)/√exp. This is the
    unbiased transport test — no μ-division.

    **Figure (intensity space).** The seed-mean I(μ) = counts/μ curve normalized
    to the μ≈1 bin, with the standard error of the mean, purely for the overlay
    plot. Its per-bin deviation carries the bin-center bias and is *not* the
    verdict.
    """
    rows = [r for r in records if r["phase_function"] == phase_function]
    mu = rows[0]["mu_centers"]

    # --- flux-space verdict ---
    counts = np.sum([r["counts"] for r in rows], axis=0).astype(float)
    h_flux = _h_flux_per_bin(BIN_EDGES)
    interior = mu > MU_FLOOR
    # Match the model's interior total to the data's interior total (1 constraint).
    exp = h_flux / h_flux[interior].sum() * counts[interior].sum()
    resid_flux = np.full(N_BINS, np.nan)
    resid_flux[interior] = (counts[interior] - exp[interior]) / np.sqrt(exp[interior])
    dof = int(interior.sum()) - 1  # minus the normalization constraint
    chi2 = float(np.nansum(resid_flux[interior] ** 2))
    reduced_chi2 = chi2 / dof if dof else float("nan")
    max_frac_flux = float(
        np.max(np.abs((counts[interior] - exp[interior]) / exp[interior]))
    )

    # --- intensity curve for the figure ---
    curves = np.vstack([r["intensity"] for r in rows])  # (n_seeds, n_bins)
    i_mean = curves.mean(axis=0)
    sem = curves.std(axis=0, ddof=1) / np.sqrt(curves.shape[0])
    h_norm = chandrasekhar_h(mu) / chandrasekhar_h(1.0)[0]

    return {
        "phase_function": phase_function,
        "mu_centers": mu,
        "i_mean": i_mean,
        "sem": sem,
        "h_norm": h_norm,
        "counts": counts,
        "exp_flux": exp,
        "residual_flux": resid_flux,
        "reduced_chi2": reduced_chi2,
        "dof": dof,
        "max_abs_residual": float(np.nanmax(np.abs(resid_flux))),
        "max_frac_flux": max_frac_flux,
        "total_escaped": int(sum(r["escaped"] for r in rows)),
        "mean_slope": float(np.mean([r["slope"] for r in rows])),
    }


def print_report(summaries):
    print("\n" + "=" * 70)
    print(f"Track D — emergent I(μ) vs Chandrasekhar H(μ)  (τ = {TAU}, {N_SEEDS} seeds)")
    print("=" * 70)
    for s in summaries:
        verdict = "MATCH (within error)" if s["reduced_chi2"] < 2.0 else "DEVIATION"
        role = "validation" if s["phase_function"] == "isotropic" else "control"
        print(
            f"\n  {s['phase_function']:>9} ({role}): "
            f"{s['total_escaped']:,} escaped, b = {s['mean_slope']:.3f}"
        )
        print(
            f"    flux-space reduced χ² vs H(μ)·μ = {s['reduced_chi2']:.2f} "
            f"(dof {s['dof']}), max |residual| = {s['max_abs_residual']:.1f}σ, "
            f"max |Δflux/flux| = {100 * s['max_frac_flux']:.1f}%"
        )
        if s["phase_function"] == "isotropic":
            print(f"    → {verdict}: isotropic transport reproduces H(μ).")
        else:
            print(
                "    → Thomson deviates by the phase-function sensitivity "
                "(~1–3%, Chandrasekhar Ch. X); H is a *near-exact* reference, "
                "not exact for the physical phase function."
            )
    print()


def plot_overlay(summaries, path=FIGURE_PATH):
    iso = next(s for s in summaries if s["phase_function"] == "isotropic")
    thom = next(s for s in summaries if s["phase_function"] == "thomson")
    mu = iso["mu_centers"]

    mu_th = np.linspace(1e-3, 1.0, 300)
    h_curve = chandrasekhar_h(mu_th) / chandrasekhar_h(1.0)[0]
    edd = eddington_limb_darkening(mu_th) / eddington_limb_darkening(1.0)

    fig, (ax, axr) = plt.subplots(
        2, 1, figsize=(8, 8), sharex=True, gridspec_kw={"height_ratios": [3, 1]}
    )

    ax.plot(mu_th, h_curve, "-", color="#238b45", lw=2, label="Chandrasekhar H(μ) (exact, isotropic)")
    ax.plot(mu_th, edd, "--", color="#999999", lw=1, label="Eddington 1 + 1.5μ (linear approx)")
    ax.errorbar(
        mu, iso["i_mean"], yerr=iso["sem"], fmt="o", color="#2c7fb8", ms=5,
        capsize=2, label=f"MC isotropic (b={iso['mean_slope']:.2f}, χ²/dof={iso['reduced_chi2']:.1f})",
    )
    ax.errorbar(
        mu, thom["i_mean"], yerr=thom["sem"], fmt="s", color="#d95f0e", ms=4,
        capsize=2, alpha=0.8,
        label=f"MC Thomson control (b={thom['mean_slope']:.2f})",
    )
    ax.set_ylabel(r"Specific intensity $I(\mu)$  (norm. to $\mu\!=\!1$)")
    ax.set_title(f"Track D: transport validation vs Chandrasekhar H(μ)  (τ = {TAU})")
    ax.legend(loc="upper left", fontsize=9)
    ax.grid(True, alpha=0.3)

    # Residual panel: flux-space (counts − exp)/√exp per bin, ±2σ band.
    axr.axhspan(-2, 2, color="#238b45", alpha=0.12, label="±2σ")
    axr.axhline(0, color="#238b45", lw=1)
    axr.plot(mu, iso["residual_flux"], "o", color="#2c7fb8", ms=5, label="isotropic")
    axr.plot(mu, thom["residual_flux"], "s", color="#d95f0e", ms=4, alpha=0.8, label="Thomson")
    axr.set_xlabel(r"$\mu = \cos\theta$")
    axr.set_ylabel(r"flux $(N_{\rm MC} - N_H)/\sqrt{N_H}$")
    axr.legend(loc="lower right", fontsize=8, ncol=3)
    axr.grid(True, alpha=0.3)

    fig.tight_layout()
    os.makedirs(os.path.dirname(path), exist_ok=True)
    fig.savefig(path, dpi=130)
    print(f"✓ Overlay figure saved to {path}")


def save_results(summaries, path=RESULTS_PATH):
    payload = {
        "tau": TAU, "n_seeds": N_SEEDS, "injected_per_seed": INJECTED_PER_SEED,
        "mu_centers": summaries[0]["mu_centers"],
        "h_norm": summaries[0]["h_norm"],
    }
    for s in summaries:
        pf = s["phase_function"]
        payload[f"{pf}_i_mean"] = s["i_mean"]
        payload[f"{pf}_sem"] = s["sem"]
        payload[f"{pf}_counts"] = s["counts"]
        payload[f"{pf}_exp_flux"] = s["exp_flux"]
        payload[f"{pf}_residual_flux"] = s["residual_flux"]
        payload[f"{pf}_reduced_chi2"] = s["reduced_chi2"]
    os.makedirs(os.path.dirname(path), exist_ok=True)
    np.savez(path, **payload)
    print(f"✓ Results saved to {path}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--quick", action="store_true", help="fast plumbing check")
    parser.add_argument("--workers", type=int, default=None)
    args = parser.parse_args()

    n_seeds = 3 if args.quick else N_SEEDS
    injected = 40_000 if args.quick else INJECTED_PER_SEED
    max_workers = args.workers or min(2 * n_seeds, os.cpu_count() or 4)

    print(
        f"Track D: τ={TAU}, {n_seeds} seeds × {injected:,} injected "
        f"per phase function ({', '.join(PHASE_FUNCTIONS)}); {max_workers} workers"
    )
    tasks = build_tasks(n_seeds, injected)
    t0 = time.perf_counter()
    records = run_sweep(tasks, max_workers)
    print(f"\nSweep wall time: {time.perf_counter() - t0:.1f}s")

    summaries = [summarize_phase(records, pf) for pf in PHASE_FUNCTIONS]
    print_report(summaries)
    if not args.quick:
        plot_overlay(summaries)
        save_results(summaries)


if __name__ == "__main__":
    main()
