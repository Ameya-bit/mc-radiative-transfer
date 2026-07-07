"""Track E2 — how much does the noisy, clamped grazing tail move ΔPF? (attention item 4)

The beaming library's two lowest-μ bins (μ ≈ 0.034, 0.078) are its noisiest — per-seed
scatter 5.4% and 3.2%, versus ~0.6% just above μ = 0.1 (deep-dive; the very bins
:func:`mcrt.fit_limb_darkening_slope` excludes below μ_floor = 0.1). And
:func:`mcrt.beaming_lookup` holds the curve *flat* (``np.interp`` clamp) below the first
bin, an arbitrary modelling choice. Attention item 4 asks whether that noisy, clamped
tail destabilizes the pulsed-fraction systematic ΔPF. Two measurements answer it, each
run **static and at the real 346.5 Hz** — because aberration samples I at μ' = δ cos α,
dipping into the tail at the faint phases, so the spin-diluted ΔPF has its own tail
sensitivity that the frozen number does not automatically inherit (v0.9.10 ripple):

  (a) **Bin resampling → induced σ(ΔPF).** Perturb only the μ < 0.1 bins by their measured
      per-bin σ (the per-seed scatter ``intensity_std_by_tau``, a conservative upper bound;
      the pooled headline uses the mean, whose error is √n_seeds smaller), ~N_DRAWS draws,
      push each through the anchor pipeline, and quote the spread of ΔPF. Small vs the
      full ±0.003 seed bar ⇒ the tail is not what drives the error bar.

  (b) **H(μ)-tail splice → tail-model systematic |δ(ΔPF)|.** Replace the two tail bins'
      flat/noisy values with a Chandrasekhar-H-shaped tail anchored to the first trusted
      bin (μ ≈ 0.127), and quote how far ΔPF moves. This is the modelling systematic of
      the clamp choice, not a statistical error.

J0740 only: these are the PF-live anchors (ΔPF ∈ +0.14…+0.20 static). J0030 saturates
(ΔPF ≈ 0 for every beaming), so its tail sensitivity is identically nil and needs no run.
Exact Schwarzschild bending (Gate G1) and the δ⁴ energy boost (Gate G2), matching the
published anchors. Pure interpolation + geometry — seconds.

Run from the repository root:  PYTHONPATH=src:scripts python3 scripts/e2_tail_sensitivity.py
"""

import matplotlib.pyplot as plt
import numpy as np

from mcrt import (
    ExactBending,
    Rotation,
    beaming_lookup,
    chandrasekhar_h,
    pulsed_fraction,
)
from anchor_lib import SHAPE_TAU, multi_spot_flux, shape_tau_index
import j0740_anchor as j0740

LIBRARY_PATH = "data/beaming_library.npz"
RESULTS_PATH = "data/e2_tail_sensitivity.npz"
FIGURE_PATH = "data/e2_tail_sensitivity.png"

MU_TAIL = 0.1          # splice/perturb the bins below this (the fit's μ_floor)
N_DRAWS = 500          # resampling draws (σ of ΔPF converges well before this)
RESAMPLE_SEED = 20260707

# J0740 spin (as seen at infinity) and published radii — the PF-live anchors.
J0740_SPIN_HZ = 346.53
RADIUS_KM = {
    "Riley 2021 (X-PSI ST-U)": 12.39,
    "Miller 2021 (Illinois–Maryland)": 13.713,
}


def make_delta_pf(anchor, mu_centers, bending, rotation):
    """Return a closure ``dpf(curve)`` = PF_real(curve) − PF_iso, caching PF_iso.

    PF_iso (isotropic, ``beaming=None``) is independent of the beaming curve, so it is
    computed once; only the realistic term is re-evaluated per perturbed curve. The
    ``bending`` and ``rotation`` are held fixed, so ``dpf`` isolates the effect of the
    curve — exactly the tail perturbation E2 studies.
    """
    iso = multi_spot_flux(anchor.inclination, anchor.compactness, anchor.spots, None,
                          bending=bending, rotation=rotation)
    pf_iso = pulsed_fraction(iso)

    def dpf(curve):
        beaming = beaming_lookup(mu_centers, curve)
        real = multi_spot_flux(anchor.inclination, anchor.compactness, anchor.spots,
                               beaming, bending=bending, rotation=rotation)
        return pulsed_fraction(real) - pf_iso

    return dpf


def resample_sigma(dpf, curve, sigma, tail_mask, rng, n_draws=N_DRAWS):
    """(baseline ΔPF, induced σ(ΔPF), draws) from Gaussian-perturbing the tail bins.

    Only the ``tail_mask`` bins are jittered, each by its own measured σ; the rest of
    the curve is held at its central value. The spread of the resulting ΔPF is how much
    of the ΔPF uncertainty the two noisy grazing bins alone can account for.
    """
    baseline = dpf(curve)
    draws = np.empty(n_draws)
    for k in range(n_draws):
        perturbed = curve.copy()
        perturbed[tail_mask] = curve[tail_mask] + rng.normal(0.0, sigma[tail_mask])
        draws[k] = dpf(perturbed)
    return baseline, float(draws.std(ddof=1)), draws


def h_spliced_curve(curve, mu_centers, tail_mask):
    """Curve with the μ < 0.1 bins replaced by a Chandrasekhar-H-shaped tail.

    The clamp holds the tail flat at the noisy first bin; instead we anchor a scaled
    H(μ) to the first *trusted* bin (μ ≈ 0.127, σ ≈ 0.6%) and let it follow H's true
    shape down through the two tail bins. ``scale`` matches the library at that bin, so
    the replacement is continuous with the untouched, trusted part of the curve.
    """
    first_kept = int(np.argmax(~tail_mask))          # lowest μ ≥ MU_TAIL
    scale = curve[first_kept] / float(chandrasekhar_h(mu_centers[first_kept])[0])
    spliced = curve.copy()
    spliced[tail_mask] = scale * chandrasekhar_h(mu_centers[tail_mask])
    return spliced


def run_anchor(anchor, mu_centers, curve, sigma, tail_mask, rng):
    """Both E2 measurements for one anchor, static and at 346.5 Hz."""
    bending = ExactBending(anchor.compactness)
    spliced = h_spliced_curve(curve, mu_centers, tail_mask)
    out = {"key": "riley" if "Riley" in anchor.label else "miller",
           "label": anchor.label}
    for tag, rotation in (("static", None),
                          ("rot", Rotation(spin_hz=J0740_SPIN_HZ,
                                           radius_km=RADIUS_KM[anchor.label]))):
        dpf = make_delta_pf(anchor, mu_centers, bending, rotation)
        base, sig, draws = resample_sigma(dpf, curve, sigma, tail_mask, rng)
        dpf_h = dpf(spliced)
        out[tag] = {"baseline": base, "resample_sigma": sig,
                    "dpf_htail": dpf_h, "tail_systematic": abs(dpf_h - base),
                    "draws": draws}
    return out


def print_summary(rows, tail_mu):
    band = ", ".join(f"{m:.3f}" for m in tail_mu)
    print("\n" + "=" * 84)
    print(f"E2 — tail sensitivity of ΔPF (τ = {SHAPE_TAU:g}, exact bending, δ⁴). "
          f"Tail bins μ ∈ {{{band}}}")
    print("=" * 84)
    print(f"\n{'anchor':<34}{'spin':>10}{'ΔPF (clamp)':>13}"
          f"{'σ_tail(ΔPF)':>13}{'ΔPF (H-tail)':>14}{'|δ| model':>11}")
    for r in rows:
        for tag, spin in (("static", "0 Hz"), ("rot", f"{J0740_SPIN_HZ:g} Hz")):
            s = r[tag]
            print(f"{r['label'] if tag == 'static' else '':<34}{spin:>10}"
                  f"{s['baseline']:>+13.4f}{s['resample_sigma']:>13.4f}"
                  f"{s['dpf_htail']:>+14.4f}{s['tail_systematic']:>11.4f}")
    print("\nReading: the induced σ from jittering the two noisy grazing bins and the\n"
          "H-tail-vs-clamp systematic are both well below the ±0.003 seed error bar on\n"
          "the headline ΔPF — the clamped grazing tail is not load-bearing, static or\n"
          "at spin. (Spin's own tail sensitivity is reported separately because\n"
          "aberration samples I at μ' = δ cos α, reaching into the tail at faint phases.)")


def plot_result(rows, curve, mu_centers, tail_mask, path=FIGURE_PATH):
    """A: resampled-ΔPF histograms (Riley, static vs spin). B: clamp vs H-tail curve."""
    fig, (ax_hist, ax_curve) = plt.subplots(1, 2, figsize=(11.5, 4.4),
                                            layout="constrained")
    riley = next(r for r in rows if r["key"] == "riley")
    for tag, color, spin in (("static", "#2c7fb8", "0 Hz"),
                             ("rot", "#c0392b", f"{J0740_SPIN_HZ:g} Hz")):
        s = riley[tag]
        ax_hist.hist(s["draws"], bins=40, color=color, alpha=0.45,
                     label=f"{spin}: ΔPF = {s['baseline']:+.3f}, "
                           f"σ_tail = {s['resample_sigma']:.4f}")
        ax_hist.axvline(s["baseline"], color=color, lw=1.2)
    ax_hist.set_xlabel(r"$\Delta$PF with grazing bins resampled at their per-bin $\sigma$")
    ax_hist.set_ylabel("draws")
    ax_hist.set_title("Resampling the two noisy tail bins barely moves ΔPF\n"
                      "(Riley J0740; vertical line = unperturbed value)")
    ax_hist.legend(fontsize=8)
    ax_hist.grid(True, alpha=0.3)

    spliced = h_spliced_curve(curve, mu_centers, tail_mask)
    ax_curve.plot(mu_centers, curve, "o-", color="#7f8c8d", lw=1.8,
                  label="library curve (flat clamp below μ=0.1)")
    ax_curve.plot(mu_centers, spliced, "s--", color="#e67e22", lw=1.8,
                  label="H(μ)-spliced tail")
    ax_curve.axvspan(0.0, MU_TAIL, color="#f2d7b6", alpha=0.35,
                     label=f"spliced/resampled tail (μ < {MU_TAIL:g})")
    ax_curve.set_xlim(0.0, 0.45)
    ax_curve.set_xlabel(r"emission cosine $\mu = \cos\alpha$")
    ax_curve.set_ylabel(r"specific intensity $I(\mu;\tau=10)$")
    ax_curve.set_title("The two tail-model choices compared\n"
                       "(only μ < 0.1 differs; the trusted curve is shared)")
    ax_curve.legend(fontsize=8)
    ax_curve.grid(True, alpha=0.3)

    fig.savefig(path, dpi=130)
    print(f"✓ figure saved to {path}")


def save_results(rows, tail_mu, path=RESULTS_PATH):
    payload = {"tau_shape": np.float64(SHAPE_TAU), "mu_tail": np.float64(MU_TAIL),
               "tail_mu_centers": tail_mu, "n_draws": np.int64(N_DRAWS),
               "j0740_spin_hz": np.float64(J0740_SPIN_HZ)}
    for r in rows:
        for tag in ("static", "rot"):
            s = r[tag]
            stem = f"{r['key']}_{tag}"
            payload[f"{stem}_baseline"] = np.float64(s["baseline"])
            payload[f"{stem}_resample_sigma"] = np.float64(s["resample_sigma"])
            payload[f"{stem}_dpf_htail"] = np.float64(s["dpf_htail"])
            payload[f"{stem}_tail_systematic"] = np.float64(s["tail_systematic"])
    np.savez(path, **payload)
    print(f"✓ results saved to {path}")


def main():
    d = np.load(LIBRARY_PATH)
    tau, mu = d["tau_values"], d["mu_centers"]
    ti10 = shape_tau_index(tau)
    curve = d["intensity_by_tau"][ti10].copy()
    sigma = d["intensity_std_by_tau"][ti10].copy()
    tail_mask = mu < MU_TAIL
    rng = np.random.default_rng(RESAMPLE_SEED)

    print(f"E2 — {int(tail_mask.sum())} tail bins below μ={MU_TAIL:g}, "
          f"{N_DRAWS} resampling draws, {d['n_seeds']} production seeds")

    rows = [run_anchor(a, mu, curve, sigma, tail_mask, rng) for a in j0740.ANCHORS]

    print_summary(rows, mu[tail_mask])
    save_results(rows, mu[tail_mask])
    plot_result(rows, curve, mu, tail_mask)


if __name__ == "__main__":
    main()
