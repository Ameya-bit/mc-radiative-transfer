"""Track E3 — is the ΔPF systematic an artifact of the Thomson slab? (robustness rows)

Every real-star number so far uses one limb-darkening law: the Monte Carlo Thomson-
scattering slab, I(μ; τ=10) from ``data/beaming_library.npz``. A referee's cheap attack
is "your slab isn't a hydrogen atmosphere." This driver defuses it by re-running the
*identical* J0740 anchor swap under two independent, analytic beaming laws already in
:mod:`mcrt.theory`:

    Eddington      I(μ) = 1 + (3/2)μ        (the classic grey-atmosphere law)
    Chandrasekhar  I(μ) ∝ H(μ)              (exact isotropic-scattering, ω = 1)

alongside the Thomson slab as the reference row. Nothing else changes — same spots,
compactness, exact Schwarzschild bending (Gate G1), and the same isotropic (I ≡ 1)
baseline — so any law that darkens the limb should reproduce the *sign and geometry-
dependence* of the effect even if the magnitude differs with its slope b. Because PF is a
ratio, ΔPF is invariant to each law's overall normalization; only its μ-shape matters.

Each law is reported at **0 Hz and 346.5 Hz** (v0.9.10 ripple): the diluted + routed spin
numbers are the quotable ones, so the robustness claim must hold there too — both the
residual ΔPF and the peak-normalized waveform-shape RMS. Deliverable: "same sign and
geometry-routing under three independent limb-darkening laws," static and at spin.

Pure interpolation + geometry — seconds.
Run from the repository root:  PYTHONPATH=src:scripts python3 scripts/e3_atmosphere_laws.py
"""

import matplotlib.pyplot as plt
import numpy as np

from mcrt import (
    ExactBending,
    Rotation,
    beaming_lookup,
    chandrasekhar_h,
    eddington_limb_darkening,
    fit_limb_darkening_slope,
    pulsed_fraction,
)
from anchor_lib import SHAPE_TAU, multi_spot_flux, shape_tau_index, waveform_shape_change
import j0740_anchor as j0740

LIBRARY_PATH = "data/beaming_library.npz"
RESULTS_PATH = "data/e3_atmosphere_laws.npz"
FIGURE_PATH = "data/e3_atmosphere_laws.png"

J0740_SPIN_HZ = 346.53
RADIUS_KM = {
    "Riley 2021 (X-PSI ST-U)": 12.39,
    "Miller 2021 (Illinois–Maryland)": 13.713,
}


def build_laws(mu_centers, library_curve):
    """Ordered list of ``(name, beaming_callable, slope_b)`` — Thomson slab + 2 analytic.

    The Thomson slab is the tabulated library row wrapped in the same interpolating
    lookup the anchors use; Eddington and H(μ) are their exact analytic callables, so
    they are evaluated continuously (no clamp) — a genuinely independent limb-darkening
    model, not a re-binning of the same curve. ``slope_b`` is the fitted b of I = a(1+bμ)
    over the library μ-grid, for context on why the magnitudes differ.
    """
    laws = [
        ("Thomson slab (τ=10)", beaming_lookup(mu_centers, library_curve),
         fit_limb_darkening_slope(mu_centers, library_curve)),
        ("Eddington 1+1.5μ", eddington_limb_darkening,
         fit_limb_darkening_slope(mu_centers, eddington_limb_darkening(mu_centers))),
        ("Chandrasekhar H(μ)", chandrasekhar_h,
         fit_limb_darkening_slope(mu_centers, chandrasekhar_h(mu_centers))),
    ]
    return laws


def anchor_metrics(anchor, beaming, bending, rotation):
    """(ΔPF, shape RMS) for one anchor under one beaming law and (bend, rot) pair."""
    iso = multi_spot_flux(anchor.inclination, anchor.compactness, anchor.spots, None,
                          bending=bending, rotation=rotation)
    real = multi_spot_flux(anchor.inclination, anchor.compactness, anchor.spots,
                           beaming, bending=bending, rotation=rotation)
    rms, _ = waveform_shape_change(iso, real)
    return pulsed_fraction(real) - pulsed_fraction(iso), rms


def sweep_laws(laws):
    """Nested dict: results[law_name][anchor_key][tag] = {dpf, rms}, over both spins."""
    results = {}
    for name, beaming, _b in laws:
        results[name] = {}
        for anchor in j0740.ANCHORS:
            key = "riley" if "Riley" in anchor.label else "miller"
            bending = ExactBending(anchor.compactness)
            row = {}
            for tag, rotation in (("0", None),
                                  ("rot", Rotation(spin_hz=J0740_SPIN_HZ,
                                                   radius_km=RADIUS_KM[anchor.label]))):
                dpf, rms = anchor_metrics(anchor, beaming, bending, rotation)
                row[tag] = {"dpf": dpf, "rms": rms}
            results[name][key] = row
    return results


def print_summary(laws, results):
    print("\n" + "=" * 92)
    print(f"E3 — J0740 ΔPF under independent limb-darkening laws (exact bending, δ⁴). "
          f"τ={SHAPE_TAU:g} slab")
    print("=" * 92)
    print(f"\n{'law':<22}{'b':>6}"
          f"{'Riley ΔPF(0)':>14}{'ΔPF(346)':>11}{'RMS(0)':>9}"
          f"{'Miller ΔPF(0)':>15}{'ΔPF(346)':>11}{'RMS(0)':>9}")
    for name, _beaming, b in laws:
        r = results[name]["riley"]
        m = results[name]["miller"]
        print(f"{name:<22}{b:>6.2f}"
              f"{r['0']['dpf']:>+14.4f}{r['rot']['dpf']:>+11.4f}{r['0']['rms']:>9.4f}"
              f"{m['0']['dpf']:>+15.4f}{m['rot']['dpf']:>+11.4f}{m['0']['rms']:>9.4f}")
    print("\nReading: all three laws give a positive, PF-live static ΔPF that the real spin\n"
          "dilutes to +0.02…+0.09 while the waveform-shape RMS stays ~0.10 — same sign, same\n"
          "geometry-routing (PF-visible at J0740, spin-diluted, shape-conserved). The Thomson\n"
          "slab is not special; magnitudes track each law's limb-darkening slope b.")


def plot_result(laws, results, path=FIGURE_PATH):
    """Two panels (Riley, Miller): grouped ΔPF bars per law, 0 Hz vs 346.5 Hz."""
    names = [n for n, _b, _s in laws]
    x = np.arange(len(names))
    w = 0.38
    fig, axes = plt.subplots(1, 2, figsize=(12.5, 4.6), sharey=True,
                             layout="constrained")
    for ax, (key, title) in zip(axes, [("riley", "Riley 2021 (u=0.494)"),
                                        ("miller", "Miller 2021 (u=0.444)")]):
        dpf0 = [results[n][key]["0"]["dpf"] for n in names]
        dpfr = [results[n][key]["rot"]["dpf"] for n in names]
        ax.bar(x - w / 2, dpf0, w, color="#2c7fb8", label="0 Hz (frozen)")
        ax.bar(x + w / 2, dpfr, w, color="#c0392b", label=f"{J0740_SPIN_HZ:g} Hz (δ⁴)")
        for xi, v in zip(x - w / 2, dpf0):
            ax.text(xi, v + 0.004, f"{v:+.3f}", ha="center", fontsize=7.5)
        for xi, v in zip(x + w / 2, dpfr):
            ax.text(xi, v + 0.004, f"{v:+.3f}", ha="center", fontsize=7.5)
        ax.axhline(0.0, color="#7f8c8d", lw=0.8)
        ax.set_xticks(x)
        ax.set_xticklabels(names, rotation=12, ha="right", fontsize=8)
        ax.set_title(title)
        ax.grid(True, axis="y", alpha=0.3)
    axes[0].set_ylabel(r"$\Delta$PF = PF$_{\rm real}$ − PF$_{\rm iso}$")
    axes[0].legend(fontsize=8, loc="upper right")
    fig.suptitle("Beaming systematic is robust across independent limb-darkening laws "
                 "(positive, PF-live, spin-diluted)", fontsize=11)
    fig.savefig(path, dpi=130)
    print(f"✓ figure saved to {path}")


def save_results(laws, results, path=RESULTS_PATH):
    slug = {"Thomson slab (τ=10)": "thomson", "Eddington 1+1.5μ": "eddington",
            "Chandrasekhar H(μ)": "chandra_h"}
    payload = {"tau_shape": np.float64(SHAPE_TAU),
               "j0740_spin_hz": np.float64(J0740_SPIN_HZ),
               "laws": np.array([slug[n] for n, _b, _s in laws])}
    for name, _beaming, b in laws:
        s = slug[name]
        payload[f"{s}_slope_b"] = np.float64(b)
        for key in ("riley", "miller"):
            r = results[name][key]
            payload[f"{s}_{key}_dpf0"] = np.float64(r["0"]["dpf"])
            payload[f"{s}_{key}_dpfrot"] = np.float64(r["rot"]["dpf"])
            payload[f"{s}_{key}_rms0"] = np.float64(r["0"]["rms"])
            payload[f"{s}_{key}_rmsrot"] = np.float64(r["rot"]["rms"])
    np.savez(path, **payload)
    print(f"✓ results saved to {path}")


def main():
    d = np.load(LIBRARY_PATH)
    tau, mu = d["tau_values"], d["mu_centers"]
    ti10 = shape_tau_index(tau)
    library_curve = d["intensity_by_tau"][ti10]

    laws = build_laws(mu, library_curve)
    print("E3 — three independent limb-darkening laws, J0740 anchors, "
          f"exact bending, δ⁴ boost\n  laws: " +
          ", ".join(f"{n} (b={b:.2f})" for n, _bm, b in laws))

    results = sweep_laws(laws)
    print_summary(laws, results)
    save_results(laws, results)
    plot_result(laws, results)


if __name__ == "__main__":
    main()
