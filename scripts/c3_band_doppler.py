"""Track C3 — the band-limited Doppler dilution of ΔPF (closing Gate G2's number).

C2 measured the Doppler coupling *bolometrically* and found the headline dilution
(+0.137 → +0.037 Riley under δ⁴). That left two open ambiguities, both now
headline-sized relative to the diluted number:

  1. **δ⁴ energy vs δ³ photon** — a ~0.017 spread from the flux convention alone.
  2. **Bolometric vs band-limited** — NICER counts photons in a fixed calibrated
     band that sits on the **Wien tail** of these soft spots (E/kT_obs ≈ 5–25),
     where the flux responds to the Doppler factor *exponentially*, i.e. steeper
     than any δⁿ.

C3 replaces the bolometric δⁿ with the exact in-band blackbody weight
(:func:`mcrt.rotating.band_boost`): for a comoving Planck spectrum the δ³ intensity
boost cancels the E′³ Planck numerator (the same collapse the SD1c validation uses)
and *all* Doppler/redshift dependence enters through one effective observed
temperature T_obs = δ√(1−u)·kT. In-band photon counts (the NICER-faithful
convention) then weight the pulse by Φ₂(δ)/Φ₂(1) with

    Φ_k(δ) = ∫_{E1}^{E2} E^k dE / (exp[E/(δ√(1−u)kT)] − 1),   k = 2 counts / 3 energy,

whose wide-band limits are exactly δ³/δ⁴ — so this measurement *contains* C2's
conventions as limits and interpolates the honest answer between/beyond them.

Inputs (published values, see docs/paper/references.md):
  - Riley 2021 ST-U:  log₁₀T = 5.988 (p) / 5.992 (s) → kT ≈ 0.0842 keV (equal-
    temperature caps to ~1%, so one kT keeps the multi-spot weights exact);
    NICER channels [30, 150) → 0.3–1.5 keV.
  - Miller 2021:  kT₁ = kT₂ = 0.094 keV; NICER channels 30–123 → 0.3–1.24 keV.

**Blackbody stand-in caveat:** the published fits use hydrogen atmospheres, not
blackbodies, and our angular I(μ; τ) stays the grey library (separability
I′(E′,μ′) = Planck × I(μ′)). C3 is therefore the *spectral-leverage* measurement —
how much the band's Wien-tail sensitivity moves the dilution — not a NICER response
simulation. kT and band-edge robustness scans quantify how much that choice matters.

Structure mirrors C2: per-seed difference-of-differences on the production 5-seed
library, exact bending, iso and real sharing the same rotation treatment. Also
re-checks **shape routing in-band** (the systematic must still live in the waveform).
Deterministic, seconds to run.

Run from the repository root:  PYTHONPATH=src:scripts python3 scripts/c3_band_doppler.py
"""

import numpy as np

from mcrt import BandSpectrum, ExactBending, Rotation, band_boost, beaming_lookup, \
    pulsed_fraction, spot_speed
from anchor_lib import SHAPE_TAU, multi_spot_flux, shape_tau_index, waveform_shape_change
from c2_doppler_coupling import J0740_SPIN_HZ, RADIUS_KM, harmonics, per_seed_delta_pf
import j0740_anchor as j0740

LIBRARY_PATH = "data/beaming_library.npz"
RESULTS_PATH = "data/c3_band_doppler.npz"

RILEY = "Riley 2021 (X-PSI ST-U)"
MILLER = "Miller 2021 (Illinois–Maryland)"

# Comoving effective temperatures (blackbody stand-in for the H-atmosphere fits):
# k_B·10^5.99 K = 0.0842 keV (Riley, both caps); Miller Table 7: kT = 0.094 keV.
KT_KEV = {RILEY: 0.0842, MILLER: 0.094}
# Calibrated NICER bands actually fit: Riley channels [30,150), Miller channels 30–123.
BAND_KEV = {RILEY: (0.3, 1.5), MILLER: (0.3, 1.24)}

# Robustness grids (Riley): band edges and kT scanned around the adopted values.
E_MIN_SCAN = (0.25, 0.3, 0.4)
E_MAX_SCAN = (1.0, 1.5, 2.0)
KT_SCAN = (0.07, 0.0842, 0.10)


def anchor_band(label: str) -> BandSpectrum:
    """The adopted BandSpectrum for a published anchor."""
    e_min, e_max = BAND_KEV[label]
    return BandSpectrum(kt_kev=KT_KEV[label], e_min_kev=e_min, e_max_kev=e_max)


def variant_rotations(label: str):
    """The (tag, Rotation) ladder C3 compares, mildest to steepest boost."""
    r_km = RADIUS_KM[label]
    return [
        ("0hz", None),
        ("d3", Rotation(J0740_SPIN_HZ, r_km, photon_flux=True)),
        ("d4", Rotation(J0740_SPIN_HZ, r_km)),
        ("band_en", Rotation(J0740_SPIN_HZ, r_km, band=anchor_band(label))),
        ("band_ph", Rotation(J0740_SPIN_HZ, r_km, photon_flux=True,
                             band=anchor_band(label))),
    ]


VARIANT_NAMES = {
    "0hz": "0 Hz (frozen)",
    "d3": "δ³ bolometric photon",
    "d4": "δ⁴ bolometric energy",
    "band_en": "band-limited energy",
    "band_ph": "band-limited photon (NICER-faithful)",
}


def effective_exponent(label: str, compactness: float, beta_eq: float) -> float:
    """d ln Φ₂/d ln δ at the anchor's blueshift extreme — the 'how steep' diagnostic.

    4 would mean the band behaves like the bolometric δ⁴; the Wien tail pushes it
    well above (deep-dive §8.2 direction argument, now a number).
    """
    delta = 1.0 + beta_eq                     # ~the maximum blueshift, δ ≈ 1 + β
    return float(np.log(band_boost(delta, compactness, anchor_band(label)))
                 / np.log(delta))


def coupling_table(seed_curves_at_tau, mu_centers):
    """ΔPF ± σ for every variant × anchor; returns keyed rows for the payload."""
    rows = {}
    print("=== C3: ΔPF(τ=10) per Doppler treatment, 5 production seeds, exact bending ===")
    for anchor in j0740.ANCHORS:
        bend = ExactBending(anchor.compactness)
        beta_eq = spot_speed(J0740_SPIN_HZ, RADIUS_KM[anchor.label], np.pi / 2.0,
                             anchor.compactness)
        n_eff = effective_exponent(anchor.label, anchor.compactness, beta_eq)
        key = "riley" if "Riley" in anchor.label else "miller"
        kt, (e1, e2) = KT_KEV[anchor.label], BAND_KEV[anchor.label]
        print(f"\n{anchor.label}: kT = {kt} keV, band {e1}–{e2} keV, "
              f"β_eq = {beta_eq:.3f}, in-band n_eff(δ=1+β) = {n_eff:.2f}")
        for tag, rot in variant_rotations(anchor.label):
            m, s, _ = per_seed_delta_pf(anchor, mu_centers, seed_curves_at_tau, bend, rot)
            rows[f"{key}_{tag}"] = (m, s)
            print(f"    {VARIANT_NAMES[tag]:<40} ΔPF = {m:+.4f} ± {s:.4f}")
        rows[f"{key}_n_eff"] = (n_eff, 0.0)
        rows[f"{key}_beta_eq"] = (beta_eq, 0.0)
    return rows


def riley_robustness(seed_curves_at_tau, mu_centers):
    """ΔPF(band photon) across band-edge and kT grids — is the C3 number stable?"""
    anchor = next(a for a in j0740.ANCHORS if "Riley" in a.label)
    bend = ExactBending(anchor.compactness)
    r_km = RADIUS_KM[anchor.label]

    print("\n=== Riley robustness: ΔPF(band photon) vs band edges (kT = 0.0842 keV) ===")
    print("    e_min\\e_max " + "".join(f"{e2:>9.2f}" for e2 in E_MAX_SCAN))
    edge_grid = np.empty((len(E_MIN_SCAN), len(E_MAX_SCAN)))
    for i, e1 in enumerate(E_MIN_SCAN):
        cells = []
        for j, e2 in enumerate(E_MAX_SCAN):
            rot = Rotation(J0740_SPIN_HZ, r_km, photon_flux=True,
                           band=BandSpectrum(KT_KEV[RILEY], e1, e2))
            m, _, _ = per_seed_delta_pf(anchor, mu_centers, seed_curves_at_tau, bend, rot)
            edge_grid[i, j] = m
            cells.append(f"{m:+9.4f}")
        print(f"    {e1:>9.2f}  " + "".join(cells))

    print("\n=== Riley robustness: ΔPF(band photon) vs kT (band 0.3–1.5 keV) ===")
    kt_row = np.empty(len(KT_SCAN))
    for i, kt in enumerate(KT_SCAN):
        rot = Rotation(J0740_SPIN_HZ, r_km, photon_flux=True,
                       band=BandSpectrum(kt, *BAND_KEV[RILEY]))
        kt_row[i], _, _ = per_seed_delta_pf(anchor, mu_centers, seed_curves_at_tau,
                                            bend, rot)
        print(f"    kT = {kt:>6.4f} keV   ΔPF = {kt_row[i]:+.4f}")
    return edge_grid, kt_row


def band_shape_routing(pooled_curve, mu_centers):
    """Shape metrics in-band: the routing claim must survive the band treatment."""
    print("\n=== In-band shape routing (pooled curve, band photon) ===")
    print(f"{'anchor':<34}{'ΔPF':>8}{'shapeRMS':>10}{'Δ(A1/A0)':>10}{'Δ(A2/A0)':>10}")
    beaming = beaming_lookup(mu_centers, pooled_curve)
    out = {}
    for anchor in j0740.ANCHORS:
        bend = ExactBending(anchor.compactness)
        rot = Rotation(J0740_SPIN_HZ, RADIUS_KM[anchor.label], photon_flux=True,
                       band=anchor_band(anchor.label))
        iso = multi_spot_flux(anchor.inclination, anchor.compactness, anchor.spots,
                              None, bending=bend, rotation=rot)
        real = multi_spot_flux(anchor.inclination, anchor.compactness, anchor.spots,
                               beaming, bending=bend, rotation=rot)
        rms, _ = waveform_shape_change(iso, real)
        (a1_i, a2_i), (a1_r, a2_r) = harmonics(iso), harmonics(real)
        key = "riley" if "Riley" in anchor.label else "miller"
        out[key] = {"dpf": pulsed_fraction(real) - pulsed_fraction(iso),
                    "shape_rms": rms, "da1": a1_r - a1_i, "da2": a2_r - a2_i}
        print(f"{anchor.label:<34}{out[key]['dpf']:>+8.4f}{rms:>10.4f}"
              f"{out[key]['da1']:>+10.4f}{out[key]['da2']:>+10.4f}")
    return out


def main():
    d = np.load(LIBRARY_PATH)
    mu = d["mu_centers"]
    tau = d["tau_values"]
    ti10 = shape_tau_index(tau)
    seed_curves = d["intensity_by_tau_seed"][ti10]
    pooled = d["intensity_by_tau"][ti10]

    print(f"C3 — band-limited Doppler dilution of ΔPF(τ={SHAPE_TAU:g}), "
          f"{seed_curves.shape[0]} production seeds\n")

    rows = coupling_table(seed_curves, mu)
    edge_grid, kt_row = riley_robustness(seed_curves, mu)
    routing = band_shape_routing(pooled, mu)

    riley_band_m, riley_band_s = rows["riley_band_ph"]
    riley_d4_m, _ = rows["riley_d4"]
    riley_d3_m, _ = rows["riley_d3"]
    print("\n" + "=" * 78)
    print("C3 verdict (Riley, NICER-faithful in-band photon counts):")
    print(f"  ΔPF = {riley_band_m:+.4f} ± {riley_band_s:.4f}   "
          f"[bolometric brackets: δ³ {riley_d3_m:+.4f}, δ⁴ {riley_d4_m:+.4f}; "
          f"0 Hz {rows['riley_0hz'][0]:+.4f}]")
    print(f"  band spread across edge scan: {edge_grid.min():+.4f} … {edge_grid.max():+.4f}; "
          f"kT scan: {kt_row.min():+.4f} … {kt_row.max():+.4f}")
    print(f"  in-band shape routing intact: shape RMS = {routing['riley']['shape_rms']:.3f}, "
          f"Δ(A1/A0) = {routing['riley']['da1']:+.4f}")

    payload = {"tau_shape": np.float64(SHAPE_TAU),
               "spin_hz": np.float64(J0740_SPIN_HZ),
               "riley_kt_kev": np.float64(KT_KEV[RILEY]),
               "miller_kt_kev": np.float64(KT_KEV[MILLER]),
               "riley_band_kev": np.array(BAND_KEV[RILEY]),
               "miller_band_kev": np.array(BAND_KEV[MILLER]),
               "riley_edge_grid": edge_grid,
               "riley_edge_e_min": np.array(E_MIN_SCAN),
               "riley_edge_e_max": np.array(E_MAX_SCAN),
               "riley_kt_scan": np.array(KT_SCAN),
               "riley_kt_row": kt_row}
    for name, (m, s) in rows.items():
        payload[f"j0740_{name}_m"] = np.float64(m)
        payload[f"j0740_{name}_s"] = np.float64(s)
    for key, r in routing.items():
        for f, v in r.items():
            payload[f"j0740_{key}_band_{f}"] = np.float64(v)
    np.savez(RESULTS_PATH, **payload)
    print(f"\n✓ C3 results saved to {RESULTS_PATH}")


if __name__ == "__main__":
    main()
