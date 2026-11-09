"""Revision-round analyses for the APSS submission (referee-driven checks).

Five independent checks, each answering one referee point raised against
``docs/paper/apss_template.tex`` (see the revision plan in that file's
history). Everything reuses the verified pipeline; nothing here touches
``mcrt`` internals.

R1 — toy compactness inversion (the "missing kilometre"):
    Fit the isotropic-emission model to the *beamed* J0740 synthetic profile
    (static, bolometric, exact bending, τ=10 library), letting only the
    compactness u float with the published geometry held fixed. The best-fit
    shift Δu, converted at fixed mass (R = 2GM/uc²), gives the radius bias in
    km that the isotropy approximation would induce in this idealized
    single-parameter fit. A deliberately narrow toy — it bounds the
    order of magnitude, not a posterior shift.

R2 — J0030 at its observed spin (205.53 Hz):
    c2_doppler_coupling already showed ΔPF stays pinned at 0 under rotation
    (saturation is geometric). Here we add the missing number: the spinning
    waveform-shape rms, so the manuscript can quote both NICER sources at
    their observed spins.

R3 — limb-darkening-slope sensitivity (b = 1.72 vs 1.79):
    The validation figure's independent Thomson control quotes b = 1.72; the
    production library's pooled extraction gives b = 1.79. Recompute the
    static J0740 ΔPF with pure linear laws I(μ) ∝ 1 + bμ at both slopes to
    show what the 4% slope discrepancy is worth in ΔPF units.

R4 — exact-bending tiling boundary:
    The phase diagram's analytic boundary Δφ_crit(θ) = f_ecl(θ) is drawn with
    the linear bending map. Recompute f_ecl with the exact Schwarzschild map
    for both regimes so the figure can overlay both and show the criterion's
    placement is bending-robust.

R5 — ΔPF error budget assembly:
    Collect the already-measured systematics (seed scatter, grazing-tail
    perturbation, linear-vs-exact bending, point-spot reduction) into one
    place, so the manuscript can state that ±0.003 is the statistical bar
    and quote the systematic terms alongside it.

Writes ``data/r1_revision_checks.npz`` and prints a summary table.
Run from the repository root:  PYTHONPATH=src:scripts python3 scripts/r1_revision_checks.py
"""

import numpy as np

from mcrt.bending import ExactBending
from mcrt.beaming import beaming_lookup
from mcrt.rotating import Rotation

import j0740_anchor as j0740
import j0030_anchor as j0030
from anchor_lib import (multi_spot_flux, load_library, shape_tau_index,
                        single_spot_eclipsed_fraction, waveform_shape_change)

RESULTS_PATH = "data/r1_revision_checks.npz"

# Published equatorial radii (km) for the fixed-mass u → R conversion and the
# rotation layer; same values as c2_doppler_coupling.RADIUS_KM.
RADIUS_KM = {
    "Riley 2021 (X-PSI ST-U)": 12.39,
    "Miller 2021 (Illinois–Maryland)": 13.713,
    "Riley 2019 (X-PSI ST+PST)": 12.71,
    "Miller 2019 (Illinois–Maryland)": 13.019,
}
J0030_SPIN_HZ = 205.53

# Phase-diagram regimes (must match phase_diagram.py REGIMES).
PD_REGIMES = {
    "j0740": dict(inclination=np.deg2rad(87.6), compactness=0.494),
    "j0030": dict(inclination=np.deg2rad(53.9), compactness=0.312),
}
THETA_AXIS_DEG = np.linspace(0.0, 180.0, 91)


def pulsed_fraction(flux):
    lo, hi = float(flux.min()), float(flux.max())
    return (hi - lo) / (hi + lo) if (hi + lo) > 0 else 0.0


def _profile_rms(a, b):
    """rms difference between two peak-normalized profiles."""
    return float(np.sqrt(np.mean((a / a.max() - b / b.max()) ** 2)))


# ---------------------------------------------------------------------------
# R1 — toy compactness inversion
# ---------------------------------------------------------------------------

def toy_inversion(anchor, beaming, u_lo=0.25, u_hi=0.56, n_coarse=42):
    # u_hi: ExactBending refuses u ≳ 0.57 (grazing deflection ≥ π, multiple
    # imaging); the isotropic fit moves toward LOWER u anyway (flatter model
    # must sharpen), so the cap does not bind.
    """Best-fit isotropic-model compactness to the beamed truth profile.

    Truth: beamed emission at the anchor's published u (exact bending).
    Model: isotropic emission, identical spot geometry, u free.
    Objective: rms between peak-normalized profiles (scale cancels; only the
    pulse shape is fit). Coarse grid then parabolic refinement about the
    minimum — the objective is smooth in u.
    """
    truth = multi_spot_flux(anchor.inclination, anchor.compactness, anchor.spots,
                            beaming, bending=ExactBending(anchor.compactness))

    def objective(u):
        iso = multi_spot_flux(anchor.inclination, u, anchor.spots, None,
                              bending=ExactBending(u))
        return _profile_rms(iso, truth)

    grid = np.linspace(u_lo, u_hi, n_coarse)
    vals = np.array([objective(u) for u in grid])
    k = int(np.argmin(vals))
    # parabolic refinement on the three points around the minimum
    if 0 < k < n_coarse - 1:
        x0, x1, x2 = grid[k - 1:k + 2]
        y0, y1, y2 = vals[k - 1:k + 2]
        denom = (y0 - 2 * y1 + y2)
        u_best = x1 + 0.5 * (x0 - x1) if denom == 0 else \
            x1 - 0.5 * ((x1 - x0) ** 2 * (y1 - y2) - (x1 - x2) ** 2 * (y1 - y0)) / \
            ((x1 - x0) * (y1 - y2) - (x1 - x2) * (y1 - y0))
    else:
        u_best = grid[k]
    rms_best = objective(u_best)
    rms_true_u = objective(anchor.compactness)

    r_true = RADIUS_KM[anchor.label]
    # u = 2GM/Rc² at fixed M  ⇒  R_fit = R_true · (u_true / u_fit)
    r_fit = r_true * anchor.compactness / u_best
    return dict(u_true=anchor.compactness, u_fit=float(u_best),
                du=float(u_best - anchor.compactness),
                r_true=r_true, r_fit=float(r_fit),
                dr_km=float(r_fit - r_true),
                rms_at_fit=rms_best, rms_at_true_u=rms_true_u)


# ---------------------------------------------------------------------------


def main():
    tau_values, mu_centers, intensity_by_tau = load_library()
    ti = shape_tau_index(tau_values)
    beaming = beaming_lookup(mu_centers, intensity_by_tau[ti])
    out = {}

    # ---- R1: toy inversion --------------------------------------------------
    print("=== R1: toy isotropic-fit compactness/radius bias (static, bolometric) ===")
    print(f"{'anchor':<34}{'u_true':>8}{'u_fit':>8}{'Δu':>8}{'R_true':>8}"
          f"{'R_fit':>8}{'ΔR (km)':>9}")
    for anchor in j0740.ANCHORS:
        r = toy_inversion(anchor, beaming)
        key = "riley" if "Riley" in anchor.label else "miller"
        for k, v in r.items():
            out[f"r1_{key}_{k}"] = np.float64(v)
        print(f"{anchor.label:<34}{r['u_true']:>8.3f}{r['u_fit']:>8.3f}"
              f"{r['du']:>+8.3f}{r['r_true']:>8.2f}{r['r_fit']:>8.2f}"
              f"{r['dr_km']:>+9.2f}")
        print(f"    rms(fit u) = {r['rms_at_fit']:.4f}   "
              f"rms(true u) = {r['rms_at_true_u']:.4f}")

    # ---- R2: J0030 spinning shape rms --------------------------------------
    print("\n=== R2: J0030 at 205.53 Hz — spinning shape rms (ΔPF stays 0) ===")
    print(f"{'anchor':<34}{'ΔPF(ν)':>9}{'rms(0)':>9}{'rms(ν)':>9}")
    for anchor in j0030.ANCHORS:
        key = "riley" if "Riley" in anchor.label else "miller"
        rot = Rotation(spin_hz=J0030_SPIN_HZ, radius_km=RADIUS_KM[anchor.label])
        rows = {}
        for tag, r in (("static", None), ("rot", rot)):
            iso = multi_spot_flux(anchor.inclination, anchor.compactness,
                                  anchor.spots, None, rotation=r)
            real = multi_spot_flux(anchor.inclination, anchor.compactness,
                                   anchor.spots, beaming, rotation=r)
            rows[f"dpf_{tag}"] = pulsed_fraction(real) - pulsed_fraction(iso)
            rows[f"rms_{tag}"], _ = waveform_shape_change(iso, real)
        for k, v in rows.items():
            out[f"r2_{key}_{k}"] = np.float64(v)
        print(f"{anchor.label:<34}{rows['dpf_rot']:>+9.4f}"
              f"{rows['rms_static']:>9.4f}{rows['rms_rot']:>9.4f}")

    # ---- R3: slope sensitivity b = 1.72 vs 1.79 ----------------------------
    print("\n=== R3: static J0740 ΔPF under linear laws I ∝ 1 + bμ ===")
    print(f"{'anchor':<34}{'b=1.72':>9}{'b=1.79':>9}{'diff':>8}{'library':>9}")
    for anchor in j0740.ANCHORS:
        key = "riley" if "Riley" in anchor.label else "miller"
        bend = ExactBending(anchor.compactness)
        iso = multi_spot_flux(anchor.inclination, anchor.compactness,
                              anchor.spots, None, bending=bend)
        dpf = {}
        for b in (1.72, 1.79):
            real = multi_spot_flux(anchor.inclination, anchor.compactness,
                                   anchor.spots, lambda mu, b=b: 1.0 + b * np.asarray(mu),
                                   bending=bend)
            dpf[b] = pulsed_fraction(real) - pulsed_fraction(iso)
        lib = multi_spot_flux(anchor.inclination, anchor.compactness,
                              anchor.spots, beaming, bending=bend)
        dpf_lib = pulsed_fraction(lib) - pulsed_fraction(iso)
        out[f"r3_{key}_dpf_b172"] = np.float64(dpf[1.72])
        out[f"r3_{key}_dpf_b179"] = np.float64(dpf[1.79])
        out[f"r3_{key}_dpf_library"] = np.float64(dpf_lib)
        print(f"{anchor.label:<34}{dpf[1.72]:>+9.4f}{dpf[1.79]:>+9.4f}"
              f"{dpf[1.79]-dpf[1.72]:>+8.4f}{dpf_lib:>+9.4f}")

    # ---- R4: exact-bending tiling boundary ---------------------------------
    print("\n=== R4: exact-bending boundary Δφ_crit(θ) = f_ecl(θ) ===")
    theta_axis = np.deg2rad(THETA_AXIS_DEG)
    out["r4_theta_axis_deg"] = THETA_AXIS_DEG
    for key, reg in PD_REGIMES.items():
        bend = ExactBending(reg["compactness"])
        lin = np.array([single_spot_eclipsed_fraction(
            reg["inclination"], t, reg["compactness"]) for t in theta_axis])
        exact = np.array([single_spot_eclipsed_fraction(
            reg["inclination"], t, reg["compactness"], bending=bend)
            for t in theta_axis])
        out[f"r4_{key}_boundary_linear"] = lin
        out[f"r4_{key}_boundary_exact"] = exact
        moved = np.nanmax(np.abs(exact - lin))
        print(f"{key}: max |Δφ_crit(exact) − Δφ_crit(linear)| = {moved:.4f} cycles")

    # ---- R5: error budget assembly -----------------------------------------
    print("\n=== R5: ΔPF error budget (J0740, static headline) ===")
    a3 = np.load("data/a3_seed_errors.npz")
    e2 = np.load("data/e2_tail_sensitivity.npz")
    b2 = np.load("data/b2_exact_bending.npz")
    fc = np.load("data/finite_cap.npz")
    tail = float(np.max([np.max(np.abs(np.atleast_1d(e2[k]).astype(float)))
                         for k in e2.files if k.endswith("tail_systematic")]))
    budget = {
        "seed_riley": float(a3["j0740_riley_dpf_std"][shape_tau_index(a3["tau_values"])]),
        "seed_miller": float(a3["j0740_miller_dpf_std"][shape_tau_index(a3["tau_values"])]),
        "tail_max": tail,
        "bending_linear_bias_riley": float(np.atleast_1d(b2["j0740_riley_shift"])[0]),
        "bending_linear_bias_miller": float(np.atleast_1d(b2["j0740_miller_shift"])[0]),
        "pointspot_riley": float(fc["j0740_riley_point"]) - float(fc["j0740_riley_converged"]),
        "pointspot_miller": float(fc["j0740_miller_point"]) - float(fc["j0740_miller_converged"]),
    }
    for k, v in budget.items():
        out[f"r5_{k}"] = np.float64(v)
        print(f"  {k:<22}{v:+.4f}")

    # ---- R6: Vinciguerra+2024 J0030 reanalysis modes ------------------------
    # Medians from Vinciguerra et al. 2024 (ApJ 961, 62) online Figure set 6
    # (the hot-region parameters are not printed in the paper's tables).
    # u = 2*1.47662*(M/Msun)/R_km; i = arccos(cos i); the secondary phase is
    # defined "plus half a cycle" (their Table 1), so
    # dphi = (phi_s + 0.5) - phi_p.
    # Both modes are PDT (protruding double temperature): a small hot cap
    # embedded off-center in a larger cool cap. That does NOT reduce cleanly to
    # a single-temperature point spot, so this check is indicative only --- it
    # asks which side of the tiling boundary the *geometry* falls on, using
    # the superseding (hot) cap as the point-spot stand-in.
    print("\n=== R6: Vinciguerra+2024 J0030 modes vs the tiling criterion ===")
    v_modes = {
        "pdtu_joint": dict(  # highest evidence (log E = -42677)
            label="PDT-U (NICER+XMM)", m=1.70, r_km=14.438, cos_i=0.26,
            theta_p=1.31, theta_s=1.48, phi_p=0.4258, phi_s=0.4659),
        "stpdt_joint": dict(  # second, also "well reproduces" the data
            label="ST+PDT (NICER+XMM)", m=1.40, r_km=11.71, cos_i=0.137,
            theta_p=1.76, theta_s=1.56, phi_p=0.4380, phi_s=0.4818),
    }
    print(f"{'mode':<22}{'u':>7}{'i(deg)':>8}{'Dphi':>7}"
          f"{'Dphi_crit':>11}{'floor':>8}{'dPF':>9}{'rms':>8}")
    for key, m in v_modes.items():
        u = 2.0 * 1.47662 * m["m"] / m["r_km"]
        incl = float(np.arccos(m["cos_i"]))
        dphi = (m["phi_s"] + 0.5) - m["phi_p"]
        dphi = min(dphi % 1.0, 1.0 - (dphi % 1.0))  # fold to [0, 0.5]
        spots = (j0740.Spot(colatitude=m["theta_p"], azimuth=0.0, weight=1.0),
                 j0740.Spot(colatitude=m["theta_s"], azimuth=dphi, weight=1.0))
        bend = ExactBending(u)
        iso = multi_spot_flux(incl, u, spots, None, bending=bend)
        real = multi_spot_flux(incl, u, spots, beaming, bending=bend)
        floor = float(iso.min() / iso.max())
        dpf = pulsed_fraction(real) - pulsed_fraction(iso)
        rms, _ = waveform_shape_change(iso, real)
        # critical separation at the mean colatitude of the two spots
        theta_mean = 0.5 * (m["theta_p"] + m["theta_s"])
        dphi_crit = single_spot_eclipsed_fraction(incl, theta_mean, u,
                                                  bending=bend)
        for name, val in (("u", u), ("incl", incl), ("dphi", dphi),
                          ("dphi_crit", dphi_crit), ("floor", floor),
                          ("dpf", dpf), ("rms", rms)):
            out[f"r6_{key}_{name}"] = np.float64(val)
        print(f"{m['label']:<22}{u:>7.3f}{np.rad2deg(incl):>8.1f}{dphi:>7.3f}"
              f"{dphi_crit:>11.3f}{floor:>8.3f}{dpf:>+9.4f}{rms:>8.4f}")
    print("  (floor > 0 => tiles the rotation => PF-sensitive regime)")

    np.savez(RESULTS_PATH, **out)
    print(f"\nwrote {RESULTS_PATH}")


if __name__ == "__main__":
    main()
