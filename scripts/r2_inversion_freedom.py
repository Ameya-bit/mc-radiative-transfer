"""How much of the toy inversion's Δu survives extra model freedom?

Revision-round-2 check behind the rewritten Section 3.1. The R1 toy
(``r1_revision_checks.toy_inversion``) fits the beamed J0740 profile with the
isotropic-emission model letting only the compactness ``u`` float, and finds
Δu ≈ −0.13. Two objections follow immediately, and this script answers both
with the same verified pipeline:

R2a — *is Δu just an artifact of freezing everything else?*
    Re-run the fit with, in turn, the observer inclination, a common
    spot-colatitude offset, an additive unpulsed background, and the relative
    weight of the two spots released — individually and jointly. Δu moves by
    at most ~0.02, so the shift is a property of the beaming mismatch, not of
    the single-parameter restriction.

R2b — *does the fit ever actually become good?*
    Track the residual at each level of freedom. It floors near 0.035 rms:
    the isotropic model is structurally unable to reproduce a limb-darkened
    waveform, and the compactness it prefers (u ≈ 0.30–0.36, i.e. R ≈ 17–20 km
    at fixed mass) lies outside every viable dense-matter equation of state
    and outside the published analyses' prior support. That is the basis for
    the manuscript's "detectable rather than quietly biasing" framing, and for
    dropping the earlier 4–6 km radius-bias claim.

R2c — sign cross-check against Zhao et al. (2025, ApJ 982, 112).
    They report that a model *more* peaked than the truth biases compactness
    *down* (δu/u ≃ −δh/h). Run both mis-specification directions through this
    pipeline. It is monotonic in the opposite sense, because the point-spot
    reduction carries no cap-area parameter (their q, s) for the fit to move —
    the reconciliation stated in the Discussion.

Writes ``data/r2_inversion_freedom.npz`` and prints a summary table.
Run from the repository root:
    PYTHONPATH=src:scripts python3 scripts/r2_inversion_freedom.py
"""

import numpy as np

from mcrt.bending import ExactBending
from mcrt.beaming import beaming_lookup

import j0740_anchor as j0740
from anchor_lib import multi_spot_flux, load_library, shape_tau_index, Spot

RESULTS_PATH = "data/r2_inversion_freedom.npz"

# Same fixed-mass u → R conversion as r1_revision_checks.
RADIUS_KM = {
    "Riley 2021 (X-PSI ST-U)": 12.39,
    "Miller 2021 (Illinois–Maryland)": 13.713,
}

# ExactBending refuses u ≳ 0.57 (grazing deflection ≥ π); the isotropic fit
# moves toward *lower* u anyway, so the upper cap never binds in R2a/R2b.
U_LO, U_HI = 0.26, 0.56


def profile_rms(a, b):
    """rms difference between two peak-normalized profiles."""
    return float(np.sqrt(np.mean((a / a.max() - b / b.max()) ** 2)))


def minimize_grid(objective, lo, hi, rounds=6, n=11, sweeps=8, seed=None):
    """Successive-refinement minimisation on a box (no scipy in this repo).

    One and two dimensions are cheap enough for a full product grid, refined
    ``rounds`` times about the incumbent (halving the span each round, so
    ``rounds=6`` resolves u to ~1e-4). Above two dimensions a product grid is
    prohibitive, so we run cyclic coordinate descent instead, seeded from
    ``seed`` (the lower-dimensional optimum) rather than the box centre — the
    objective has a long flat valley in u and a centre start does not reliably
    find it. The caller checks that the joint residual is no worse than the
    nested subsets; if it is, the descent has stalled and the result is not
    reported.
    """
    lo, hi = np.array(lo, float), np.array(hi, float)
    if len(lo) <= 2:
        best_x, best_f = None, np.inf
        for _ in range(rounds):
            axes = [np.linspace(a, b, n) for a, b in zip(lo, hi)]
            grid = np.array(np.meshgrid(*axes, indexing="ij")).reshape(len(lo), -1).T
            for point in grid:
                value = objective(point)
                if value < best_f:
                    best_f, best_x = value, point.copy()
            span = (hi - lo) / (n - 1)
            lo, hi = best_x - span, best_x + span
        return best_x, best_f

    best_x = 0.5 * (lo + hi) if seed is None else np.clip(np.array(seed, float), lo, hi)
    best_f = objective(best_x)
    span = 0.5 * (hi - lo)
    for _ in range(sweeps):
        improved = False
        for k in range(len(lo)):
            for value in np.linspace(best_x[k] - span[k], best_x[k] + span[k], n):
                trial = best_x.copy()
                trial[k] = np.clip(value, lo[k], hi[k])
                f = objective(trial)
                if f < best_f - 1e-12:
                    best_f, best_x, improved = f, trial, True
        span = span / 2.0 if improved else span / 4.0
    return best_x, best_f


def isotropic_model(anchor, u, incl=None, dtheta=0.0, weight_ratio=None,
                    background=0.0):
    """Isotropic-emission profile at compactness ``u`` with optional freedom.

    ``dtheta`` shifts both spot colatitudes by a common offset, ``weight_ratio``
    re-weights the second spot relative to the first, and ``background`` adds a
    phase-independent floor as a fraction of the model's own peak. Every knob
    defaults to the published value, so calling with ``u`` alone reproduces the
    R1 toy bit-for-bit.
    """
    spots = anchor.spots
    if dtheta:
        spots = tuple(s._replace(
            colatitude=float(np.clip(s.colatitude + dtheta, 1e-3, np.pi - 1e-3)))
            for s in spots)
    if weight_ratio is not None:
        first = spots[0]
        spots = (first, spots[1]._replace(weight=first.weight * weight_ratio))
    flux = multi_spot_flux(anchor.inclination if incl is None else float(incl),
                           float(u), spots, None, bending=ExactBending(float(u)))
    return flux + background * flux.max() if background else flux


def main():
    tau_values, mu_centers, intensity_by_tau = load_library()
    beaming = beaming_lookup(mu_centers, intensity_by_tau[shape_tau_index(tau_values)])
    out = {}

    print("=== R2a/R2b: does Δu survive extra freedom, and does the fit ever fit? ===")
    for anchor in j0740.ANCHORS:
        key = "riley" if "Riley" in anchor.label else "miller"
        u0, i0 = anchor.compactness, anchor.inclination
        r_true = RADIUS_KM[anchor.label]
        truth = multi_spot_flux(i0, u0, anchor.spots, beaming,
                                bending=ExactBending(u0))

        baseline = profile_rms(isotropic_model(anchor, u0), truth)
        out[f"r2_{key}_rms_no_freedom"] = np.float64(baseline)

        # (label, parameter box, unpacker → kwargs for isotropic_model)
        cases = (
            ("u only", [U_LO], [U_HI],
             lambda p: dict(u=p[0])),
            ("u + inclination", [U_LO, np.deg2rad(45)], [U_HI, np.pi / 2],
             lambda p: dict(u=p[0], incl=np.clip(p[1], 0.3, np.pi / 2))),
            ("u + colatitude", [U_LO, -0.4], [U_HI, 0.4],
             lambda p: dict(u=p[0], dtheta=float(np.clip(p[1], -0.6, 0.6)))),
            ("u + background", [U_LO, 0.0], [U_HI, 3.0],
             lambda p: dict(u=p[0], background=max(float(p[1]), 0.0))),
            ("u + weight ratio", [U_LO, 0.05], [U_HI, 4.0],
             lambda p: dict(u=p[0], weight_ratio=max(float(p[1]), 0.02))),
            ("u + all four", [U_LO, np.deg2rad(45), -0.4, 0.0, 0.05],
             [U_HI, np.pi / 2, 0.4, 3.0, 4.0],
             lambda p: dict(u=p[0], incl=np.clip(p[1], 0.3, np.pi / 2),
                            dtheta=float(np.clip(p[2], -0.6, 0.6)),
                            background=max(float(p[3]), 0.0),
                            weight_ratio=max(float(p[4]), 0.02))),
        )

        print(f"\n{anchor.label}   u_true={u0:.3f}  R_true={r_true:.2f} km")
        print(f"    isotropic model at the true u, no freedom: rms = {baseline:.4f}")
        print(f"    {'freedom':<20}{'u_fit':>8}{'Δu':>9}{'R_fit':>9}{'ΔR (km)':>10}{'rms':>9}")
        u_only, best_subset_rms = None, np.inf
        for name, lo, hi, unpack in cases:
            def objective(p, unpack=unpack):
                kw = unpack(p)
                kw["u"] = float(np.clip(kw["u"], U_LO, U_HI))
                return profile_rms(isotropic_model(anchor, **kw), truth)

            # Seed the joint descent at the u-only optimum with the remaining
            # axes neutral (published inclination, no offset/background,
            # published weight ratio), not at the box centre.
            seed = None
            if len(lo) > 2 and u_only is not None:
                seed = [u_only, i0, 0.0, 0.0, anchor.spots[1].weight / anchor.spots[0].weight]
            point, rms = minimize_grid(objective, lo, hi, seed=seed)
            if len(lo) > 2 and rms > best_subset_rms + 1e-6:
                raise RuntimeError(
                    f"joint descent stalled for {anchor.label}: rms {rms:.4f} is worse "
                    f"than the best nested subset ({best_subset_rms:.4f})")
            best_subset_rms = min(best_subset_rms, rms)
            u_fit = float(np.clip(unpack(point)["u"], U_LO, U_HI))
            if u_only is None:
                u_only = u_fit
            r_fit = r_true * u0 / u_fit
            tag = name.replace(" ", "_").replace("+", "plus")
            out[f"r2_{key}_{tag}_u"] = np.float64(u_fit)
            out[f"r2_{key}_{tag}_rms"] = np.float64(rms)
            out[f"r2_{key}_{tag}_dr_km"] = np.float64(r_fit - r_true)
            print(f"    {name:<20}{u_fit:>8.3f}{u_fit - u0:>+9.3f}"
                  f"{r_fit:>9.2f}{r_fit - r_true:>+10.2f}{rms:>9.4f}")

    # ---- R2c: sign cross-check against Zhao et al. (2025) -------------------
    print("\n=== R2c: sign of the compactness shift for both mis-specifications ===")
    print("    Zhao et al. (2025): model MORE peaked than truth ⇒ u biased DOWN.")
    strong = lambda mu: 1.0 + 4.0 * np.asarray(mu)   # noqa: E731 (b=4, well above the slab's 1.79)
    directions = (
        ("beamed truth  -> isotropic model", beaming, None),
        ("isotropic truth -> beamed model", None, beaming),
        ("isotropic truth -> strong b=4 model", None, strong),
        ("beamed truth  -> strong b=4 model", beaming, strong),
    )
    print(f"    {'direction':<38}{'u_fit':>8}{'Δu':>9}{'rms':>9}   sense")
    for anchor in j0740.ANCHORS:
        key = "riley" if "Riley" in anchor.label else "miller"
        u0 = anchor.compactness
        print(f"  {anchor.label}  (u_true = {u0:.3f})")
        for name, truth_beaming, model_beaming in directions:
            truth = multi_spot_flux(anchor.inclination, u0, anchor.spots,
                                    truth_beaming, bending=ExactBending(u0))
            point, rms = minimize_grid(
                lambda p: profile_rms(
                    multi_spot_flux(anchor.inclination,
                                    float(np.clip(p[0], U_LO, U_HI)), anchor.spots,
                                    model_beaming,
                                    bending=ExactBending(float(np.clip(p[0], U_LO, U_HI)))),
                    truth),
                [U_LO], [U_HI])
            u_fit = float(np.clip(point[0], U_LO, U_HI))
            tag = name.split("->")[1].strip().replace(" ", "_").replace("=", "")
            out[f"r2c_{key}_{tag}_u"] = np.float64(u_fit)
            sense = "u DOWN (R up)" if u_fit < u0 else "u UP (R down)"
            print(f"    {name:<38}{u_fit:>8.3f}{u_fit - u0:>+9.3f}{rms:>9.4f}   {sense}")

    np.savez(RESULTS_PATH, **out)
    print(f"\nwrote {RESULTS_PATH}")


if __name__ == "__main__":
    main()
