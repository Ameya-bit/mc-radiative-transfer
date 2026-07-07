"""Track C4 — caveat re-audit against the *diluted* ΔPF (v0.9.10 close-out).

C2/C3 shrank the J0740 headline from +0.137/+0.195 to +0.02–0.06
(convention/band-dependent). Two geometry-only effects were dismissed as
negligible *relative to 0.14*; C4 re-audits them against the diluted number, so
the paper's caveat sentences carry numbers sized to the claim they guard.

1. **Photon light-travel-time delay** (Bogdanov eq. 18–19). The observed-phase
   warp Δφ = ν·Δt(α) is beaming-independent and *exactly* PF-invariant for a
   single spot (a pure phase reparametrization preserves max and min). For a
   **multi-spot sum** the warp does not commute with the summation, so its
   cancellation in ΔPF is only second-order. C4 measures the residual directly:
   per-seed ΔPF with and without the warp applied to each spot's profile before
   the roll-and-add (warp-then-roll is exact — the delay depends only on the
   spot's own emission angle, a function of φ − φ_spot), at 346.5 Hz under both
   the δ⁴ bolometric and the band-photon (C3 headline) treatments. The shift is
   evaluated **per seed, paired**, so seed noise cancels in the difference.
   J0030 is skipped: its PF is pinned at 1 (dark gap) under any phase warp.

2. **Rotational oblateness** (AlGendy & Morsink 2014). Not modelled exactly (no
   oblate ray-tracing — the OS-level extension); instead the *leading point-spot
   effects* are bounded by direct perturbation. AM14 shape:
   R(θ) = R_eq·(1 + o₂ cos²θ), o₂ = Ω̄²(−0.788 + 1.030·x), with x = GM/(R_eq c²)
   = u/2 and Ω̄² = Ω²R_eq³/(GM). At J0740, Ω̄² ≈ 0.033 → o₂ ≈ −1.7% pole
   flattening — but the anchors' spots are near-equatorial (cos²θ ≤ 0.10), so the
   spot-local radius shifts by ≲ 0.2% and the surface normal tilts poleward by
   |η| = |o₂ sin 2θ| ≲ 0.5°. C4 applies both perturbations (per-spot colatitude
   tilt θ → θ + η; worst-spot common rescale of u and the rotation radius) and
   reports the paired per-seed ΔPF shift. The unmodelled remainder (oblate
   bending/Jacobian differences) is an OS-level correction of the same order —
   quoted as the residual caveat.

Verdict rule per effect: |shift| ≤ seed σ ⇒ the caveat sentence stands, with this
measured bound attached; larger ⇒ escalate (model it properly).

Deterministic, seconds to run. Run from the repository root:
    PYTHONPATH=src:scripts python3 scripts/c4_caveat_audit.py
"""

import numpy as np

from mcrt import (
    ExactBending,
    Rotation,
    beaming_lookup,
    compute_profile,
    pulsed_fraction,
    travel_time_delay,
)
from anchor_lib import N_PHASE, SHAPE_TAU, Spot, multi_spot_flux, shape_tau_index
from c2_doppler_coupling import J0740_SPIN_HZ, RADIUS_KM
from c3_band_doppler import anchor_band
import j0740_anchor as j0740

LIBRARY_PATH = "data/beaming_library.npz"
RESULTS_PATH = "data/c4_caveat_audit.npz"

# AM14 o₂ coefficients (their eq. 20 fit): o₂ = Ω̄²(A + B·x).
AM14_A, AM14_B = -0.788, 1.030


def delayed_multi_spot_flux(anchor, beaming, bending, rotation, radius_km, spin_hz):
    """`multi_spot_flux` + the eq. 18–19 phase warp applied per spot before the roll."""
    total = np.zeros(N_PHASE)
    uniform = np.arange(N_PHASE) / N_PHASE
    for spot in anchor.spots:
        prof = compute_profile(anchor.inclination, spot.colatitude, anchor.compactness,
                               beaming=beaming, n_phase=N_PHASE, bending=bending,
                               rotation=rotation)
        delay = travel_time_delay(prof.cos_alpha, anchor.compactness, radius_km)
        obs = (prof.phase / (2.0 * np.pi) + spin_hz * delay) % 1.0
        order = np.argsort(obs)
        xs, ys = obs[order], prof.flux[order]
        warped = np.interp(uniform, np.concatenate([xs - 1.0, xs, xs + 1.0]),
                           np.concatenate([ys, ys, ys]))
        shift = int(round(spot.azimuth * N_PHASE)) % N_PHASE
        total += spot.weight * np.roll(warped, shift)
    return total


def per_seed_dpf(anchor, mu_centers, seed_curves, flux_fn) -> np.ndarray:
    """Per-seed ΔPF array under an arbitrary flux builder (paired-shift statistics)."""
    pf_iso = pulsed_fraction(flux_fn(anchor, None))
    dpf = np.empty(seed_curves.shape[0])
    for si in range(seed_curves.shape[0]):
        real = flux_fn(anchor, beaming_lookup(mu_centers, seed_curves[si]))
        dpf[si] = pulsed_fraction(real) - pf_iso
    return dpf


def delay_audit(seed_curves, mu_centers):
    """Measured multi-spot time-delay shift of ΔPF, per treatment, paired per seed."""
    print("=== C4.1 light-travel-time delay: ΔPF(warp) − ΔPF(no warp), paired per seed ===")
    print(f"{'anchor':<34}{'treatment':>12}{'ΔPF(no warp)':>15}{'shift':>10}{'|shift|/σ':>11}")
    rows = {}
    for anchor in j0740.ANCHORS:
        bend = ExactBending(anchor.compactness)
        r_km = RADIUS_KM[anchor.label]
        key = "riley" if "Riley" in anchor.label else "miller"
        for tag, rot in (("d4", Rotation(J0740_SPIN_HZ, r_km)),
                         ("band_ph", Rotation(J0740_SPIN_HZ, r_km, photon_flux=True,
                                              band=anchor_band(anchor.label)))):
            def plain(a, b):
                return multi_spot_flux(a.inclination, a.compactness, a.spots, b,
                                       bending=bend, rotation=rot)

            def warped(a, b):
                return delayed_multi_spot_flux(a, b, bend, rot, r_km, J0740_SPIN_HZ)

            base = per_seed_dpf(anchor, mu_centers, seed_curves, plain)
            with_delay = per_seed_dpf(anchor, mu_centers, seed_curves, warped)
            paired = with_delay - base
            sigma = base.std(ddof=1)
            ratio = abs(paired.mean()) / sigma if sigma > 0 else np.inf
            rows[f"{key}_delay_{tag}"] = (base.mean(), paired.mean(), sigma, ratio)
            print(f"{anchor.label:<34}{tag:>12}"
                  + f"{base.mean():+.4f}".rjust(15)
                  + f"{paired.mean():+.5f}".rjust(10)
                  + f"{ratio:.2f}σ".rjust(11))
    return rows


def oblate_perturbed(anchor, radius_km):
    """(perturbed anchor, perturbed radius, diagnostics) under AM14 first order.

    Per-spot poleward normal tilt η = o₂ sin 2θ (o₂ < 0 tilts toward the nearer
    pole on both hemispheres); worst-spot |ΔR/R| = |o₂| max(cos²θ) applied as a
    common rescale of u (∝ 1/R) and the rotation radius. First-order bound, not an
    oblate ray-trace.
    """
    x = anchor.compactness / 2.0
    omega = 2.0 * np.pi * J0740_SPIN_HZ
    # Ω̄² = Ω²R³/GM with GM = u·R·c²/2 (c in km/s, R in km) ⇒ Ω̄² = 2Ω²R²/(u·c²).
    ombar2 = 2.0 * omega**2 * radius_km**2 / (anchor.compactness * 299792.458**2)
    o2 = ombar2 * (AM14_A + AM14_B * x)
    spots = [Spot(colatitude=s.colatitude + o2 * np.sin(2.0 * s.colatitude),
                  azimuth=s.azimuth, weight=s.weight) for s in anchor.spots]
    dr_over_r = o2 * max(np.cos(s.colatitude) ** 2 for s in anchor.spots)  # signed, worst
    u_pert = anchor.compactness / (1.0 + dr_over_r)
    r_pert = radius_km * (1.0 + dr_over_r)
    pert = anchor._replace(spots=spots, compactness=u_pert)
    tilts_deg = [np.rad2deg(o2 * np.sin(2.0 * s.colatitude)) for s in anchor.spots]
    return pert, r_pert, o2, ombar2, tilts_deg, dr_over_r


def oblateness_audit(seed_curves, mu_centers):
    """Paired per-seed ΔPF shift under the AM14 first-order spot perturbations (δ⁴)."""
    print("\n=== C4.2 oblateness (AM14 first order): perturbed − spherical, paired per seed ===")
    print(f"{'anchor':<34}{'o₂':>8}{'tilts (°)':>16}{'ΔR/R':>9}{'shift':>10}{'|shift|/σ':>11}")
    rows = {}
    for anchor in j0740.ANCHORS:
        r_km = RADIUS_KM[anchor.label]
        key = "riley" if "Riley" in anchor.label else "miller"
        pert, r_pert, o2, ombar2, tilts_deg, dr = oblate_perturbed(anchor, r_km)
        bend0, bend1 = ExactBending(anchor.compactness), ExactBending(pert.compactness)
        rot0, rot1 = Rotation(J0740_SPIN_HZ, r_km), Rotation(J0740_SPIN_HZ, r_pert)

        def spherical(a, b):
            return multi_spot_flux(a.inclination, a.compactness, a.spots, b,
                                   bending=bend0, rotation=rot0)

        def oblate(_, b):
            return multi_spot_flux(pert.inclination, pert.compactness, pert.spots, b,
                                   bending=bend1, rotation=rot1)

        base = per_seed_dpf(anchor, mu_centers, seed_curves, spherical)
        shifted = per_seed_dpf(anchor, mu_centers, seed_curves, oblate)
        paired = shifted - base
        sigma = base.std(ddof=1)
        ratio = abs(paired.mean()) / sigma if sigma > 0 else np.inf
        rows[f"{key}_oblate"] = (paired.mean(), sigma, ratio, o2, ombar2, dr)
        tilt_str = "/".join(f"{t:+.2f}" for t in tilts_deg)
        print(f"{anchor.label:<34}{o2:>8.4f}{tilt_str:>16}{dr:>9.5f}"
              + f"{paired.mean():+.5f}".rjust(10) + f"{ratio:.2f}σ".rjust(11))
    return rows


def main():
    d = np.load(LIBRARY_PATH)
    mu = d["mu_centers"]
    ti10 = shape_tau_index(d["tau_values"])
    seed_curves = d["intensity_by_tau_seed"][ti10]

    print(f"C4 — caveat re-audit vs the diluted ΔPF(τ={SHAPE_TAU:g}), "
          f"{seed_curves.shape[0]} production seeds, exact bending, 346.53 Hz\n")

    delay_rows = delay_audit(seed_curves, mu)
    oblate_rows = oblateness_audit(seed_curves, mu)

    print("\n" + "=" * 78)
    print("C4 verdict (all shifts paired per seed; σ = seed error of ΔPF):")
    worst_delay = max(abs(v[1]) for v in delay_rows.values())
    worst_delay_ratio = max(v[3] for v in delay_rows.values())
    worst_oblate = max(abs(v[0]) for v in oblate_rows.values())
    worst_oblate_ratio = max(v[2] for v in oblate_rows.values())
    print(f"  time delay: worst |shift| = {worst_delay:.5f} ({worst_delay_ratio:.2f}σ) — "
          f"{'caveat STANDS' if worst_delay_ratio <= 1.0 else 'ESCALATE'}")
    print(f"  oblateness (first order): worst |shift| = {worst_oblate:.5f} "
          f"({worst_oblate_ratio:.2f}σ) — "
          f"{'caveat STANDS' if worst_oblate_ratio <= 1.0 else 'ESCALATE'}"
          "\n  (residual: oblate bending/Jacobian differences are OS-level, same order)")

    payload = {"tau_shape": np.float64(SHAPE_TAU), "spin_hz": np.float64(J0740_SPIN_HZ)}
    for name, (base_m, shift, sigma, ratio) in delay_rows.items():
        payload[f"j0740_{name}_base"] = np.float64(base_m)
        payload[f"j0740_{name}_shift"] = np.float64(shift)
        payload[f"j0740_{name}_sigma"] = np.float64(sigma)
        payload[f"j0740_{name}_ratio"] = np.float64(ratio)
    for name, (shift, sigma, ratio, o2, ombar2, dr) in oblate_rows.items():
        payload[f"j0740_{name}_shift"] = np.float64(shift)
        payload[f"j0740_{name}_sigma"] = np.float64(sigma)
        payload[f"j0740_{name}_ratio"] = np.float64(ratio)
        payload[f"j0740_{name}_o2"] = np.float64(o2)
        payload[f"j0740_{name}_ombar2"] = np.float64(ombar2)
        payload[f"j0740_{name}_dr_over_r"] = np.float64(dr)
    np.savez(RESULTS_PATH, **payload)
    print(f"\n✓ C4 results saved to {RESULTS_PATH}")


if __name__ == "__main__":
    main()
