"""Validate the Doppler + aberration layer against NICER code-comparison SD1c (v0.9.10).

`code_comparison.py` reproduced test **SD1a** (1 Hz) — the pure light-bending limit,
where the rotational velocity is negligible. This script closes the loop on the new
rotational layer (`mcrt.rotating`) by reproducing test **SD1c**, which is *the same
geometry as SD1a but spinning at 200 Hz* (Bogdanov et al. 2019, ApJL 887 L26,
Table 1):

  - point spot (angular radius 0.01 rad), isotropic Planck emission, kT = 0.35 keV;
  - i = θ_s = 90° (equatorial spot, equatorial observer → maximal Doppler);
  - M = 1.4 M_sun, R = 12 km ⇒ u = 2GM/Rc² ≈ 0.3445;
  - ν = 200 Hz ⇒ β = 2πνR/(c√(1−u)) ≈ 0.062.

Because SD1a and SD1c share geometry, **SD1c − SD1a isolates the pure Doppler +
aberration effect** — and the reference pair lets us check that our S+D layer
reproduces it, not just that it is internally consistent.

The observable is the monochromatic photon flux at 1 keV (photons/cm²/s/keV). For a
blackbody with isotropic beaming the Bogdanov eq. (20) flux, written for the shape
(all φ-independent constants dropped), collapses cleanly: the δ³ boost cancels the
E'³ of the Planck head, leaving

    f(φ) ∝ cos α'(φ) · D(φ) / (exp[E'(φ)/kT] − 1),   E'(φ) = 1 keV / (δ(φ) √(1−u)),

with cos α' = δ cos α the aberrated projection (eq. 13) and D = d cosα/d cosψ the
lensing factor. All of the rotation enters through δ(φ): the aberration in cos α'
and the *spectral* shift in E'(φ) (steep on the Wien tail, since 1 keV ≫ kT), which
is what tilts the 200 Hz pulse relative to the symmetric 1 Hz one. This uses the
*same* `spot_speed`/`cos_xi`/`doppler_factor` primitives as the production flux.

One extra ingredient is needed *here* that the production layer deliberately omits:
the **photon light-travel-time delay** (Bogdanov eq. 18–19). The reference SD1c
waveform includes it, and at 200 Hz it warps the phase by up to ~0.008 cycles —
enough to leave several-percent residuals at the steep eclipse edges if ignored. It
is a *geometry-only* term (independent of the beaming law), so it stays out of the
production flux; the quadrature lives in `mcrt.rotating.travel_time_delay` as a
driver-level utility (used here and by the C4 caveat audit, which bounds its
second-order multi-spot effect on ΔPF). With it applied, our only remaining approximation vs. the IM
code is the Beloborodov (2002) linear bending map (≈1% at u ≈ 0.34), so — as for
SD1a — agreement is expected at the ~1% level, worst at the grazing eclipse edge.

Reference data: ``data/l26_reference/SD1c_test_IM.txt`` (and SD1a for contrast) —
two columns: rotational phase (cycles) and photon flux (photons/cm²/s/keV).

Run from the repository root:  python3 scripts/c1_doppler_validate.py
"""

import numpy as np
import matplotlib.pyplot as plt

from mcrt import bend, cos_psi, cos_xi, doppler_factor, spot_speed, travel_time_delay

# Shared SD1 parameters (Bogdanov et al. 2019, Table 1 + §4.2 preamble).
M_SUN_KM = 1.47662                 # GM_sun/c² in km
NS_MASS_MSUN = 1.4
NS_RADIUS_KM = 12.0
INCLINATION = np.deg2rad(90.0)     # colatitude of observer ζ
COLATITUDE = np.deg2rad(90.0)      # colatitude of spot center θ_c
COMPACTNESS = 2.0 * M_SUN_KM * NS_MASS_MSUN / NS_RADIUS_KM   # u ≈ 0.3445
KT_KEV = 0.35                      # comoving blackbody temperature
E_OBS_KEV = 1.0                    # observed energy at which the flux is tabulated

SPIN_SD1A_HZ = 1.0                 # SD1a: the near-static reference
SPIN_SD1C_HZ = 200.0               # SD1c: same geometry, rotating

REF_SD1A = "data/l26_reference/SD1a_test_IM.txt"
REF_SD1C = "data/l26_reference/SD1c_test_IM.txt"
FIGURE_PATH = "data/pulse_profile_doppler_sd1c.png"


def _emission_flux(phi: np.ndarray, spin_hz: float) -> np.ndarray:
    """Collapsed Bogdanov eq. (20) monochromatic flux at *emission* phases φ (rad).

    Blackbody + isotropic spot; zero where the spot is eclipsed (cos α < 0). Does
    NOT include the light-travel-time phase warp — that is applied downstream in
    :func:`monochromatic_flux`.
    """
    cos_psi_vals = cos_psi(phi, INCLINATION, COLATITUDE)
    cos_alpha = bend(cos_psi_vals, COMPACTNESS)            # linear Beloborodov map
    lensing = 1.0 - COMPACTNESS                            # D = d cosα/d cosψ (linear)
    visible = cos_alpha >= 0.0

    beta = spot_speed(spin_hz, NS_RADIUS_KM, COLATITUDE, COMPACTNESS)
    delta = doppler_factor(beta, cos_xi(phi, cos_alpha, cos_psi_vals, INCLINATION))

    cos_alpha_comoving = delta * cos_alpha                 # aberrated projection (eq. 13)
    e_prime = E_OBS_KEV / (delta * np.sqrt(1.0 - COMPACTNESS))   # comoving energy (eq. 15)
    planck_tail = 1.0 / np.expm1(e_prime / KT_KEV)         # 1/(exp(E'/kT) − 1)

    flux = cos_alpha_comoving * lensing * planck_tail
    return np.where(visible, flux, 0.0)


def monochromatic_flux(phase_cycles: np.ndarray, spin_hz: float,
                       with_time_delay: bool = True) -> np.ndarray:
    """Observed 1 keV photon-flux shape at the reference phases (cycles).

    Builds the emission-phase flux on a dense grid, applies the light-travel-time
    phase lag ``Δφ = ν Δt(α)`` [cycles] (eq. 19; toggle with ``with_time_delay``),
    and resamples periodically onto ``phase_cycles``. Set ``with_time_delay=False``
    to see the residual the omitted delay would leave — the production layer's regime.
    """
    n = 4096
    phi = np.linspace(0.0, 2.0 * np.pi, n, endpoint=False)
    f_emit = _emission_flux(phi, spin_hz)

    obs_cycles = phi / (2.0 * np.pi)
    if with_time_delay:
        cos_alpha = bend(cos_psi(phi, INCLINATION, COLATITUDE), COMPACTNESS)
        obs_cycles = obs_cycles + spin_hz * travel_time_delay(
            cos_alpha, COMPACTNESS, NS_RADIUS_KM)                          # Δφ[cyc] = ν Δt

    # Resample onto the requested phases, periodic in one cycle.
    order = np.argsort(obs_cycles % 1.0)
    xs = (obs_cycles % 1.0)[order]
    ys = f_emit[order]
    xs_ext = np.concatenate([xs - 1.0, xs, xs + 1.0])
    ys_ext = np.concatenate([ys, ys, ys])
    return np.interp(phase_cycles % 1.0, xs_ext, ys_ext)


def _residual_stats(ours: np.ndarray, ref: np.ndarray):
    """Peak-normalize both curves and return (ours_n, ref_n, max|Δ|, RMS Δ)."""
    ours_n = ours / ours.max()
    ref_n = ref / ref.max()
    resid = ours_n - ref_n
    return ours_n, ref_n, float(np.abs(resid).max()), float(np.sqrt((resid ** 2).mean()))


def _asymmetry(flux_n: np.ndarray) -> float:
    """Fore-aft asymmetry: RMS of f(φ) − f(−φ) over the cycle (0 for a static pulse).

    The pure Doppler signature. Computed by mirroring the peak-normalized curve about
    phase 0 (index k → −k on the periodic grid) and differencing; a frozen, symmetric
    pulse gives ~0, a rotating one a finite value that grows with spin.
    """
    mirrored = flux_n[(-np.arange(flux_n.size)) % flux_n.size]
    return float(np.sqrt(np.mean((flux_n - mirrored) ** 2)))


def main():
    for label, spin, path in [("SD1a (1 Hz)", SPIN_SD1A_HZ, REF_SD1A),
                              ("SD1c (200 Hz)", SPIN_SD1C_HZ, REF_SD1C)]:
        ref = np.loadtxt(path)
        phase, f_ref = ref[:, 0], ref[:, 1]
        beta = spot_speed(spin, NS_RADIUS_KM, COLATITUDE, COMPACTNESS)
        print(f"{label}:  β = {beta:.5f}")
        for tag, wd in [("Doppler+aberration+delay", True), ("delay omitted (prod. regime)", False)]:
            ours_n, ref_n, max_abs, rms = _residual_stats(
                monochromatic_flux(phase, spin, with_time_delay=wd), f_ref)
            print(f"    {tag:<30}  max|Δ| = {100*max_abs:5.3f}%   RMS = {100*rms:5.3f}%   "
                  f"asym: ours {100*_asymmetry(ours_n):5.2f}% / ref {100*_asymmetry(ref_n):5.2f}%")

    # Figure: overlay both spins + the SD1c residual.
    ref_a = np.loadtxt(REF_SD1A)
    ref_c = np.loadtxt(REF_SD1C)
    phase = ref_c[:, 0]
    a_ours, a_ref, _, _ = _residual_stats(monochromatic_flux(phase, SPIN_SD1A_HZ), ref_a[:, 1])
    c_ours, c_ref, c_max, c_rms = _residual_stats(monochromatic_flux(phase, SPIN_SD1C_HZ), ref_c[:, 1])

    fig, (ax, ax_r) = plt.subplots(
        2, 1, figsize=(8, 6.8), height_ratios=[3, 1], sharex=True, layout="constrained")

    ax.plot(phase, a_ref, "-", lw=1.6, color="#7f7f7f", label="IM reference — SD1a (1 Hz)")
    ax.plot(phase, c_ref, "-", lw=2.0, color="#1f77b4", label="IM reference — SD1c (200 Hz)")
    ax.plot(phase[::3], c_ours[::3], "o", ms=4.2, color="#d62728",
            label="this work — SD1c (Doppler + aberration)")
    ax.set_ylabel("normalized flux  f(φ) / f_max")
    ax.set_title("SD1c Doppler validation: 200 Hz point spot vs. NICER code comparison "
                 f"(u = {COMPACTNESS:.3f}, i = θ_s = 90°)")
    ax.legend(loc="upper center", fontsize=9)
    ax.margins(x=0.01)

    ax_r.axhline(0.0, color="#bbbbbb", lw=0.8)
    ax_r.plot(phase, 100.0 * (c_ours - c_ref), "-", color="#555555", lw=1)
    ax_r.set_xlabel("rotational phase  φ / 2π  [cycles]")
    ax_r.set_ylabel("Δ  [% of peak]")
    ax_r.set_title(f"SD1c: this work − reference   "
                   f"(max |Δ| = {100*c_max:.2f}%, RMS = {100*c_rms:.2f}% — worst at eclipse edge)")

    fig.savefig(FIGURE_PATH, dpi=130)
    print(f"saved {FIGURE_PATH}")


if __name__ == "__main__":
    main()
