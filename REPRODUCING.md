# Reproducing the paper

Every number and figure in the paper ("The Systematic in the isotropic assumption
for Neutron stars", `docs/paper/paper_template.tex`) is produced by a script in
`scripts/` and saved to an artifact in `data/`. This file gives the execution
order (script → artifact → paper location). For the reverse direction — "which
file backs this specific claim?" — see `docs/paper/claims-evidence.md`.

## Setup

```bash
pip install -r requirements.txt   # pinned: bit-for-bit MC repro needs these exact versions
```

All scripts run from the repository root with:

```bash
PYTHONPATH=src:scripts python3 scripts/<script>.py
```

Each script's module docstring states what it computes, what it writes, and any
`--quick` plumbing-check mode. The Monte Carlo pipeline is seeded end-to-end, so
reruns reproduce the shipped npz values exactly (same numpy version, see
`requirements.txt`).

## Third-party reference data (not redistributed)

The two code-comparison validations (`code_comparison.py`, `c1_doppler_validate.py`)
test against reference waveforms from Bogdanov et al. (2019, ApJL 887, L26) — AAS
supplementary material whose redistribution is not granted, so it is **not** in this
repository. To obtain it: download the supplementary archive `apjlab5968.tar.gz`
from the article page (doi:10.3847/2041-8213/ab5968) and extract
`SD1a_test_IM.txt` and `SD1c_test_IM.txt` into `data/l26_reference/`.
Everything else runs without it.

## Execution order

**Step 0 — the Monte Carlo beaming library** (everything downstream reads it):

| Script | Writes | Used for |
|---|---|---|
| `tau_sweep.py` | `data/beaming_library.npz` | the production I(μ; τ) library: 6 τ values × 5 seeds, ~4×10⁵ escaped photons each (§2.1) |

**Step 1 — validation** (§2.3, Figure 2):

| Script | Writes | Paper |
|---|---|---|
| `d_isotropic_validate.py` | `d_isotropic_h.npz`, `d_isotropic_h_overlay.png` | H(μ) flux-space χ² = 0.70 (Fig. 2 left) |
| `code_comparison.py` † | `code_comparison.npz`, `pulse_profile_code_comparison.png` | SD1a: 0.11% exact / 0.81% linear (Fig. 2 middle) |
| `c1_doppler_validate.py` † | `c1_doppler_validate.npz`, `pulse_profile_doppler_sd1c.png` | SD1c: 1.35% at 200 Hz (Fig. 2 right) |

† needs `data/l26_reference/` (see above).

**Step 2 — the two anchors** (§2.4, §3.1–3.2, Figures 3–4):

| Script | Writes | Paper |
|---|---|---|
| `j0740_anchor.py` | `j0740_anchor.npz`, `pulse_profile_j0740.png` | static ΔPF +0.137/+0.195; tiling floors 0.71/0.63; ~14%/~23% eclipses |
| `j0030_anchor.py` | `j0030_anchor.npz`, `pulse_profile_j0030.png` | PF saturation; shape RMS 0.061/0.058; background re-exposure +0.03…+0.14 |

**Step 3 — error bars and bending** (§3.1, §2.2):

| Script | Writes | Paper |
|---|---|---|
| `a3_seed_errors.py` | `a3_seed_errors.npz` | seed-to-seed errors ±0.003/±0.005 |
| `b2_exact_bending.py` | `b2_exact_bending.npz` | linear-map bias −0.0066 at u = 0.494 |

**Step 4 — spin, band, and geometry dilution** (§3.4, Figure 6, Table 2):

| Script | Writes | Paper |
|---|---|---|
| `c2_doppler_coupling.py` | `c2_doppler_coupling.npz` | ladder rows 1–2 (δ⁴); Δ(A1/A0) +0.033→+0.029 |
| `c3_band_doppler.py` | `c3_band_doppler.npz` | ladder row 3 (NICER-band photons) |
| `c4_caveat_audit.py` | `c4_caveat_audit.npz` | ladder rows 4–5 (light-travel warp); oblateness bound; warp-inclusive shape RMS ≈ 0.10 |
| `doppler_figures.py` | `doppler_dilution_ladder.png` (+2 more) | Figure 6 |

**Step 5 — the phase diagram / tiling criterion** (§3.3, Figure 5):

| Script | Writes | Paper |
|---|---|---|
| `phase_diagram.py` | `phase_diagram.npz`, `phase_diagram.png` | the regime map, Figure 5 |
| `validate_phase_diagram.py` | `validate_phase_diagram.png` | engine eclipse onset ≡ analytic condition (0 mismatches) |

**Step 6 — robustness checks** (§2.1, §4):

| Script | Writes | Paper |
|---|---|---|
| `e2_tail_sensitivity.py` | `e2_tail_sensitivity.npz`, `e2_tail_sensitivity.png` | grazing-tail perturbations shift ΔPF ≤ 0.006 |
| `e3_atmosphere_laws.py` | `e3_atmosphere_laws.npz`, `e3_atmosphere_laws.png` | ΔPF tracks b across Eddington/H(μ)/Thomson |
| `finite_cap.py` | `finite_cap.npz`, `finite_cap.png` | point-spot reduction bias ≤ 0.004 (linear map) |
| `convergence_study_v2.py` | `convergence_v2_results.npz` (+2 png) | N^(−1/2) convergence, no systematic floor (long run; `--quick` for plumbing) |

**Figure 1** (`docs/paper/geometry_figure.pdf`) is a native TikZ redraw — build
`docs/paper/geometry_figure.tex` with pdflatex.

## Tests

```bash
PYTHONPATH=src:scripts python3 -m pytest tests/   # 127 tests
```

`tests/test_rotating.py` includes the from-scratch 4-vector Lorentz cross-check
quoted in §2.3 (machine-precision agreement).

## Building the paper

```bash
cd docs/paper
pdflatex paper_template.tex && bibtex paper_template && pdflatex paper_template.tex && pdflatex paper_template.tex
```

Figures are pulled from `data/` via `\graphicspath`.
