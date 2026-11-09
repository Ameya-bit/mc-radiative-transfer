# Reproducing the paper

Every number and figure in the paper is produced by a script in `scripts/` and
saved to an artifact in `data/`. This file gives the execution order (script →
artifact → paper location). For the reverse direction — "which file backs this
specific claim?" — see `docs/paper/claims-evidence.md`.

**Manuscript sources are not in this repository.** The LaTeX, bibliography,
and submission material for the *Astrophysics and Space Science* version
(`docs/paper/apss_template.tex`, Springer `sn-jnl` class) and the earlier
AASTeX drafting copy (`docs/paper/paper_template.tex`) are kept privately and
ignored by git. What is tracked is everything needed to regenerate the numbers
and figures: the scripts, the `data/*.npz` artifacts, the figure generator in
`docs/paper/poster/updated_graphs/`, and the paper-prep documentation
(`docs/paper/*.md`). Figures for the APSS version are the redesigned set in
`docs/paper/poster/updated_graphs/` (see that folder's `README.md`); figures
for the AASTeX version are the originals in `data/`.

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

## Revision round (APSS submission)

These were added while preparing `docs/paper/apss_template.tex` in response to
a pre-submission referee-style review. Run **after** Steps 0–6, since the script
reads the artifacts they produce.

| Script | Writes | Paper |
|---|---|---|
| `r1_revision_checks.py` | `r1_revision_checks.npz` | six checks, listed below |
| `nicer_counts.py` | reads `nicer_counts_j0030.csv`, `nicer_counts_j0740.csv` | adopted NICER count levels (§2.4): background / spot counts per bin digitised from Miller et al. 2019 Fig. 10 and 2021 Fig. 6 |

```bash
PYTHONPATH=src:scripts python3 scripts/r1_revision_checks.py
```

Reads `beaming_library.npz`, `a3_seed_errors.npz`, `e2_tail_sensitivity.npz`,
`b2_exact_bending.npz`, `finite_cap.npz`. Runtime a few minutes (R1's grid
rebuilds an `ExactBending` map per trial compactness).

| Check | Answers | Result used in the paper |
|---|---|---|
| **R1** toy compactness inversion | "what is the bias worth at the inference level?" — fits the isotropic model to the beamed J0740 profile with only `u` free | Δu ≈ −0.13 both fits ⇒ radius overestimated by **+4.6 km** (Riley) / **+5.9 km** (Miller); rms mismatch drops 0.107→0.039 and 0.125→0.035 (§3.1, §5, abstract) |
| **R2** J0030 at 205.53 Hz | the spinning J0030 numbers the abstract previously implied but never computed | ΔPF stays 0.0000 (saturation is geometric); shape rms 0.061→**0.060** (Riley), 0.058→**0.057** (Miller) (§3.2, abstract) |
| **R3** slope sensitivity | whether the b = 1.72 vs 1.79 discrepancy matters | ΔPF differs by only **0.004** between the two linear laws (§2.3, replaces the old footnote) |
| **R4** exact-bending boundary | whether the tiling criterion's placement survives exact bending across the whole grid | boundary shifts by ≤ **0.15 cyc** (J0740-like) / **0.10 cyc** (J0030-like); no anchor changes regime (§3.3, Figure 5's second boundary curve) |
| **R5** error budget | separates the statistical bar from the systematics | seed ±0.0032/±0.0052; tail 0.0059; bending −0.0066/−0.0060; point-spot +0.0036/+0.0017 ⇒ total systematic ≲0.01 (§3.1) |
| **R6** Vinciguerra+2024 J0030 modes | whether the criterion predicts the regime of the *reanalyzed* J0030 geometry | both favored modes **tile** (floor 0.60 / 0.52) ⇒ PF-sensitive regime with ΔPF **+0.24 / +0.25** — the same star crosses the boundary, as the criterion predicts (§3.3, §5) |

R6's input geometry is hard-coded in the script from Vinciguerra et al. (2024,
ApJ 961, 62) online **Figure set 6** — the hot-region medians are *not* printed
in that paper's tables. Compactness is formed as `u = 2·1.47662·(M/M☉)/R_km`
and the separation as `Δφ = (φ_s + 0.5) − φ_p`, the half-cycle following their
Table 1 definition of the secondary phase. Because both favored models are of
protruding double-temperature type (a small hot cap offset inside a larger cool
one), the check uses the hot component as a point-spot stand-in and is
indicative for the ΔPF value, though the regime assignment follows from
geometry alone.

## Figures for the APSS submission

```bash
python3 docs/paper/poster/updated_graphs/make_updated_figs.py
```

Regenerates all five plotted paper figures as vector PDF plus 600 dpi PNG,
reading the same `data/*.npz` artifacts (including `r1_revision_checks.npz` for
Figure 5's exact-bending boundary). Rationale, the figure audit, and the
surveyed APSS figure conventions are in that folder's `README.md`.

## Tests

```bash
PYTHONPATH=src:scripts python3 -m pytest tests/   # 127 tests
```

`tests/test_rotating.py` includes the from-scratch 4-vector Lorentz cross-check
quoted in §2.3 (machine-precision agreement).

## Building the paper

The commands below apply to a checkout that has the private manuscript sources
in `docs/paper/` (they are git-ignored; see the note at the top of this file).

AASTeX drafting copy:

```bash
cd docs/paper
pdflatex paper_template.tex && bibtex paper_template && pdflatex paper_template.tex && pdflatex paper_template.tex
```

APSS submission version (needs `sn-jnl.cls` and `sn-basic.bst` alongside it):

```bash
cd docs/paper
pdflatex apss_template.tex && bibtex apss_template && pdflatex apss_template.tex && pdflatex apss_template.tex
```

Figures are pulled from `data/` and `docs/paper/poster/updated_graphs/` via
`\graphicspath`. `apss_template.tex` renders revision changes as markup —
additions in red, deletions struck through — controlled by one switch near the
top of the preamble:

```latex
\showmarkuptrue    % review copy: additions red, deletions struck out
\showmarkupfalse   % clean copy: additions plain black, deletions removed
```

Flip it to `\showmarkupfalse` to produce the version for submission.
