# Updated paper figures (APSS drafts)

Redesigns of the three load-bearing figures of the manuscript, drafted for the
Astrophysics and Space Science (APSS) submission. Regenerate with:

```bash
python make_updated_figs.py    # writes fig_*_updated.{pdf,png}
```

Everything is recomputed from the repo's own pipeline (`mcrt` +
`scripts/anchor_lib.py` + the committed `data/*.npz` evidence base) — no
hand-copied numbers. PDFs are vector (fonts embedded, TrueType); PNGs are
600 dpi previews. Sized for the Springer sn-jnl text block: full width
170 mm, all text ≥ ~7 pt at print size.

## Why redesign: audit of the current figures

| Current figure | Journal-readiness problems |
|---|---|
| `pulse_profile_j0740.png` / `j0030.png` | Conversational in-figure titles that editorialize ("both teams stay positive; not a saturation artifact"); literal underscores (`PF_iso=0.17`); three panels whose content is **not described by the one-line caption in the .tex** (caption promises "pulse profiles", figure also shows ΔPF(τ) and background robustness) |
| `doppler_dilution_ladder.png` | Internal track codenames leak ("C2 → C3 → C4"); `PF_real − PF_iso` with raw underscores; jargon rung names ("time-delay warp") |
| `phase_diagram.png` | Two-line conversational title duplicating the caption; anchor labels collide; regime names only in the title, not on the map |
| `d_isotropic_h_overlay.png` | Internal codename in title ("Track D") |
| `pulse_profile_doppler_sd1c.png` | **Title clipped at the right edge** ("i = θ_s = 9…"); raw underscores (`f_max`) |
| all | In-figure titles (journals put the message in the caption); default matplotlib color cycle (red/green issues aside, red–blue chosen is fine but not colorblind-optimal); no panel letters |

## Design contract used here

One message per figure; the caption states it. No in-figure titles, no
internal codenames, no raw underscores. Okabe–Ito colorblind-safe palette
with one consistent encoding across every figure:

- **blue #0072B2** = Riley fits, **vermillion #D55E00** = Miller fits,
  **gray dashed** = isotropic reference, **sky blue #56B4E9** = analytic
  boundary;
- panel letters (a)/(b)/(c); direct labeling preferred over legends;
- ticks inward on all four sides; minimal grid; vector output.

Sources (verify/browse):

- Rougier, Droettboom & Bourne (2014), "Ten Simple Rules for Better
  Figures", PLoS Comput Biol 10(9):e1003833 —
  <https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1003833>
- Okabe & Ito colorblind-safe palette — <https://jfly.uni-koeln.de/color/>
- Springer artwork instructions (line art 600–1200 dpi, halftones 300 dpi,
  EPS/TIFF preferred, standard fonts, color free of charge) — APSS
  submission guidelines: <https://link.springer.com/journal/10509/submission-guidelines>

## Revealed APSS figure norms (survey of 12 papers, 2023–2026)

A survey of recent APSS papers (full per-paper notes and URLs in the
research report accompanying this folder) found:

1. 6–12 main figures is the comfortable norm for a numerical paper;
   appendix figure dumps (20–45 total) are tolerated, not penalized.
2. Clean matplotlib output is the *revealed standard* — roughly half of
   published figures are default-styled. Consistent line-style/legend
   systems mark the stronger papers; custom design is not expected.
3. Parameter sweeps appear as color-coded curve families with explicit
   legends; parameter-space results as 2D heatmaps with colorbars
   (log-scaled where the quantity spans decades).
4. Validation is shown by overlaying observational constraints/reference
   curves directly on model plots; separate residual sub-panels are rare
   (ours are a strength worth keeping, not a liability).
5. Caption style splits between telegraphic and 100–200-word
   self-contained; the stronger papers use the latter. Panels are cited as
   Top/Middle/Bottom or Left/Right or (a)/(b).
6. Tables carry setup parameters and summary numbers (0–3 per paper);
   figures carry the results.

Closest published analogs (worth skimming for tone): a Be-star Monte Carlo
radiative-transfer prediction-grid paper (APSS 370, DOI
10.1007/s10509-025-04512-w), a dark-matter-admixed strange-star structure
paper (arXiv:2501.06379), and an SPH + MC radiative-transfer binaries
paper (DOI 10.1007/s10509-025-04494-9).

## The three drafts

### `fig_j0740_updated` — replaces Fig. 3 (`pulse_profile_j0740.png`)

(a)/(b) peak-normalized pulse profiles for the two published J0740 fits,
isotropic (gray dashed) vs. realistic beaming (color), shaded difference,
ΔPF stamped per panel; (c) ΔPF vs. slab optical depth with seed-scatter
error bars. Exact Schwarzschild bending, τ = 10 library — the paper's
production configuration. The background-robustness panel of the old
figure is dropped (it is one sentence in Sect. 3.1).

Suggested caption: *"Isotropic-emission bias for the two published
PSR J0740+6620 geometries in the static, bolometric limit. Panels (a) and
(b): peak-normalized two-spot pulse profiles for the Riley et al. (2021)
and Miller et al. (2021) fits, computed with exact Schwarzschild bending
and the τ = 10 Monte Carlo beaming library (solid) and with isotropic
surface emission (dashed); the shaded region is the waveform difference.
Realistic beaming deepens the modulation, raising the pulsed fraction by
ΔPF = +0.137 (Riley) and +0.195 (Miller). Panel (c): ΔPF versus slab
optical depth τ; error bars show the seed-to-seed scatter of the
atmosphere library. The bias saturates for τ ≳ 3, confirming the τ = 10
library samples the optically thick limit."*

### `fig_routing_updated` — replaces Fig. 6 (`doppler_dilution_ladder.png`)

The paper's central claim in one figure. (a) ΔPF dilution ladder with
physically-named rungs matching Table 2 row-for-row (values printed beside
each marker); (b) the waveform-shape rms on the same rungs — flat at
≈ 0.10–0.12 while ΔPF collapses from +0.14/+0.20 to +0.02/+0.07. A
skimming reader gets "spin kills the PF bias but not the shape bias"
without reading the text.

Suggested caption: *"Rotation redistributes the isotropic-emission bias
between observables for PSR J0740+6620. Panel (a): pulsed-fraction bias
ΔPF as successive physical effects are added to the static baseline (top
to bottom; rungs correspond to the rows of Table 2), for the Riley et al.
(2021) (circles) and Miller et al. (2021) (squares) geometries. Error bars
show Monte Carlo seed scatter; light-travel-delay rungs carry the
convention spread of Sect. 3.4. Panel (b): the peak-normalized
waveform-shape rms difference on the same rungs. While Doppler harmonic
structure dilutes ΔPF by a factor of ~4–7, the waveform-shape difference
is nearly unchanged at ≈ 0.10 rms: the bias moves, it does not vanish."*

### `fig_phase_diagram_updated` — replaces Fig. 5 (`phase_diagram.png`)

Same magma ΔPF maps, with the regimes labeled **on the map** ("PF
saturated — bias moves to waveform shape" / "PF sensitive — bias visible
in ΔPF"), the analytic tiling boundary in sky blue over the numerical
eclipse contour (white dashed), published anchor configurations as stars,
and a shared colorbar. Anchor labels de-collided; x-axis clipped to the
physical Δφ ∈ [0, 0.5] range (the stored boundary uses 1.0 as a
"no-boundary" sentinel, now masked).

Suggested caption: *"Where the isotropic-emission bias appears, as a
function of hot-spot geometry. Each map shows the static pulsed-fraction
bias ΔPF (τ = 10 library) versus azimuthal spot separation Δφ and spot
colatitude θ_s, for (a) a J0740-like star (u = 0.49, i = 87.6°) and (b) a
J0030-like star (u = 0.31, i = 53.9°). Black regions are saturated:
the summed two-spot flux reaches zero at some phase, PF is pinned at
unity for any beaming, and the bias is expressed only in the waveform
shape. The sky-blue curve is the analytic tiling boundary derived in
Sect. 3.3 (linear bending); the white dashed contour is the numerical
eclipse onset from the full engine. Stars mark the published fits: both
J0740 solutions tile the rotation (PF-sensitive regime), both J0030
solutions do not (saturated regime)."*

## Not yet redrafted (recommended next)

- **Validation gallery (Fig. 2)**: fix "Track D" title, clipped SD1c
  title, and underscores; consider one 3-panel figure with shared style —
  the residual sub-panels are already good practice.
- **Geometry schematic (Fig. 1)**: fine as is (TikZ, vector).
- **J0030 counterpart (Fig. 4)**: same treatment as `fig_j0740_updated`
  (profiles + shape-rms vs. τ instead of ΔPF).
- Regenerate final submission versions at 600 dpi from this script's PDFs
  (vector already satisfies Springer line-art requirements).
