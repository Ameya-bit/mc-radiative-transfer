# References & Data Sources

> The project's single bibliography. Deep dives and the README link **into** this
> file by anchor (e.g. `[Beloborodov (2002)](../references.md#beloborodov-2002)`)
> instead of repeating full citations. Add a source here once, cite it everywhere.
>
> The pulse-profile / NICER bibcodes below were checked against ADS/IOP on
> 2026-06-09 — all five confirmed (Beloborodov 2002 ApJ 566 L85; Poutanen &
> Beloborodov 2006 MNRAS 373 836; Bogdanov 2019 ApJL 887 L26; Riley 2019 ApJL 887
> L21; Miller 2019 ApJL 887 L24). When you add a source, include a stable link
> (ADS bibcode or DOI) and note which version first used it.

---

## Monte Carlo radiative transfer (the engine)

### Cashwell & Everett (1959)
*A Practical Manual on the Monte Carlo Method for Random Walk Problems.* Pergamon.
— The single-photon random-walk method our engine actually implements.
*Used in:* v0.1.0–v0.2.0 (transport), and the paper's methods framing.

### Whitney (2011)
*Monte Carlo radiative transfer.* Bulletin of the Astronomical Society of India, 39, 101.
[ADS: 2011BASI...39..101W](https://ui.adsabs.harvard.edu/abs/2011BASI...39..101W) —
modern review of the method; the correct citation in place of "Lucy 1999 photon
packets" (we do single-photon transport, not Lucy packet splitting).
*Used in:* v0.2.0 onward (method), paper methods.

---

## Classical radiative transfer & limb darkening (the validation targets)

### Chandrasekhar (1960)
*Radiative Transfer.* Dover. — The H-function for a semi-infinite scattering
atmosphere; the exact limb-darkening law our beaming function is validated
against. Two caveats to state honestly in the paper:
(1) **scattering law** — the engine uses the Rayleigh/Thomson pattern ¾(1+μ²),
while the H compared against is the *scalar isotropic* (p=1) solution
(Chandrasekhar's Rayleigh case is a polarization-coupled matrix problem); the
shapes agree closely but are not strictly the same physics. A clean apples-to-
apples fix is an isotropic-scattering mode that should reproduce scalar H(μ) to
within noise.
(2) **finite vs. semi-infinite** — H is the τ→∞ limit, so finite slabs (τ≲1)
deviate; the rigorous finite-τ benchmark is Chandrasekhar's **X- and Y-functions**.
*Used in:* v0.5.1, v0.6.0–v0.6.1 (`src/mcrt/theory.py`, beaming validation).

### Eddington limb darkening
The classical `I(μ) ∝ 1 + 1.5 μ` two-stream result (textbook; see Chandrasekhar
1960 / Mihalas 1978). The linear approximation our fitted slope `b` is compared to.
*Used in:* v0.5.1 onward (`eddington_limb_darkening`).

---

## Pulse profiles & gravitational light bending (the pulse module)

### Beloborodov (2002)
ApJ, 566, L85. [DOI: 10.1086/339511](https://doi.org/10.1086/339511) ·
[ADS: 2002ApJ...566L..85B](https://ui.adsabs.harvard.edu/abs/2002ApJ...566L..85B)
— The linear light-bending approximation `cos α = u + (1 − u) cos ψ` with constant
Jacobian `(1 − u)`; accurate to ~1% for `u ≲ 0.5`. The analytic-check closed form
and the core of `bend()`.
*Used in:* v0.8.0 (`bend`, `analytic_isotropic_pf`), v0.8.1 (code-comparison check).

### Poutanen & Beloborodov (2006)
MNRAS, 373, 836. [DOI: 10.1111/j.1365-2966.2006.11088.x](https://doi.org/10.1111/j.1365-2966.2006.11088.x) ·
[ADS: 2006MNRAS.373..836P](https://ui.adsabs.harvard.edu/abs/2006MNRAS.373..836P)
— The bolometric point-spot flux `F ∝ (1 − u) I(cos α) cos α` (slow-rotation form).
*Used in:* v0.8.0 (`point_spot_flux`).

---

## NICER pulse-profile modeling & code comparison (the benchmarks)

### Bogdanov et al. (2019) — "Paper II", L26
ApJL, 887, L26. [DOI: 10.3847/2041-8213/ab5968](https://doi.org/10.3847/2041-8213/ab5968) ·
[ADS: 2019ApJ...887L..26B](https://ui.adsabs.harvard.edu/abs/2019ApJ...887L..26B)
— The pulse-profile **code-comparison** suite (test problems SD1/SD2/OS1). The
Illinois–Maryland (IM) reference waveforms are our **code-comparison** benchmark;
codes agree to ≲ 0.1%. **Note:** this is Paper *II* (L26), not Paper I (L25, the data set).
The reference data is a separate, gitignored download — see
[Reference data sets](#reference-data-sets) below.
*Used in:* v0.8.1 (code-comparison check, test SD1a).

### Riley et al. (2019) — J0030, L21
ApJL, 887, L21. [DOI: 10.3847/2041-8213/ab481c](https://doi.org/10.3847/2041-8213/ab481c) ·
[ADS: 2019ApJ...887L..21R](https://ui.adsabs.harvard.edu/abs/2019ApJ...887L..21R)
— Amsterdam (X-PSI) M–R and spot geometry for PSR J0030+0451. **Values used (preferred
ST+PST model, Table 2):** M = 1.34 M⊙, R_eq = 12.71 km (compactness GM/Rc² = 0.156 ⇒
u = 0.312); inclination i = 0.94 rad (53.9°); ST spot colatitude Θp = 2.23 rad (127.8°),
PST crescent Θs = 2.91 rad (166.7°); both in the same hemisphere, azimuthal separation
≈ 0.45 cyc. The two spots are equal-temperature (log₁₀T = 6.11).
*Used in:* v0.9.1 (real-star anchor — `scripts/j0030_anchor.py`).

### Miller et al. (2019) — J0030, L24
ApJL, 887, L24. [DOI: 10.3847/2041-8213/ab50c5](https://doi.org/10.3847/2041-8213/ab50c5) ·
[ADS: 2019ApJ...887L..24M](https://ui.adsabs.harvard.edu/abs/2019ApJ...887L..24M)
— Illinois–Maryland M–R for PSR J0030+0451 and quoted uncertainties. **Values used
(three-oval-spot model, Table 8 medians):** M = 1.443 M⊙, R_e = 13.019 km (compactness
GM/Rc² = 0.163 ⇒ u = 0.326); observer θobs = 0.878 rad (50.3°); two main spots at
colatitudes θc1 = 2.270 rad (130.0°), θc2 = 2.417 rad (138.5°), longitude separation
Δφ2 = 0.460 cyc, oval sizes Δθ·f and temperatures kT ≈ 0.115–0.117 keV (for the
area × T⁴ weighting); the third spot is negligible and dropped.
*Used in:* v0.9.1 (real-star anchor — `scripts/j0030_anchor.py`).

### Riley et al. (2021) — J0740+6620, L27
ApJL, 918, L27. [DOI: 10.3847/2041-8213/ac0a81](https://doi.org/10.3847/2041-8213/ac0a81) ·
[ADS: 2021ApJ...918L..27R](https://ui.adsabs.harvard.edu/abs/2021ApJ...918L..27R) ·
arXiv:2105.06980 — Amsterdam (X-PSI) M–R and spot geometry for the high-mass MSP
PSR J0740+6620. **Values used (preferred ST-U model, two single-temperature spots,
headline NICER × XMM run, Table 4 medians + CI68):** M = 2.072₋₀.₀₆₆⁺⁰·⁰⁶⁷ M⊙,
R_eq = 12.39₋₀.₉₈⁺¹·³⁰ km ⇒ compactness GM/Rc² = 0.247, u = 2GM/Rc² = 0.494 (very
compact → strong light bending). Observer: cos i = 0.0424 ⇒ i ≈ 1.528 rad (87.57°,
almost edge-on). Spot p (primary): Θp = 1.35₋₀.₃₉⁺⁰·⁴⁶ rad (77.3°), angular radius
ζp = 0.147 rad, log₁₀T = 5.988. Spot s (secondary): Θs = 1.89₋₀.₄₆⁺⁰·⁴⁰ rad (108.3°),
angular radius ζs = 0.146 rad, log₁₀T = 5.992 (abstract rounds both to log₁₀T = 5.99).
Two spots, both retained, near-equal small caps; phases bimodal. Strong ~4-fold
geometric degeneracy because i is near the equator and the spots do not overlap.
**Eclipse check (cos(i+Θ) ≥ −u/(1−u) = −0.976): both spots ALWAYS VISIBLE** —
p: i+Θp = 164.9°, cos = −0.966; s: i+Θs = 195.9°, cos = −0.962. Neither spot is
ever eclipsed; high compactness lets the observer see around the back. **This is the
non-eclipsing geometry.**
*Used in:* v0.9.2 (second real-star anchor — primary, non-eclipsing fit;
`scripts/j0740_anchor.py`).

### Miller et al. (2021) — J0740+6620, L28
ApJL, 918, L28. [DOI: 10.3847/2041-8213/ac089b](https://doi.org/10.3847/2041-8213/ac089b) ·
[ADS: 2021ApJ...918L..28M](https://ui.adsabs.harvard.edu/abs/2021ApJ...918L..28M) ·
arXiv:2105.06979 — Illinois–Maryland independent M–R and spot geometry for
PSR J0740+6620; the second-team cross-check of Riley 2021. **Values used (two uniform
circular spots, headline NICER + XMM run, Table 7 medians ± 1σ):** M = 2.062₋₀.₀₉₁⁺⁰·⁰⁹⁰ M⊙,
R_e = 13.713₋₁.₅₀₄⁺²·⁶¹³ km ⇒ GM/Rc² = 0.222, u = 0.444. θobs = 1.527 rad (87.49°).
Spot 1: θc1 = 1.600 rad (91.67°), Δθ1 = 0.098 rad, kT1 = 0.094 keV. Spot 2:
θc2 = 1.612 rad (92.36°), Δθ2 = 0.096 rad, kT2 = 0.094 keV (equal temp). Phase
separation Δφ2 = 0.558 cyc (anti-phased). Agrees with Riley on M, R, i, anti-phasing,
and the two-circular-spot model, but places the spots essentially **on the equator**
(Θ ≈ 92°), the least-constrained parameter. **Eclipse note:** with equatorial spots and
an equatorial observer each spot's *center* is hidden ~21% of the cycle — yet because
the two are anti-phased the *combined* pulse never reaches zero, so PF stays unsaturated
and ΔPF is live (≈ +0.23). The takeaway: what saturates PF is not single-spot eclipse
but whether the two spots fail to *tile* the rotation (as J0030's same-hemisphere spots
do, and J0740's anti-phased spots do not).
*Used in:* v0.9.2 (second-team cross-check; `scripts/j0740_anchor.py`).

### Choudhury et al. (2024) — J0437−4715, L20
ApJL, 971, L20. [DOI: 10.3847/2041-8213/ad5a6f](https://doi.org/10.3847/2041-8213/ad5a6f) ·
[ADS: 2024ApJ...971L..20C](https://ui.adsabs.harvard.edu/abs/2024ApJ...971L..20C) ·
arXiv:2407.06789 — Amsterdam (X-PSI) M–R and spot geometry for the nearest/brightest
MSP PSR J0437−4715. **Values used (preferred CST+PDT model, Table 1 medians + CI68):**
M = 1.418 ± 0.037 M⊙ (radio-prior dominated), R_eq = 11.36₋₀.₆₃⁺⁰·⁹⁵ km ⇒ compactness
GM/Rc² = 0.184, u = 0.369. Observer: cos i = −0.7373 ⇒ i = 2.39993 rad (137.506°,
viewed from the southern hemisphere). Two hot regions, non-antipodal:
- **Primary p** — single-temperature ring around the north pole: center colatitude
  Θp = 0.146₋₀.₀₂₂⁺⁰·⁰²⁸ rad (8.37°); emitting-cap radius ζp = 0.433 rad with a
  concentric masking ("omitting") sub-cap ζ_o,p = 0.139 rad; log₁₀Tp = 6.101; phase
  ϕp = 0.4429 cyc.
- **Secondary s** — dual-temperature ("protruding") spot between equator and south
  pole: superseding component center Θs = 2.348₋₀.₀₇₀⁺⁰·⁰⁶⁰ rad (134.5°), radius
  ζs = 0.0299 rad, log₁₀Ts = 6.219; ceding component center Θc,s = 2.307 rad (132.2°),
  radius ζc,s = 0.197 rad, log₁₀Tc,s = 5.752; phase ϕs = 0.4704 cyc; component
  azimuthal offset χs = 0.0041 rad. Spots are nearly anti-phased (Δϕ ≈ 0.028 cyc).
**Eclipse check (cos(i+Θ) ≥ −u/(1−u) = −0.584):** primary p **ECLIPSES**
(i+Θp = 145.9°, cos = −0.828 < −0.584 → goes behind the star each rotation); the
secondary components stay marginally visible (i+Θ ≈ 270°, cos ≈ 0). So J0437 is **not**
the all-visible case — its polar primary disappears behind the limb.
*Used in:* evaluated as a candidate for v0.9.2 but **not used** — it eclipses, so it is
not the non-eclipsing star we wanted (J0740 is).

> **Eclipse summary for the pulse-profile model.** Using the point-spot criterion
> "always visible if cos(i + Θ_spot) ≥ −u/(1−u)" (light bending lets you see slightly
> around the far side): **J0740+6620 is the non-eclipsing star — under Riley's fit BOTH
> spots remain visible for the entire rotation** (near-edge-on i ≈ 87.6° but extreme
> compactness u ≈ 0.494). **J0437−4715 DOES eclipse**: its near-polar primary spot
> (Θ ≈ 8°, i ≈ 137.5°) passes behind the star each cycle. Pick **J0740+6620** when you
> need a geometry whose hot spots never go behind the star.
>
> **Caveat on the single-spot criterion (learned in v0.9.2).** The criterion above is
> about *one* spot. What actually saturates the pulsed fraction is whether the *summed*
> multi-spot pulse touches zero, i.e. whether the spots **tile** the rotation. J0030's
> two spots hug the same far hemisphere and each is hidden ~45% of the cycle, so even
> anti-phased they leave a dark gap → flux hits zero → PF pins at 1. J0740's two spots
> are anti-phased in *opposite* regions, so even Miller's median (whose spot centers do
> eclipse ~21% each) tiles the cycle and keeps the combined flux off zero → PF stays
> unsaturated. Tiling, not single-spot eclipse, is the discriminator.

---

## Beaming systematics & pulse-profile regime classification (paper novelty positioning)

> Surfaced in the v0.9.2 → v1.0.0 literature review (2026-06-25), confirmed by full-text reads.
> These are the prior-work citations the paper must engage: the canonical maps we **validate
> against**, the closest prior work we must **distinguish from**, and the qualitative antecedents
> of the effect. Each note records *why it is here* so the novelty framing is not lost at drafting.
>
> **One-line novelty summary.** No prior paper (a) reports a single ΔPF number for the
> isotropic→realistic beaming swap, (b) maps where that systematic lands (pulsed-fraction
> amplitude vs waveform shape) over geometry — in particular over **azimuthal spot separation**
> — or (c) states the **tiling** criterion (combined two-spot flux touches zero ⇔ PF saturates
> at 1, set by whether anti-phased spots cover the rotation, *not* by single-spot eclipse).

### Zhao, Psaltis, Özel & Beklen (2024) — degeneracy/model paper
arXiv:2412.12283 (ApJ submitted). — Approximate analytic Schwarzschild+Doppler light-curve model
for **two antipodal** equal spots, MCMC-fit to synthetic NICER data. Isotropic emission (h̄ = 0)
leaves one parameter unconstrained and biases inferred **compactness** by ≈ −0.02 in u; a
reparameterization removes the worst degeneracies. Works entirely in **absolute, energy-resolved
Fourier harmonics Aₙ(E)** — never defines a pulsed fraction. Fig. 7 is an (inclination, colatitude)
map carrying Poutanen–Beloborodov visibility classes I–IV, but argues the compactness bias is
*geometry-independent*. **Relation to us (closest prior work):** we differ on (i) the observable —
pulsed fraction, not inferred u; (ii) varying azimuthal spot separation (they lock antipodal);
(iii) the PF-vs-shape decomposition; (iv) the tiling/PF-saturation criterion. The only mention of
the tiling mechanism is *their own* Appendix A aside ("occultation of one spot coincides with the
maximum of the other"), undeveloped — and the Bauböck (2015) paper it cites does **not** actually
contain it (verified full-text). NB: Zhao's bibliography mislabels "2015a"/"2015b" as the same ApJ 811, 144.
*Verify yourself:* Fig. 7 + Ctrl-F "pulsed fraction" (absent).

### Zhao, Psaltis & Özel (2024) — counterintuitive radius–beaming Letter
arXiv:2412.12284 (ApJ submitted). — Companion Letter to 2412.12283. Headline: fitting a
too-strongly-beamed (deep-heated) model to a more isotropic emitter biases inferred **compactness
downward**, (δu/u) ≈ −(δh/h); worked example u_fit ≈ 0.313 vs u_syn = 0.35 (~10%). Beaming is the
linear law I ∝ [1 + h(E,T) cos α], h = beaming factor (h = 0 isotropic). Central observable is the
**fractional amplitude A₁/A₀** (Eq. 29), energy-resolved — the *cousin* of pulsed fraction, equal
only for a pure sinusoid; harmonic-fundamental vs true min/max, no saturation ceiling. Two antipodal
spots, fixed; Fig. 1 is a (θ, ζ) compactness-bias map. **Relation to us:** same villain (wrong
beaming), different crime scene (inferred radius, not the observed PF); A₁/A₀ ≠ our PF precisely in
the double-peaked two-spot tiling regime our result lives in.

### Riffert & Mészáros (1988)
ApJ, 325, 207. [ADS: 1988ApJ...325..207R](https://ui.adsabs.harvard.edu/abs/1988ApJ...325..207R)
— "Gravitational light bending near neutron stars." Origin of the qualitative fact that **more
isotropic beaming → smaller pulsation amplitude** (cited as such in Zhao 2412.12284's intro).
**Relation to us:** prior art for the *direction* of the effect; our novelty is the *quantified*
ΔPF under a realistic limb-darkening swap, not the direction.

### Annala & Poutanen (2010)
A&A, 520, A76. [DOI: 10.1051/0004-6361/200912773](https://doi.org/10.1051/0004-6361/200912773) ·
arXiv:1008.2270 — Applied Beloborodov's four visibility classes (cf.
[Poutanen & Beloborodov (2006)](#poutanen--beloborodov-2006) Fig. 5) to **124 real X-ray pulsars**,
tying class → single/double-peaked morphology and constraining compactness + magnetic obliquity
(≲ 40°). **Relation to us:** precedent that the (i, Θ) visibility-class phase diagram exists and has
been applied to data — we **extend** it (non-antipodal azimuth axis, PF-saturation labeling), not invent it.

### Sotani & Miyamoto (2018)
Phys. Rev. D, 98, 044017. [DOI: 10.1103/PhysRevD.98.044017](https://doi.org/10.1103/PhysRevD.98.044017) ·
arXiv:1807.09071 — Extended the (i, Θ) visibility classification into the ultra-compact regime
(M/R > 0.284, every surface point visible); max/min flux ratio grows with compactness. **Relation to
us:** second precedent for a geometry-space pulse-profile phase diagram — same "extend, not invent" framing.

### Bauböck, Psaltis & Özel (2015)
ApJ, 811, 144. [DOI: 10.1088/0004-637X/811/2/144](https://doi.org/10.1088/0004-637X/811/2/144) ·
arXiv:1505.00780 — Effects of **hot-spot size** on radius measurements from **single-spot** pulse
profiles; (i, θ_s) contour maps of where spot-size bias matters; single-spot self-eclipse sharpens
the profile. **No beaming, no pulsed fraction, no two-spot tiling** (verified full-text). **Relation
to us:** the natural *foil* — the prior single-spot eclipse analysis our two-spot tiling generalizes.
This is the paper Zhao 2412.12283 Appendix A cites for the tiling aside, but it does not contain the
two-spot fill-in mechanism.

### Salmi et al. (2023)
ApJ, 956, 138. [DOI: 10.3847/1538-4357/acf49d](https://doi.org/10.3847/1538-4357/acf49d) ·
arXiv:2308.09319 — Atmosphere-model systematics (composition, ionization, empirical beaming) on real
J0030/J0740 NICER fits; atmosphere choice matters for J0030's radius, less for J0740. Works at the
posterior M–R level, no PF isolation. **Relation to us:** establishes that beaming/atmosphere
knowledge matters at the ~5% level; we recast that systematic onto the directly observable PF.

### Psaltis, Özel & Chakrabarty (2014)
ApJ, 787, 136. [DOI: 10.1088/0004-637X/787/2/136](https://doi.org/10.1088/0004-637X/787/2/136) ·
arXiv:1311.1571 — Prospects for M–R from pulse-profile modeling; bolometric profiles alone cannot
break the M–R degeneracy; constant-rms contours ∝ sin i · sin Θ; harmonic amplitudes carry the
information. **Relation to us:** background for amplitude/harmonic observables and the geometric
sin i · sin Θ scaling. (Distinct from Psaltis & Özel (2014) ApJ 792, 87, the Hartle–Thorne paper
already listed below under deferred citations.)

---

## Deferred / future-work citations (not used yet)

- **Ho & Lai (2001)**, MNRAS, 327, 1081 — magnetic atmospheres. The draft's
  "magnetic anisotropy" framing; deferred (the engine is Thomson-only).
- **Morsink et al. (2007)**, ApJ, 663, 1244 — oblate-Schwarzschild; only if a
  reviewer raises oblateness (oblateness/Doppler are out of project scope).
- **Psaltis & Özel (2014)**, ApJ, 792, 87 — special-relativistic / oblateness
  corrections; same out-of-scope note.

---

## Reference data sets

### L26 code-comparison waveforms (SD1a benchmark)
- **What:** ASCII reference pulse profiles from the IM code for the SD1/SD2/OS1
  test problems, distributed as the supplementary archive `apjlab5968.tar.gz`
  attached to [Bogdanov et al. (2019)](#bogdanov-et-al-2019--paper-ii-l26) on
  IOPscience. We use `SD1a_test_IM.txt` (phase in cycles; photon flux at 1 keV).
- **How to get it:** download the "Supplementary data" archive from the article
  page ([10.3847/2041-8213/ab5968](https://doi.org/10.3847/2041-8213/ab5968)) and
  extract it to `data/l26_reference/`. Then `python scripts/code_comparison.py`
  and the code-comparison test (`tests/test_pulse.py`) pick it up automatically.
- **Why it is not committed:** the archive is third-party AAS material. AAS holds
  copyright on pre-2021-10-11 articles (L26 is 2019) unless gold-OA, and the
  archive's ReadMe grants *use* ("to facilitate testing of other independently
  developed codes") but not redistribution. So `data/` is gitignored for
  non-figure files and the data stays a local download; the test skips cleanly
  when it is absent. (AAS licensing policy:
  [journals.aas.org/article-charges-and-copyright](https://journals.aas.org/article-charges-and-copyright/).)
