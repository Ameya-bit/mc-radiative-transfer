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
