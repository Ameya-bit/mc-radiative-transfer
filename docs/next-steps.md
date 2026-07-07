# Next Steps — The Paper-Ready Sprint (v0.9.8 → v1.0.0-rc)

> Living handoff/roadmap doc. Snapshot as of **2026-07-06**, after v0.9.7 (convergence redo).
> Purpose: a **one-day, parallelized sprint plan** that closes every open attention point from
> [`pipeline-physics.md`](pipeline-physics.md) — by simulation change, by measurement, or by
> a drafted paper explanation — so that the research paper can be **started the next morning**
> with nothing load-bearing left unverified.
> Bibliography + novelty notes: [`references.md`](references.md). Attention-point audit:
> [`pipeline-physics.md`](pipeline-physics.md) §"Things that need attention". Latest deep dive:
> [`deep-dives/v0.9.7-convergence-redo.md`](deep-dives/v0.9.7-convergence-redo.md).

---

## The science in one paragraph

Most pulse-profile models assume hotspots emit **isotropically**; realistic limb-darkened
**beaming** I(μ) sharpens the pulse and raises the **pulsed fraction** (PF). The systematic is
**ΔPF = PF_real − PF_iso ≈ +0.15** (converged value; the old +0.16 was a 1σ-high realization).
Whether it shows up *in the pulsed fraction* or *hides in the waveform shape* depends on
geometry: if two anti-phased spots **tile** the rotation (combined flux never hits zero), PF
stays unsaturated and the systematic lands in PF; if the spots leave a dark gap (PF pins at 1),
it hides in shape. **Tiling — not single-spot eclipse — is the discriminator.** And (v0.9.10)
**spin is a second router**: at the real 346.5 Hz the Doppler boost saturates the min/max PF
statistic, draining ΔPF to +0.03–0.07 (all geometry terms included) while the waveform-shape
difference (~10% RMS) survives intact — the systematic moves, it does not shrink.

## Current status

- **v0.9.0–v0.9.7 complete:** systematic established (v0.9.0), J0030 shape-only anchor (v0.9.1),
  J0740 live-PF anchor (v0.9.2), phase diagram (v0.9.4), Beloborodov/PB06 validation +
  Zhao positioning (v0.9.5), finite-cap bound (v0.9.6), convergence redo (v0.9.7:
  ΔPF = **+0.1486 ± 0.0026** at τ = 10, Riley; b(τ) is a plateau, not a peak).
- **Attention-point audit done** (`pipeline-physics.md`): 8 items, dispositions assigned below.
- **This sprint = the last engineering milestone before v1.0.0 (the paper).**

## Definition of "paper-ready" (the finish-line checklist)

The sprint is done when every box below is checked:

- [x] Production beaming library regenerated: escape-matched N, 5 seeds, per-bin σ, b ± σ stored
  — **done (v0.9.7.1)**, [deep-dive](deep-dives/v0.9.7.1-production-library.md)
- [x] All downstream results re-run off it; headline quoted as **ΔPF ± σ** (Riley + Miller)
  — **done (v0.9.7.2)**: Riley **+0.144 ± 0.003**, Miller **+0.201 ± 0.005** (τ = 10, J0740);
  [deep-dive](deep-dives/v0.9.7.2-downstream-rerun.md)
- [x] Exact Schwarzschild bending implemented + validated + **ΔPF shift measured** — **B1+B2 done
  (v0.9.9)**, [deep-dive](deep-dives/v0.9.9-exact-bending.md): shift −0.0066 (2.1σ) at Riley
  u = 0.494 → Gate G1 **switch to exact**; corrected headline **+0.137 ± 0.003 (Riley) /
  +0.195 ± 0.005 (Miller)**
- [x] Doppler + aberration implemented, validated, and the beaming coupling **measured** at 346.5 Hz
  — **C1–C4 done (v0.9.10)**, [deep-dive](deep-dives/v0.9.10-doppler.md): Gate G2 fired at 37σ;
  bolometric ΔPF +0.137 → +0.037 (Riley, δ⁴), in-band photon +0.019, all-geometry (warp
  included, C4) **+0.030 (Riley) / +0.056 (Miller)**; the systematic **re-routes into waveform
  shape** (RMS ≈ 0.10, invariant to spin/band/warp), it does not vanish — framing: spin is a
  second router alongside geometry; oblateness bounded ≤ 0.001
- [x] Engine transport validated **exactly** against H(μ) (isotropic-scattering mode) — **Track D
  done (v0.9.11)**, [deep-dive](deep-dives/v0.9.11-isotropic-transport.md): opt-in isotropic phase
  function; τ = 10, 5 seeds, 4.1×10⁵ escaped → flux-space reduced **χ² = 0.70** (dof 17) vs H(μ)·μ;
  Thomson control χ² = 1.49 (2.8%) = the ~1–3% phase-function sensitivity ⇒ H is **near-exact**,
  not exact, for the real atmosphere
- [x] First-bin / clamped-tail sensitivity of ΔPF quantified (perturbation + H-tail splice) — **Track E2
  done (v0.9.12)**, [deep-dive](deep-dives/v0.9.12-tail-and-robustness.md): resample σ + H-tail splice on the
  two μ<0.1 bins shift ΔPF by ≤ 0.006 (≤ the ±0.003 seed bar), static and at 346.5 Hz — tail not load-bearing
- [x] Atmosphere-law robustness row: ΔPF repeated with an analytic limb-darkening law — **Track E3 done
  (v0.9.12)**: Eddington (b=1.50) and Chandrasekhar H(μ) (b=1.68) reproduce the positive, spin-diluted ΔPF
  (+0.13…+0.20 static, +0.03…+0.06 at 346.5 Hz); magnitude tracks slope b — same sign + geometry-routing
- [ ] Scope paragraphs drafted: bolometric-vs-band (item 6), boundary semantics (item 7)
- [ ] `SHAPE_TAU` comment + all docs reworded "plateau", stale +0.16/+0.23 numbers propagated out
- [ ] Full test suite green, incl. new tests for bending/Doppler/isotropic modes
- [ ] Advisor literature-verification request (re-)sent — the only item that can stay pending

## Attention-point → workstream map

| # (pipeline-physics.md) | Item | Disposition | Closed by |
|---|---|---|---|
| 1 | Beloborodov bending at u = 0.494 | **measure** | Track B |
| 2 | Isotropic-H vs Thomson conflation | **measure + caveat wording** | Track D |
| 3 | b(τ) non-monotonic | resolved (v0.9.7) | Track A (implements) |
| 4 | Noisy clamped first bin | **measure** | Track E2 |
| 5 | Slow-rotation limit at 346 Hz | **measure** | Track C |
| 6 | Grey / bolometric scope | **explain in paper** | Track F |
| 7 | Bottom-boundary semantics | **explain in paper** (backed by Track D) | Track F |
| 8a/8b | Normalization + bin-center bias | **fix in code** | Track A |
| 8c–8e | mfp check, weight units, τ-grid docstring | explain / leave (comments exist) | Track F (one pass) |

---

## The one-day sprint

**Principle:** the only long job is the 80M-photon library run (~30–90 min on 16 workers).
Everything else is either deterministic single-core code (bending, Doppler, robustness rows —
runs *while* the library computes) or seconds-cheap interpolation that *needs* the new library
(anchors, sensitivity — runs *after*). So: **launch the library first, build the deterministic
physics in parallel, integrate in the afternoon.**

```
 morning                                  afternoon
 ────────────────────────────────────────┼──────────────────────────────────────
 A1 rework tau_sweep ─► A2 RUN (bg, ~1h) ─► A3 downstream re-runs ─► E1 numbers
 B1 bend_exact + tests ──────────────────► B2 anchor re-run (exact) ─► gate G1
 C1 doppler layer + tests ───────────────► C2 anchor re-run (346Hz) ─► gate G2
 D1 iso-scatter sampler ─► D2 RUN (short, after A2 frees cores) ─► D3 H(μ) overlay
                                          ► E2 tail sensitivity (needs A3)
                                          ► E3 robustness row (needs A3)
 F  paper-scope paragraphs + rewording (any idle moment; no dependencies)
```

Tracks B, C, D1, and F have **zero dependency** on Track A and on each other — they can be
developed in parallel sessions/worktrees. A3/E1–E3 are blocked only on the library landing.

### Track A — Production beaming library (the long pole; start immediately)

**A1. Rework `scripts/tau_sweep.py`** (~45 min of code) — ✅ **DONE (v0.9.7.1)**
- Escape-matched injected N per (τ, seed) targeting **4×10⁵ escaped photons each**, 5 seeds
  (v0.9.7 §3.4 budget, ≈ 80M injected total):

  | τ | 0.1 | 0.3 | 1.0 | 3.0 | 10 | 30 |
  |---|---|---|---|---|---|---|
  | injected/seed | 4.4×10⁵ | 5.0×10⁵ | 7.2×10⁵ | 1.3×10⁶ | 3.4×10⁶ | 9.5×10⁶ |

- Per-(τ, seed) `SeedSequence(BASE_SEED, spawn_key=...)` streams (kill the sequential RNG
  that couples rows); reuse the v2 study's process-parallel worker pool.
- **Fold in fixes 8a + 8b** while touching the file: normalize I(μ) by the fitted `a`
  (not the single μ ≈ 0.975 bin), and record each bin's **mean escaped μ** (divide counts by
  that instead of the bin center — unbiases the low-μ bins).
- Store in `beaming_library.npz`: pooled `intensity_by_tau`, **per-bin seed std**,
  `b_of_tau` ± σ, per-seed curves (E2 needs them), escape fractions, provenance metadata.

**A2. Launch the run in the background** (~30–90 min wall, 16 workers). Do not wait on it —
switch to Track B/C immediately. — ✅ **DONE (v0.9.7.1):** ran 71.9 min wall (30 tasks, 79.3M
photons); b(τ) monotone to a ~1.8 plateau, τ = 10 ≡ τ = 30, every row hit 4×10⁵ escaped.

**A3. Downstream re-runs** (after A2 lands; pure interpolation + geometry, ~minutes total): — ✅ **DONE (v0.9.7.2)**
> ΔPF ± σ: Riley **+0.144 ± 0.003**, Miller **+0.201 ± 0.005** (τ = 10, J0740; gate 1.2σ vs
> v0.9.7 +0.1486, PASS). J0030 shape-RMS ~6% (PF saturated). finite-cap / phase-diagram re-run
> clean; verdicts unchanged. Driver: `scripts/a3_seed_errors.py`. See
> [v0.9.7.2 deep-dive](deep-dives/v0.9.7.2-downstream-rerun.md).
- `j0740_anchor.py` + `j0030_anchor.py` once per seed-library → **ΔPF ± σ** (expect ≈ +0.15
  Riley / Miller shifts proportionally); `phase_diagram.py`; `finite_cap.py`.
- Sanity gates: b(τ) monotone within error; ΔPF(τ = 10) within ~2σ of the v0.9.7 converged
  +0.1486 ± 0.0026; τ = 10 I(μ) still lands between Eddington and H(μ).

### Track B — Exact Schwarzschild bending (attention item 1; parallel with A2)

**B1. Implement + validate** (~2–3 h; deterministic, no MC, no core rewrites) — ✅ **DONE (v0.9.9)**
- `src/mcrt/bending.py`: `deflection_angle` (numpy Gauss–Legendre quadrature of the ray integral),
  `ExactBending(u)` (tabulates ψ(α), inverts to cos α(cos ψ), exports the exact lensing Jacobian
  D(ψ) = d cos α/d cos ψ), and `bend_exact`. Wired into `pulse.point_spot_flux`/`compute_profile`
  via an optional `bending=` arg — linear default **bit-for-bit unchanged**; the flux uses D(ψ)
  in place of the constant (1 − u). Photon-sphere guard (u < 2/3).
- Validated (`tests/test_bending.py`, 24 tests, suite 82 green): u→0 gives ψ = α; independent
  fine quadrature to <1e-4; Jacobian vs an analytic derivative to ~1e-4 incl. grazing; D(cosψ=1)
  = (1−u) exactly; Beloborodov agreement (gap grows with u). **Free self-check: the SD1a
  code-comparison residual shrinks 0.80% → 0.11%** vs the independent IM code — the ~1% gap was
  the linear approximation. [deep-dive](deep-dives/v0.9.9-exact-bending.md).

**B2. Measure** (production library, per seed) — ✅ **DONE (v0.9.9).** Driver
`scripts/b2_exact_bending.py` (reproduces the linear branch bit-for-bit as its control) →
`data/b2_exact_bending.npz`. Result at τ = 10:

| anchor | u | ΔPF_linear ± σ | ΔPF_exact ± σ | shift | \|shift\|/σ |
|---|---|---|---|---|---|
| Riley  | 0.494 | +0.1437 ± 0.0032 | **+0.1371 ± 0.0027** | −0.0066 | **2.1σ** |
| Miller | 0.444 | +0.2012 ± 0.0052 | **+0.1952 ± 0.0053** | −0.0060 | 1.1σ |

Resolution-stable to 5 dp (n_alpha 2048…16384); vanishes in the flat limit. J0030 (eclipsing)
**unmoved** (ΔPF ≡ 0, floor 0 under both maps); phase-diagram star markers unmoved (coords are
bending-independent; ΔPF value shifts ≤ 0.007, stays on the same side of the boundary).

**Gate G1 verdict — SWITCH to exact, ENACTED (v0.9.9.1).** Riley's 2.1σ > σ → the honest headline
for the J0740 anchors is the bending-corrected **ΔPF = +0.137 ± 0.003 (Riley) / +0.195 ± 0.005
(Miller)**. Now live in the pipeline: `j0740_anchor.py` `run()` + `phase_diagram.py` star markers
build `ExactBending(u)` per anchor; the broad phase-diagram sweep stays linear (tolerance: gap
≤ 0.007 at the compact edge); `a3_seed_errors.py` stays linear as the convergence cross-check. Same
size/sign/geometry-routing — a tightening, not a reversal — converting the referee's best attack
into a table row. **Bonus:** exact bending exposes that linear's past-grazing extrapolation hid
Riley's ≈14% per-spot eclipse; the anti-phased pair still tiles (floor 0.70), so the tiling-not-
eclipse thesis sharpens. **+0.137 is the 0 Hz baseline**; Track C/G2 (Doppler) stacks on top.

### Track C — Doppler + aberration layer (attention item 5; parallel with A2 and B)

**C1. Implement + validate** — ✅ **DONE (v0.9.10).** `src/mcrt/rotating.py`: β(θ_s) =
2πνR sinθ_s/(c√(1−u)) — **β_eq ≈ 0.127 Riley / 0.134 Miller** (this plan's earlier
0.088/0.10 had dropped the √(1−u) redshift enhancement) — Doppler δ(φ), aberration
μ' = δ cos α before the I(μ) lookup, δⁿ boost (n = 4 energy / 3 photon), opt-in
`rotation=` (ν → 0 bit-for-bit). Validated against the 200 Hz SD1c code-comparison
waveform (max|Δ| 1.36%, the linear-bending floor) + a from-scratch 4-vector boost
cross-check to machine precision.

**C2. Measure + shape routing** — ✅ **DONE (v0.9.10).** `scripts/c2_doppler_coupling.py`
→ `data/c2_doppler_coupling.npz` + four-pulse figure
(`data/pulse_profile_doppler_routing.png`). The **δ⁴ boost — not aberration — dominates**
(the naïve boost-cancels expectation was wrong; PF is a nonlinear max/min statistic):
ΔPF(0 → 346.5 Hz) = **+0.137 → +0.037 (Riley, 37σ)**, **+0.195 → +0.061 (Miller, 25σ)**;
J0030 unmoved (rotation cannot un-eclipse a spot). **But the systematic re-routes, it
does not vanish**: the iso-vs-real waveform-shape difference (peak-normalized RMS
0.107 → 0.105 Riley) and the fundamental-harmonic gap Δ(A1/A0) (+0.033 → +0.029) are
nearly spin-invariant — spin pumps a large *common* second harmonic into both models
(A2/A0: 0.02 → 0.19 iso) and saturates the min/max PF statistic as a probe of beaming.

**C3. Band-limited spectral variant** — ✅ **DONE (v0.9.10).** `scripts/c3_band_doppler.py`
+ `mcrt.rotating.BandSpectrum`/`band_boost` (wide-band limits recover δ³/δ⁴ exactly;
7 new tests, suite **108 green**). Blackbody stand-in at the published kT (0.0842 keV
Riley / 0.094 keV Miller), calibrated NICER bands (0.3–1.5 / 0.3–1.24 keV). The band
sits on the **Wien tail** (in-band effective exponent n_eff ≈ 6.1 Riley / 5.5 Miller),
so the dilution **deepens**: **ΔPF(band photon) = +0.019 ± 0.002 (Riley), +0.045 ± 0.003
(Miller)**; band-edge/kT scans move Riley across +0.009…+0.024 — the residual PF number
is both small *and* band/kT-sensitive, so it cannot carry a headline. In-band shape
routing intact (RMS 0.106, Δ(A1/A0) +0.024).

**Gate G2 verdict — FIRED (37σ), framing = ROUTING REFRAME (decided 2026-07-07).** The
paper reports the beaming systematic three-part: (i) **+0.137/+0.195** is the
slow-rotation bolometric *mechanism*; (ii) **geometry routes it** (tiling → PF,
eclipse/dark-gap → shape); (iii) **spin is a second router** — at the real 346.5 Hz the
min/max PF retains only +0.02…+0.06 of it (convention/band-dependent) while the
waveform-level cost of assuming isotropy (~10% RMS, most of Δ(A1/A0)) is spin-invariant.
PF is a spin-degraded probe of beaming; the shape difference is the conserved quantity.

**C4. Caveat re-audit vs the diluted number** — ✅ **DONE (v0.9.10).**
`scripts/c4_caveat_audit.py` (paired per seed; `mcrt.rotating.travel_time_delay`, eq. 18
quadrature generalized out of the C1 driver, flat-space closed form pinned in tests).
- **Light-travel delay: ESCALATED — now included, not caveated.** Exactly PF-invariant for
  a single spot (measured <2e-6) but 2–5σ against the diluted numbers for the two-spot
  sums, and *positive* (de-synchronizing the spots' extremes partially restores the gap).
  **All-geometry rotating residuals: Riley +0.041 (δ⁴) / +0.030 (band photon); Miller
  +0.069 / +0.056.** Shape-routing metrics unchanged by the warp.
- **Oblateness: caveat STANDS, with a measured bound.** AM14 first-order perturbation
  (o₂ ≈ −1.7%/−2.5%; near-equatorial spots ⇒ ΔR/R ≤ 0.2%, normal tilts ≤ 0.6°): ΔPF shift
  +0.0011 (0.5σ) Riley, −0.00002 Miller. Residual oblate-bending remainder is OS-level,
  same order — one caveat sentence with these numbers.

### Track D — Exact transport validation (attention item 2; small) — ✅ **DONE (v0.9.11)**

`src/mcrt/utils.py` `sample_isotropic_angle` + `monte_carlo.Simulation(phase_function=...)`;
driver `scripts/d_isotropic_validate.py` → `data/d_isotropic_h.npz` + overlay figure;
`tests/test_transport.py` (8 tests, suite **119 green**).
[deep-dive](deep-dives/v0.9.11-isotropic-transport.md).

**D1. ✅** Opt-in **isotropic-scattering** phase function (P(μ) = 1/2, μ uniform on [−1, 1]) via a
`SCATTER_SAMPLERS` dispatch; default `"thomson"` is **bit-for-bit unchanged** (same sampler call,
same RNG stream — asserted in tests), invalid names raise.
**D2. ✅** One thick slab τ = 10, 5 seeds × 7×10⁵ injected (**4.1×10⁵ escaped**), ~80 s on 10 workers.
**D3. ✅** The match is **exact within error** — but in **flux space**, not intensity space. The
naïve I(μ) = counts/μ overlay fails (χ² = 11.7) on the attention-fix-8b **bin-center bias**
(dividing ∫Iμ dμ by the bin center, a systematic that grows with N and masks the signal). Testing
the raw escaping distribution against N·∫H(μ)μ dμ directly gives **flux-space reduced χ² = 0.70
(dof 17, max residual 1.6σ)** — isotropic transport reproduces H(μ). The **Thomson control** (same
machinery, physical dipole) sits at χ² = 1.49 (2.8% max deviation) = the ~1–3% phase-function
sensitivity (Chandrasekhar Ch. X). Paper wording (Track F, item 7): the Thomson-vs-isotropic-H
overlay is a **near-exact** reference; never write "matches the exact solution" without that caveat.

### Track E — Post-library quantifications (afternoon; each ≤ 1 h, all interpolation-cheap)

**E1. Number propagation:** new ΔPF ± σ through README, deep dives, figures; reword
`SHAPE_TAU = 10` in `anchor_lib.py` to "on the saturated b(τ) plateau" (not "where b(τ) peaks").
**Now runs LAST (v0.9.10 ripple):** it propagates the three-part routing framing (mechanism +
geometry router + spin router), so it is blocked on C4 and the E2/E3 rotation-on rows.
**E2. Tail sensitivity (item 4)** — ✅ **DONE (v0.9.12).** `scripts/e2_tail_sensitivity.py` →
`data/e2_tail_sensitivity.npz` + figure; [deep-dive](deep-dives/v0.9.12-tail-and-robustness.md).
(a) Resample the two μ<0.1 bins within their measured per-bin σ (500 draws) → induced σ(ΔPF);
(b) splice a Chandrasekhar-H-shaped tail below μ=0.1 in place of the clamp → |δ(ΔPF)| tail systematic.
**Both effects ≤ 0.006** — at or below the ±0.003 seed bar — static and at 346.5 Hz; the clamped
grazing tail is not load-bearing. Run rotation-on as the ripple demanded: the one place the tail
bites is **Miller at spin** (|δ|=0.006), where the grazing spots let aberration (μ'=δ cos α) pull
the tail into the min/max PF statistic. Both drivers reproduce the +0.137/+0.195 → +0.037/+0.061
headline exactly as their baseline.
**E3. Atmosphere-law robustness row** — ✅ **DONE (v0.9.12).** `scripts/e3_atmosphere_laws.py` →
`data/e3_atmosphere_laws.npz` + figure. Repeated the J0740 swap under **Eddington 1+1.5μ** (b=1.50)
and **Chandrasekhar H(μ)** (b=1.68) alongside the Thomson slab (b=1.79). All three give a positive,
PF-live static ΔPF (+0.128…+0.196) diluting to +0.026…+0.061 at 346.5 Hz, with waveform-shape RMS
~0.10–0.13 throughout — **same sign and geometry-routing under independent laws**, magnitude tracking
slope b. Defuses "your slab isn't a hydrogen atmosphere": the effect is a property of limb darkening.

### Track F — Paper-scope drafting (no dependencies; fill idle moments)

Draft the four guard paragraphs **now**, while the context is loaded — they gate the paper more
than any number does:
1. **Bolometric-vs-band scope (item 6):** the paper quantifies the *mechanism* (limb-darkened
   beaming inflates PF where geometry leaves headroom) at bolometric level; hydrogen-atmosphere
   beaming in the NICER band is likewise limb-darkened, so sign and geometry-routing carry
   over; band-exact magnitude = future work. Position explicitly against NSATMOS/X-PSI.
2. **Boundary semantics (item 7):** I(μ; τ) is the transmission beaming of a passive scattering
   layer over a Lambertian photosphere with an absorbing base, conditioned on escape —
   converging to the semi-infinite result at large τ (now demonstrated by Track D).
3. **Framing guardrail (the big one):** the NICER teams do **not** assume isotropy — X-PSI and
   Illinois–Maryland use beamed model atmospheres. Frame every claim as "what the isotropy
   assumption costs / where geometry routes the beaming systematic," **never** "published
   results are biased by +0.15." Identical math; opposite referee outcome.
4. **8c–8e sweep:** one pass confirming the existing comments (mfp check tolerance, per-anchor
   weight units, τ-grid no-interpolation docstring) say what the paper needs; fix the
   `beaming_library` docstring that implies τ-interpolation.
5. **PF-statistic scope (new, v0.9.10):** the spin dilution is specific to the min/max
   PF = (F_max−F_min)/(F_max+F_min) — the Doppler boost pumps a common second harmonic that
   saturates the extremes. An A1-based pulsed amplitude dilutes only ~12% (Δ(A1/A0)
   +0.033 → +0.029), and the peak-normalized waveform difference not at all. The paper must
   say which statistic collapses and why, or readers translating "ΔPF" into their preferred
   pulsed-fraction definition will get a different answer. NB the framing paragraphs (items
   1 and 3) now carry the three-part claim: mechanism +0.137/+0.195 → routed by geometry →
   diluted in PF (not in shape) by spin.

### External (not completable today)

- **Advisor literature verification** (Dr. ud-Doula): human check of the two Zhao 2024 papers,
  Bauböck 2015, and the PB06/Annala/Sotani lineage before novelty claims are finalized.
  Action today: (re-)send the request with the updated numbers. Drafting can start before the
  reply; the novelty *sentences* stay provisional until it lands.

---

## Suggested version mapping

| Version | Contents | Track |
|---|---|---|
| v0.9.8 | Production library + downstream re-runs + 8a/8b fixes + SHAPE_TAU rewording | A, E1 |
| v0.9.9 | Exact bending: implementation, validation, measured ΔPF shift (gate G1) | B |
| v0.9.10 | Doppler/aberration: implementation, validation, coupling + **shape routing** + **band-limited variant** + **caveat audit** (gate G2 + routing framing) | C1–C4 |
| v0.9.11 | Isotropic-scattering H(μ) validation (**Track D done**) | D ✅ |
| v0.9.12 | Tail sensitivity (E2) + atmosphere-law robustness (E3), both **incl. rotation-on variants** | E2 ✅, E3 ✅ |
| v1.0.0-rc | Scope paragraphs merged, all numbers propagated, suite green → **start the paper** | F |

Each lands with the usual README progress-log entry + deep-dive; commit refs backfilled per
convention. Original priority **A > B > C > E2 > D > E3**; with A, B, C (incl. C4), **D**, **E2**,
and **E3** now landed the only remaining sprint item is **E1** (number propagation, in progress) —
it propagates the routing framing off the numbers E2/E3 produced — then Track F scope paragraphs.

## Novelty position (unchanged; pending human verification)

Closest prior work: **Zhao, Psaltis & Özel 2024** (arXiv:2412.12283/12284) — wrong-beaming
systematic routed into *radius* via A₁/A₀; never pulsed fraction, never azimuthal separation,
never tiling. **Bauböck et al. 2015** is the foil (single-spot spot-size; Zhao's tiling
attribution to it is not in its text). Phase-diagram lineage (PB06; Annala & Poutanen 2010;
Sotani & Miyamoto 2018) is antipodal-only, never PF-saturation-labeled. Our three contributions:
1. A single **ΔPF ± σ number** for the isotropic→realistic beaming swap.
2. A **geometry phase diagram with azimuthal separation as an axis**, labeled PF-vs-shape.
3. The **tiling / PF-saturation criterion** (combined flux touches zero ⇔ PF pins at 1).

## Key files

| Path | What |
|---|---|
| `scripts/tau_sweep.py` | **(rework — A1)** escape-matched, 5-seed production library |
| `src/mcrt/bending.py` | **(new — B1)** exact ray integral + lensing Jacobian |
| `scripts/b2_exact_bending.py` | **(new — B2)** per-seed linear-vs-exact ΔPF shift + Gate G1 verdict |
| `src/mcrt/pulse.py` | gains optional `bending=` / rotational arguments (defaults untouched) |
| `src/mcrt/monte_carlo.py` + `utils.py` | gain opt-in isotropic-scattering mode (D1) |
| `scripts/anchor_lib.py` | `SHAPE_TAU` rewording (E1); robustness-row helpers (E3) |
| `scripts/j0740_anchor.py` / `j0030_anchor.py` | re-run per seed → ΔPF ± σ (A3) |
| `docs/pipeline-physics.md` | the attention-point audit this sprint closes |
| `docs/deep-dives/v0.9.7-convergence-redo.md` | budget table + method for Track A |

## Working notes

- **Fact-Forcing Gate hook:** intercepts file edits and some Bash calls; quote the current
  instruction verbatim + state file facts, then retry. Comply and retry.
- **Commits:** straight to `main`; version entries in `README.md` with `<pending>` ref,
  backfilled via `docs: backfill vX.Y.Z commit ref`.
- **Compute:** 16-worker process pool (v0.9.7 setup). The library run owns it for ~an hour;
  Tracks B/C/F are single-core and unaffected; schedule D2 after A2 completes.
- **Numbers to retire on sight:** +0.16 (Riley) / +0.23 (Miller) — superseded by the
  production-library ΔPF ± σ (expected ≈ +0.15 Riley) the moment A3 lands.
