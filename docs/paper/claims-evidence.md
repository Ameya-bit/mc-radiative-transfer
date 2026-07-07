# Claims & Evidence — what the paper says and what backs it

> The **last link of the provenance chain** (code → equation → reference →
> **claim**). Every statement the paper will make, each pinned to: the number and
> its error bar, the script and figure that produce it, the
> [`pipeline-physics.md`](pipeline-physics.md) equations behind it, the
> [`references.md`](references.md) source, and a status flag. Read this to defend
> the paper one sentence at a time — and use it as the outline to draft prose from.
>
> *Numbers are the v0.9.12 production values; where a number is superseded by a
> later, more complete one, the newer replaces it (never quote the retired
> figures +0.16 / +0.23 or the pre-bending +0.144 / +0.201).*

## How to read a claim

Each claim below carries the same fields:

- **Number** — the result and its ± seed error (or "exact to grid" if deterministic).
- **From** — the driver script and the figure it produces (`data/…`).
- **Physics** — the equations in [`pipeline-physics.md`](pipeline-physics.md) that
  the number comes out of (e.g. "eqs (1)–(7)").
- **Refs** — the [`references.md`](references.md) sources.
- **Status** — one of:
  - **SOLID** — measured and validated; no caveat changes the claim.
  - **CAVEATED** — true only with a condition attached; the caveat must travel
    with the number wherever it is quoted.
  - **NOVELTY — PENDING** — believed to be new; awaiting Dr. ud-Doula's
    literature check before the novelty *sentence* is finalized.
  - **GUARDRAIL** — not a result but a constraint on *how* a result may be
    stated. These are the referee-facing traps; getting one wrong flips the
    paper's reception.

---

## Part A — The primary results

These four are the spine of the paper. A1 is the mechanism; A2–A4 are the three
things that happen to it.

### A1 — Realistic beaming raises the pulsed fraction

**Claim.** Replacing the usual isotropic-emission assumption ($I \equiv 1$) with a
scattering-atmosphere beaming function $I(\mu;\tau)$ makes a hot spot's pulse
*sharper*, raising the pulsed fraction. For the J0740 geometry (where the spots
tile — see A2), the effect lands in PF directly.

- **Number:** $\Delta\mathrm{PF} = +0.137 \pm 0.003$ (Riley fit) /
  $+0.195 \pm 0.005$ (Miller fit) — static, bolometric, exact bending, $\tau = 10$.
- **From:** `scripts/j0740_anchor.py` · figure `data/pulse_profile_j0740.png`.
- **Physics:** eqs (1)–(7) (the deterministic pulse + exact bending) fed by
  eqs (15)–(18) (the beaming library).
- **Refs:** geometry [Poutanen & Beloborodov (2006)](references.md#poutanen--beloborodov-2006);
  fits [Riley et al. (2021)](references.md#riley-et-al-2021--j07406620-l27) /
  [Miller et al. (2021)](references.md#miller-et-al-2021--j07406620-l28).
- **Status:** **SOLID** for the static/bolometric mechanism. Note the number
  *dilutes* once spin and the instrument band are added — that is claim A4, not a
  weakening of A1.

### A2 — Geometry routes the systematic (tiling vs. shape)

**Claim.** Whether the beaming systematic shows up *in the pulsed fraction* or
*hides in the waveform shape* is decided by geometry. If two anti-phased spots
**tile** the rotation (their combined flux never drops to zero), PF is
unsaturated and the systematic lands in PF (the J0740 case, A1). If the spots
leave a dark gap, PF pins at 1 (saturated), $\Delta\mathrm{PF} \approx 0$, and the
same systematic moves into the *shape* of the curve instead.

- **Number:** J0030 (does not tile) → $\Delta\mathrm{PF} \approx 0$, systematic
  appears as waveform-shape RMS ~6%. J0740 (tiles) → the live A1 number.
- **From:** `scripts/j0030_anchor.py` (figure `data/pulse_profile_j0030.png`);
  the map over all geometries is `scripts/phase_diagram.py`
  (figure `data/phase_diagram.png`).
- **Physics:** eqs (1)–(6); the saturation condition is eq (2) hitting zero on the
  spot sum eq (3).
- **Refs:** visibility-class lineage [Poutanen & Beloborodov (2006)](references.md#poutanen--beloborodov-2006).
- **Status:** **SOLID.** The routing is a measured, reproduced behavior. The
  *criterion* that governs it is the novelty claim A3.

### A3 — The tiling criterion (the discriminator)

**Claim.** PF saturates at 1 exactly when the *summed* multi-spot flux touches
zero — i.e. when the spots fail to tile the rotation — and **not** when any single
spot is eclipsed. This is the sharp line between the two regimes in A2.

- **Evidence:** J0740-Miller's two spots each hide ~21% of the cycle (each is
  individually eclipsed for part of the rotation), yet the anti-phased pair still
  tiles (combined flux floor 0.70) and PF stays live. Exact bending (A-side,
  eq (7)) revealed that the linear map had hidden Riley's ~14% per-spot eclipse;
  the pair still tiles, so the criterion *sharpened* rather than broke.
- **From:** the analytic boundary in `scripts/phase_diagram.py` /
  `scripts/validate_phase_diagram.py`, built from `single_spot_eclipsed_fraction`
  in `scripts/anchor_lib.py`.
- **Physics:** eq (3) touching zero; boundary $\cos(i+\theta) < -u/(1-u)$ from eq (6).
- **Refs:** phase-diagram lineage is antipodal-only and never PF-saturation-labeled —
  [Poutanen & Beloborodov (2006)](references.md#poutanen--beloborodov-2006),
  [Annala & Poutanen (2010)](references.md#annala--poutanen-2010),
  [Sotani & Miyamoto (2018)](references.md#sotani--miyamoto-2018).
- **Status:** **NOVELTY — PENDING** (contribution 3 of 3; see [Part B](#part-b--the-novelty-position)).

### A4 — Spin is a second router (PF drains, shape survives)

**Claim.** At the real spin (J0740 = 346.5 Hz) the Doppler boost pumps a large,
*common* second harmonic (a twice-per-rotation Fourier component of the pulse —
see [`pipeline-physics.md`](pipeline-physics.md)) into both the isotropic and
realistic curves, saturating the min/max PF statistic. So the PF systematic drains to a small,
convention-dependent residual — **but the waveform-shape difference is nearly
spin-invariant.** The systematic moves; it does not vanish.

- **Numbers (the dilution ladder, Riley / Miller):**

  | configuration | Riley | Miller | source |
  |---|---|---|---|
  | 0 Hz static, bolometric (= A1) | $+0.137$ | $+0.195$ | `b2_exact_bending.py` |
  | 346.5 Hz, bolometric $\delta^4$ | $+0.037$ | $+0.061$ | `c2_doppler_coupling.py` |
  | 346.5 Hz, NICER-band photons | $+0.019$ | $+0.045$ | `c3_band_doppler.py` |
  | 346.5 Hz, **all geometry** (incl. light-travel warp), $\delta^4$ | $+0.041$ | $+0.069$ | `c4_caveat_audit.py` |
  | 346.5 Hz, **all geometry**, band photons | $+0.030$ | $+0.056$ | `c4_caveat_audit.py` |

  The honest spinning-PF range is **$+0.02\ldots+0.06$** (convention/band-dependent).
  Meanwhile the peak-normalized waveform-shape difference holds at **RMS ≈ 0.10**
  across spin, band, and warp, and most of the fundamental-harmonic gap
  $\Delta(A_1/A_0)$ — the once-per-rotation amplitude relative to the mean —
  survives ($+0.033 \to +0.029$).
- **From:** `scripts/c2_doppler_coupling.py`, `c3_band_doppler.py`,
  `c4_caveat_audit.py` · figures `data/pulse_profile_doppler_routing.png`,
  `data/doppler_dpf_vs_spin.png`, `data/doppler_dilution_ladder.png`,
  `data/doppler_band_amplification.png`.
- **Physics:** eqs (8)–(14) (the whole rotation layer); the band steepening is eq (14).
- **Refs:** [Bogdanov et al. (2019) §2](references.md#bogdanov-et-al-2019--paper-ii-l26).
- **Status:** **CAVEATED** — the caveat *is* the claim: quote the PF residual only
  alongside the spin-invariant shape difference, or the result reads as "the effect
  disappears" (it doesn't). Gate G2 (the pre-registered rule to report the anchors
  with spin included if it moved the result) fired at 37σ; the agreed framing is the
  three-part statement in [Part B](#part-b--the-novelty-position).

---

## Part B — The novelty position

Three contributions, **all pending Dr. ud-Doula's literature verification** before
the novelty sentences are finalized. Closest prior work is
[**Zhao, Psaltis & Özel (2024)**](references.md#zhao-psaltis--özel-2024--counterintuitive-radiusbeaming-letter)
(the wrong-beaming systematic routed into *radius* via $A_1/A_0$ — never pulsed
fraction, never azimuthal separation, never tiling) and
[**Bauböck et al. (2015)**](references.md#bauböck-psaltis--özel-2015)
(single-spot spot-size), both detailed in [`references.md`](references.md).

1. **A single $\Delta\mathrm{PF} \pm \sigma$ number** for the isotropic→realistic
   beaming swap. — **NOVELTY — PENDING**
2. **A geometry phase diagram with azimuthal separation as an axis**, labeled
   PF-vs-shape (not the antipodal-only diagrams of the PB06 lineage). — **NOVELTY — PENDING**
   (`scripts/phase_diagram.py`, `data/phase_diagram.png`)
3. **The tiling / PF-saturation criterion** = claim A3. — **NOVELTY — PENDING**

**Action:** the request to Dr. ud-Doula (human check of the two Zhao 2024 papers,
Bauböck 2015, and the PB06 / Annala–Poutanen / Sotani–Miyamoto lineage) can be
sent with the current numbers; drafting proceeds, but the novelty *wording* stays
provisional until he replies.

---

## Part C — Why the machinery is trustworthy

Validation claims — each says "this piece of the pipeline reproduces an
independent truth." These are what let a referee trust A1–A4.

### C1 — The transport engine is exact

**Claim.** In its isotropic-scattering mode the Monte Carlo reproduces
Chandrasekhar's exact $H(\mu)$ solution.

- **Number:** flux-space reduced $\chi^2 = 0.70$ (dof 17, max residual 1.6σ) vs
  $H(\mu)\cdot\mu$, at $\tau = 10$, 5 seeds, $4.1\times10^5$ escaped. The Thomson
  control sits at $\chi^2 = 1.49$ (2.8% max deviation) — the physical
  phase-function sensitivity, so $H$ is **near-exact, not exact**, for the real atmosphere.
- **From:** `scripts/d_isotropic_validate.py` · figure `data/d_isotropic_h_overlay.png`;
  `tests/test_transport.py`.
- **Physics:** eqs (16)–(22), isotropic mode (eq 21 with $P = \tfrac12$).
- **Refs:** [Chandrasekhar (1960)](references.md#chandrasekhar-1960).
- **Status:** **SOLID**, with the "near-exact" wording mandatory (see guardrail E5).

### C2 — The emergent brightness is physical

**Claim.** The $\tau = 10$ beaming curve lands between the two classical thick-atmosphere predictions.

- **Number:** slope $b \approx 1.79$, bracketed by
  [Eddington's](references.md#eddington-limb-darkening) $b = 1.5$ below and
  [Chandrasekhar's $H(\mu)$](references.md#chandrasekhar-1960) ($b \approx 1.68$)
  above; ≤ 2.8% from $H$.
- **From:** `scripts/tau_sweep.py` · figures `data/beaming_tau_curves.png`,
  `data/beaming_slope_vs_tau.png`.
- **Physics:** eqs (15)–(18).
- **Status:** **SOLID.**

### C3 — The pulse + bending machinery matches an independent code

**Claim.** At the Bogdanov SD1a benchmark our pulse profile matches the
independent Illinois–Maryland reference waveform.

- **Number:** residual 0.80% with linear bending → **0.11% with exact bending**.
- **From:** `scripts/code_comparison.py` · figure `data/pulse_profile_code_comparison.png`;
  `tests/test_bending.py`.
- **Physics:** eqs (2)–(7).
- **Refs:** [Bogdanov et al. (2019)](references.md#bogdanov-et-al-2019--paper-ii-l26)
  (SD1 suite); bending [Pechenick, Ftaclas & Cohen (1983)](references.md#pechenick-ftaclas--cohen-1983).
- **Status:** **SOLID.** The shrinking residual is independent evidence the exact
  map is right (the ~1% gap *was* the linear approximation).

### C4 — The Doppler + aberration layer matches an independent code

**Claim.** At the 200 Hz SD1c benchmark the rotation layer reproduces the
reference waveform.

- **Number:** max $|\Delta| = 1.36\%$ (the linear-bending floor); plus a
  from-scratch 4-vector boost cross-check to machine precision.
- **From:** `scripts/c1_doppler_validate.py` · figure `data/pulse_profile_doppler_sd1c.png`.
- **Physics:** eqs (8)–(13).
- **Refs:** [Bogdanov et al. (2019)](references.md#bogdanov-et-al-2019--paper-ii-l26).
- **Status:** **SOLID.**

### C5 — Exact bending was necessary, and its effect is measured

**Claim.** The linear bending approximation biased the headline enough to matter,
so the anchors (the two real pulsars J0030 and J0740, used as test cases) use exact
bending (pre-registered Gate G1).

- **Number:** shift $-0.0066$ (2.1σ) at Riley's $u = 0.494$; $-0.0060$ (1.1σ) at
  Miller's $u = 0.444$. Resolution-stable to 5 dp; vanishes in the flat limit.
- **From:** `scripts/b2_exact_bending.py` (reproduces the linear branch bit-for-bit
  as its own control).
- **Physics:** eq (6) vs. eq (7).
- **Status:** **SOLID.** Gate G1 **enacted** (v0.9.9.1) — a tightening, not a
  reversal; same size/sign/routing, and it converts the referee's best attack into
  a table row.

---

## Part D — Robustness (the result doesn't hinge on a modeling choice)

Each of these answers a specific "but what if…" a referee will raise.

### D1 — The grazing tail is not load-bearing

**Claim.** The noisy, clamped $\mu < 0.1$ tail of the beaming curve does not carry the result.

- **Number:** perturbing the two grazing bins within their seed σ, and splicing an
  $H$-shaped tail under the clamp, both shift $\Delta\mathrm{PF}$ by
  **≤ 0.006** (at or below the ±0.003 seed bar), static and at 346.5 Hz. The one
  place it bites is Miller-at-spin ($|\delta| = 0.006$), where grazing spots let
  aberration pull the tail into the PF statistic.
- **From:** `scripts/e2_tail_sensitivity.py` · figure `data/e2_tail_sensitivity.png`;
  `tests/test_sensitivity.py`.
- **Physics:** eqs (15)–(16).
- **Status:** **SOLID** (with the Miller-at-spin note recorded, not hidden).

### D2 — It is not an artifact of our specific slab

**Claim.** The sign and geometry-routing survive under independent limb-darkening laws.

- **Number:** Eddington ($b = 1.50$) and Chandrasekhar $H(\mu)$ ($b = 1.68$)
  alongside the Thomson slab ($b = 1.79$) all give a positive, PF-live static
  $\Delta\mathrm{PF}$ ($+0.128\ldots+0.196$) diluting to $+0.026\ldots+0.061$ at
  346.5 Hz; shape RMS ~0.10–0.13 throughout. **Magnitude tracks the slope $b$.**
- **From:** `scripts/e3_atmosphere_laws.py` · figure `data/e3_atmosphere_laws.png`.
- **Physics:** eq (15) swapped for the analytic laws; eq (18) for $b$.
- **Refs:** [Eddington](references.md#eddington-limb-darkening),
  [Chandrasekhar (1960)](references.md#chandrasekhar-1960).
- **Status:** **SOLID.** Defuses "your slab isn't a hydrogen atmosphere" — the
  effect is a property of limb darkening itself.

### D3 — The point-spot idealization is a small, measured bias

**Claim.** Treating each spot as a point (rather than a finite cap) is a known,
bounded approximation.

- **Number:** area-weighted ring tiling of the published finite cap moves
  $\Delta\mathrm{PF}$ by ~$-0.003$.
- **From:** `scripts/finite_cap.py` · figure `data/finite_cap.png`.
- **Physics:** eqs (3)–(4).
- **Status:** **CAVEATED** — quote the ~$-0.003$ bias when the point-spot reduction is stated.

### D4 — The photon budget is sufficient

**Claim.** The Monte Carlo is converged; the error scales as $N^{-1/2}$ with no systematic floor.

- **Number:** error $\propto N^{-1/2}$; production run 5 seeds × $4\times10^5$
  escaped per $\tau$ row.
- **From:** `scripts/convergence_study.py` / `convergence_study_v2.py` ·
  figure `data/convergence_v2_dpf.png`.
- **Physics:** eqs (19)–(22).
- **Status:** **SOLID.**

### D5 — Deliberately-bounded geometry terms (delay + oblateness)

**Claim.** The two omitted geometry effects are bounded and beaming-independent.

- **Light-travel delay: ESCALATED — now included, not caveated.** PF-invariant for
  a single spot (< 2e-6) but 2–5σ and *positive* for the two-spot sums (it
  de-synchronizes the spots' extremes, partially restoring the gap) — hence it is
  in the A4 all-geometry rows.
- **Oblateness: caveat STANDS, with a bound.** AlGendy–Morsink first-order
  perturbation → $\Delta\mathrm{PF}$ shift $+0.0011$ (0.5σ) Riley, $-0.00002$
  Miller — one caveat sentence with these numbers.
- **From:** `scripts/c4_caveat_audit.py`; `mcrt.rotating.travel_time_delay` (eq 18).
- **Refs:** [AlGendy & Morsink (2014)](references.md#algendy--morsink-2014).
- **Status:** **CAVEATED** (oblateness) / folded-in (delay).

---

## Part E — Framing guardrails (how claims must be stated)

Not results — constraints. Each is a way the paper can be *technically correct but
fatally mis-framed*. These matter more than any single number.

### E1 — Never say "published results are biased"

The NICER teams (X-PSI, Illinois–Maryland) do **not** assume isotropy — they use
beamed model atmospheres. Frame **every** claim as *"what the isotropy assumption
costs"* / *"where geometry routes the beaming systematic,"* never *"published
results are off by +0.15."* Identical math, opposite referee outcome.
**Status: GUARDRAIL (the big one).**

### E2 — "ΔPF" means the min/max statistic specifically

The spin dilution (A4) is specific to
$\mathrm{PF} = (F_{\max}-F_{\min})/(F_{\max}+F_{\min})$ — the Doppler boost pumps a
common second harmonic that saturates the *extremes*. An $A_1$-based pulsed
amplitude (the fundamental Fourier component of the pulse, defined in
[`pipeline-physics.md`](pipeline-physics.md)) dilutes only ~12%, and the
peak-normalized waveform difference not at all. The paper must say *which* PF
statistic collapses and why, or a reader using a different PF definition will get
a different number.
**Status: GUARDRAIL.**

### E3 — The bolometric-vs-band scope

The mechanism (A1) is quantified at **bolometric** (grey, all-energy) level.
Hydrogen-atmosphere beaming in the NICER band is likewise limb-darkened, so sign
and geometry-routing carry over (demonstrated by D2), but the band-exact
*magnitude* is future work. Position explicitly against the hydrogen-atmosphere
models (e.g. NSATMOS) and NICER inference codes (X-PSI, Illinois–Maryland); do not
claim a band-exact number.
**Status: GUARDRAIL.**

### E4 — Boundary semantics of $I(\mu;\tau)$

$I(\mu;\tau)$ is the *transmission* beaming of a passive scattering layer over a
Lambertian photosphere with an absorbing base, conditioned on escape — converging
to the semi-infinite result at large $\tau$ (demonstrated by C1). State this; do
not imply it is a full self-consistent atmosphere.
**Status: GUARDRAIL.**

### E5 — $H(\mu)$ is near-exact, not exact

For the real (Thomson) atmosphere the $H(\mu)$ overlay agrees to 2.8% — the
phase-function sensitivity (C1). Never write "matches the exact solution" without
that qualifier.
**Status: GUARDRAIL.**

---

## Open items

- **Advisor literature verification** (Dr. ud-Doula) — the only item that can stay
  pending at submission. Gates the novelty *wording* of A3 and Part B, nothing else.
- **Number propagation (E1 in the roadmap)** — sweeping the retired +0.16 / +0.23
  and pre-bending +0.144 / +0.201 out of every doc and figure caption. Track in
  [`next-steps.md`](../next-steps.md).
