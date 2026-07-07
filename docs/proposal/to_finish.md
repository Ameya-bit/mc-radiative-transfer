# To Finish — Pulse-Profile Synthesis Through Paper Completion

> **Purpose.** A self-contained brief covering everything left between here (beaming-function
> phase complete, v0.7.0) and a submitted paper. It can be picked up cold in a new chat or handed
> to a mentor. It supersedes the "pulse-profile is next phase" pointers scattered across
> `README.md`, `proposal.tex`, and the task list (`docs/monte_carlo_nicer.pdf`). The governing
> principle is: **verify the machinery on a setting everyone agrees on, *then* let differences
> carry meaning.**

> **Scope (decided 2026-06-08, with user).** This is **one paper** covering the beaming-function
> library **and** its validation **and** pulse-profile synthesis with the isotropic-vs-realistic
> comparison. Magnetic fields are deferred (out of scope; future work).

---

## 0. Where the project is

Done and validated:

- Monte Carlo transport engine (Thomson scattering, plane-parallel slab), reproducible via an
  explicit `numpy.random.Generator` (v0.6.2).
- Beaming function `I(μ)` extraction (flux→intensity, validated vs. Eddington / Chandrasekhar H),
  v0.5.1.
- τ-swept beaming **library** `data/beaming_library.npz` — `I(μ; τ)` for τ ∈ {0.1…30} (v0.6.0),
  thin-τ injection defect fixed (v0.6.1).
- Convergence study justifying photon counts; the low-μ tail is the slowest observable (v0.7.0).

What remains, at a glance: build a verified pulse-profile model that consumes the library, run the
isotropic-vs-realistic comparison, anchor it to a real star differentially, and write the paper.

---

## 1. The core idea, in one paragraph

A neutron star has a hot spot; the star spins; we measure brightness vs. rotation **phase** — the
**pulse profile**. Its strength is summarized by the **pulsed fraction**
`PF = (F_max − F_min)/(F_max + F_min)`. The observed profile is three ingredients multiplied
together: **(1)** geometry + gravitational **light bending** (gravity curves rays so we partly see
the far side, which smooths the pulse and encodes mass/radius), **(2)** the **beaming function**
`I(μ)` — how surface brightness depends on emission angle `μ = cos θ` (our physics; the alternative
to the usual *isotropic* assumption), and **(3)** the observing regime (we are monochromatic /
bolometric; NICER is energy-resolved). If we skip verification and a real-star comparison
disagrees, we cannot tell whether ingredient (2) is real physics or a bug in ingredient (1). Hence
the ladder below.

---

## 2. The verification → comparison ladder (the acceptance criteria)

Rungs A–B are **pure verification** (no new physics, just "is the machine correct?"). Rungs C–D are
the **science**, and they only mean something *because* A–B passed.

| Rung | Compare against | Target | What it proves |
|------|-----------------|--------|----------------|
| **A** | Beloborodov (2002) analytic isotropic-spot profile | < 1% | geometry + bending + integration are bug-free |
| **B** | Bogdanov et al. (2019, L26) published code-comparison test case (low spin) | their stated ~% | our machinery matches the community-standard codes |
| **C** | isotropic vs. our `I(μ; τ)`, **fixed geometry** | — (this is the result) | the beaming systematic: ΔPF and Δ(shape) |
| **D** | a real star's published geometry, **differential** ΔPF vs. NICER's uncertainty | — | the systematic is large enough for NICER to see |

**Rung A is the load-bearing verification** — it depends on no one else's code, only on a
closed-form formula. **Rung B strengthens the credibility** ("agrees with X-PSI / the
Illinois–Maryland code to ~%") but is best-effort: if the L26 test tables are hard to obtain, A
alone is sufficient to verify and B becomes a stretch goal.

---

## 3. The physics/math the module needs

All in the **Schwarzschild + slow-rotation** regime (no oblateness, no Doppler) — adequate for the
low-spin verification cases and for a differential result. Oblateness/Doppler are explicitly out of
scope (§7).

**Geometry — angle between spot normal and line of sight, as the star rotates:**

```
cos ψ(φ) = cos i · cos θ_s + sin i · sin θ_s · cos φ
```

- `i` = observer inclination, `θ_s` = spot colatitude, `φ` = rotational phase (0 → 2π).

**Light bending — Beloborodov (2002) approximation** maps the emission angle `α` (from the local
radial normal) to the geometric angle `ψ`:

```
cos α = u + (1 − u) · cos ψ ,     u ≡ R_s / R = 2GM / (R c²)   (compactness)
```

- `u` is the only relativity parameter; bigger `u` = more bending = more smoothing.
- The bending **Jacobian** `d(cos α)/d(cos ψ) = (1 − u)` is constant — this is what makes the
  approximation clean and cheap.

**Visibility:** the spot is visible while `cos α ≥ 0`, i.e. `cos ψ ≥ −u/(1 − u)`. Because of
bending this admits `ψ > 90°` — that is "seeing around" the star.

**Observed (bolometric) flux from a point spot** — proportional, slow-rotation form (see Poutanen &
Beloborodov 2006 for the full expression):

```
F(φ) ∝ (1 − u) · I(μ = cos α) · cos α        when visible, else 0
```

- The local emission direction *is* `μ = cos α`, so this is exactly where the **beaming library is
  looked up**: `I(cos α)` from `data/beaming_library.npz` at a chosen τ.
- **Isotropic baseline:** replace `I(cos α)` with a constant. The `cos α` projection factor stays;
  only the beaming term changes. That isolated swap is Rung C.

**A finite spot** (task list §4 allows one or two circular spots): tile the spot into surface
elements, apply the above per element with each element's own `ψ`, sum. Start with a **point spot**
(closed-form, validates against Rung A), then optionally extend to a finite cap.

---

## 4. Milestones (versioning is a suggestion — pick the next free version)

- **v0.8.0 — Pulse-profile machinery + Rung A.** New module `src/mcrt/pulse.py`
  (geometry, Beloborodov bending, visibility, point-spot flux, `pulsed_fraction`). Unit-tested,
  reproduces the Beloborodov analytic isotropic profile to < 1%. Deep dive + README entry.
- **v0.8.1 — Rung B (best-effort).** Reproduce a low-spin Bogdanov (2019, L26) test case to its
  stated tolerance; record the agreement. If tables unavailable, document as attempted and lean on
  Rung A.
- **v0.9.0 — Rung C (the core result).** Isotropic vs. `I(μ; τ)` at fixed geometry, swept over a
  few `(θ_s, i)` operating points and over τ. Headline figure: ΔPF (and profile-shape change) vs.
  geometry/τ. **This is the figure the paper is built on.**
- **v0.9.1 — Rung D (real-star anchor).** Use PSR J0030+0451's published geometry (Riley 2019 /
  Miller 2019) as a realistic operating point; report the *differential* ΔPF and compare its size to
  NICER's reported PF / M–R uncertainty. **Relative change only — not a fit (see §6).**
- **v1.0.0 — Paper.** Reframe the draft, assemble figures, write methods incl. the convergence
  appendix and the verification ladder, submit. Citations live in `docs/paper/references.md`; the
  paper-fix and physics-caveat checklist is **§10** below. Key items:
  - Drop "Lucy 1999 photon-packet" framing → cite Whitney 2011 / Cashwell & Everett 1959 (the code
    is a single-photon random walk).
  - Remove "magnetic anisotropy (Ho & Lai 2001)" — the engine has no magnetic physics (deferred).
  - Remove the fabricated "∼15% pulsed-fraction" number from the draft; replace with the **actual**
    ΔPF measured in Rung C/D.
  - Make the Chandrasekhar H-function comparison the centerpiece of the validation section.
  - Fix the draft's "ApJL 887 L25" → **L26** for the code-comparison claim (L25 = Paper I, the data set).
  - Address the H-function caveats (§10) at the chosen depth (an open decision, §10).

---

## 5. Reuse (don't reimplement)

- Beaming library + lookup: `data/beaming_library.npz` (`tau_values`, `mu_centers`,
  `intensity_by_tau`, `b_of_tau`) — interpolate `I(μ)` at `μ = cos α`.
- Library/sweep + figure-style pattern: `scripts/tau_sweep.py`.
- Reproducible RNG plumbing (if any stochastic sampling is added): `Simulation(rng=...)`,
  `numpy.random.default_rng` / `SeedSequence` — same pattern as the convergence study.
- Pure-helper + unit-test discipline: mirror `src/mcrt/convergence.py` /
  `tests/test_convergence.py` (keep geometry/bending helpers pure, test them directly).
- Progress-log + deep-dive workflow: the "How to update the progress log" section at the bottom of
  `README.md`.

---

## 6. Two tensions in the existing docs this phase must settle

**(1) The draft paper claims physics the code does not do.** `docs/RNAA_draft.pdf` attributes the
deviation from isotropy to **magnetic anisotropy (Ho & Lai 2001)** and reports a "∼15%" pulsed-
fraction change that does not yet exist. The engine is **Thomson scattering only**; magnetic
effects are deferred. The verifiable, defensible result is **scattering-induced limb darkening**
(Rung C), which is also the only version with an analytic benchmark (Rungs A/B). **Action:** when
writing v1.0.0, reframe around limb darkening from photon transport and substitute the real ΔPF.
This is itemized in the paper-fix checklist (§10).

**(2) Two planning docs disagree on the real-star comparison.**
`week_remaining_convergence_study.md` (line ~111) calls real-NICER comparison a "later stretch …
needs GR ray-tracing tuned to a source." The task list (`docs/monte_carlo_nicer.pdf` §5) treats it
as core but says "**emphasize relative changes rather than full parameter fitting.**" Rung D follows
the task list:

- **In scope (Rung D):** fixed published geometry, swap beaming, report the *differential* ΔPF, and
  compare its magnitude to NICER's error bar.
- **Out of scope (the "later stretch"):** fitting real data / tuning a ray-tracer to recover M–R.
  That is essentially re-running the full NICER inference and is a separate project.

The caution in the week-doc was about the *hard* version; Rung D is the *light* version §5 asks for.
This document is the tie-breaker: real-star comparison is **differential, not a fit.**

---

## 7. Explicitly NOT in scope

- **Oblateness & Doppler** (oblate-Schwarzschild / special-relativistic boosting). Slow-rotation
  Schwarzschild is sufficient for the low-spin verification cases and a differential result; note as
  a limitation. (Morsink et al. 2007; Psaltis & Özel 2014 if a reviewer asks.)
- **Energy-resolved / spectral** pulse profiles. We are monochromatic/bolometric; NICER is
  energy-dependent. State as a limitation; compare bolometric shape/PF only.
- **Full parameter fitting / M–R inference** from real data (see §6).
- **Magnetic / polarization physics** (already deferred project-wide; see README "What we defer").
- The three post-paper directions in `future_directions_after_completion.md` (speed work, mech
  interp, SBI) — only after v1.0.0.

---

## 8. Verification / acceptance checklist

1. `PYTHONPATH=src python3 -m pytest -q` — existing suite stays green; new pure geometry/bending
   helpers have direct unit tests.
2. **Rung A:** point-spot, isotropic, slow-rotation profile matches the Beloborodov (2002) closed
   form to < 1% across phase, for at least two `(θ_s, i, u)` settings including one with
   `ψ_max > 90°` (visible "around the back").
3. **Rung B (best-effort):** a low-spin Bogdanov (2019, L26) case reproduced to its stated
   tolerance, or documented as attempted.
4. **Rung C:** isotropic and `I(μ; τ)` runs share *identical* geometry code (only the brightness
   term differs); ΔPF reported with the sign/trend explained physically (limb darkening sharpens or
   softens the pulse — say which and why).
5. **Rung D:** differential ΔPF at J0030 geometry stated alongside NICER's quoted uncertainty, with
   the monochromatic/model-dependent-geometry caveats written out.
6. Figures force-added under `data/` (the `data/*` PNG-tracking convention from prior versions).

---

## 9. Required reading for this phase

- **Beloborodov (2002), ApJ, 566, L85** — the `cos α = u + (1−u)cos ψ` bending approximation (Rung A).
- **Poutanen & Beloborodov (2006), MNRAS, 373, 836** — full bolometric pulse-profile flux formula.
- **Bogdanov et al. (2019), ApJL, 887, L26** — pulse-profile **code-comparison** test suite (Rung B).
  *(Note: the task list already lists Bogdanov L25 = paper I; L26 = paper II is the one needed here.)*
- **Riley et al. (2019), ApJL, 887, L21** and **Miller et al. (2019), ApJL, 887, L24** — J0030
  published geometry/parameters and quoted uncertainties (Rung D).
- *(Optional, only if oblateness comes up)* Morsink et al. (2007), ApJ, 663, 1244; Psaltis & Özel
  (2014), ApJ, 792, 87.
- For paper-writing: citations are in **`docs/paper/references.md`**; the draft-fix and physics-caveat
  checklist is folded into **§10** below.

---

## 10. Paper grounding (folded in from the old grounding doc)

Paper-prose / citation / caveat work — no code is wrong here. Citations: `docs/paper/references.md`.

**Draft fixes (the current `RNAA_draft.pdf` mis-describes the code):**

| Draft claim | Reality | Action |
|---|---|---|
| "photon-packet Monte Carlo (Lucy 1999)" | single-photon random walk, not Lucy packets | cite Whitney 2011 / Cashwell & Everett 1959 |
| "magnetic anisotropy (Ho & Lai 2001)" | no magnetic physics in the code | remove (deferred to future work) |
| "Schwarzschild light bending" pulse profiles | done via Beloborodov (2002), Rungs A/B | report the real results |
| "pulsed fractions modified by ∼15%" | placeholder, no such result yet | replace with the actual ΔPF from Rung C/D |
| (Chandrasekhar H not mentioned) | the single strongest validation | make it the centerpiece |
| cites "ApJL 887 L25" for code comparison | L25 = Paper I (data set) | change to **L26** (Paper II) |

**Physics caveats to state honestly (idealizations, not bugs):**

1. **Scattering law ≠ benchmark law.** The engine uses Rayleigh/Thomson ¾(1+μ²);
   the H-function compared against is the *scalar isotropic* (p=1) solution.
   Shapes agree closely but it isn't strictly the same physics. *Optional fix:*
   add an isotropic-scattering mode and show it reproduces scalar H(μ) to within
   noise (apples-to-apples), then present the Thomson run as the physics case.
2. **H is the semi-infinite (τ→∞) limit; the slabs are finite.** The rigorous
   finite-τ benchmark is Chandrasekhar's **X- and Y-functions** — which is why
   τ=10 matches H but τ≲1 deviates. *Minimum:* state this; *optional:* implement X/Y.
3. **Idealizations:** grey/monochromatic (no frequency dependence), conservative
   (no true absorption; the base just removes photons), unpolarized. Real NS
   atmospheres (NSX, McPhac, Suleimanov) include all three.

**Open decisions:**
- **Venue:** RNAA short note (matches current scope, editor-reviewed) vs. short
  ApJL/MNRAS paper (peer-reviewed, more grounding/figures) — drives reference count and rigor.
- **Caveat depth:** written caveats only, or also implement the isotropic-mode
  validation (caveat 1) and/or the X/Y functions (caveat 2).
