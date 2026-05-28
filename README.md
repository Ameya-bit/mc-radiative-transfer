# Monte Carlo Radiative Transfer

We simulate X-ray photons bouncing through the thin layer of plasma that sits on top of a
neutron star. Each photon is followed one scatter at a time until it either escapes into space
or is reabsorbed by the surface. By recording the directions of the photons that escape, we
measure the star's **beaming function** `I(μ)` — how its brightness depends on viewing angle.
That beaming function is the input NASA's NICER telescope needs to turn a pulsing X-ray signal
into a measurement of a neutron star's mass and radius. Most models assume the surface glows
equally in all directions; this project tests how wrong that assumption is.

> **For the full math** behind every step, each progress-log entry links to a matching
> **deep dive** in [`docs/deep-dives/`](docs/deep-dives/) — plain-language, figure-by-figure
> derivations of that week's physics.

---

## Progress Log

*Newest first. Each entry: a one-line headline, why it matters, a figure, and the technical
details tucked underneath. Week numbers follow the 10-week project plan.*

### Week 5 — The engine is physically validated
*2026-03-14 · commits `f21738d`–`a3abc18`*

**The simulation now reproduces a textbook result: the escaping radiation is "limb-darkened"
exactly the way classical theory predicts.**

This is the milestone that says the engine can be trusted. We let photons scatter through a
thick atmosphere, histogrammed the angles at which they escaped, and compared that to the
century-old Eddington prediction that a scattering atmosphere should look brightest face-on and
dimmest at the edge (the same effect that makes the edge of the Sun look dim). The Monte Carlo
result lands right on the predicted `1 + 1.5μ` curve. Two independent sanity checks also pass,
so we know no photons are being silently lost or created.

![Beaming function: Monte Carlo result vs. the Eddington 1 + 1.5μ prediction](data/beaming_function.png)

📐 **Full derivation:** [Week 5 — Validation & the Beaming Function](docs/deep-dives/week-5-validation.md)

<details>
<summary>Technical details</summary>

- **Beaming function** `I(μ)`: histogram of escape angles `μ = cos θ` for a thick slab
  (τ_total = 10), overlaid on the semi-infinite Eddington limit `I(μ) ∝ 1 + 1.5μ` — they match.
- **Energy/photon conservation:** every injected photon ends as either *escaped* or *absorbed*;
  5000/5000 accounted for, exactly.
- **Mean free path:** total path length ÷ total scatters ≈ **1.0** optical depth (measured
  ~1.00–1.02), confirming the `−ln(U)` step sampling is correct.
- Code: `scripts/validate_engine.py`.
</details>

**Next:** extract beaming functions across a range of atmosphere thicknesses (Weeks 6–7).

---

### Weeks 3–4 — Photons travel through the atmosphere end-to-end
*2026-03-14 · commit `f21738d`*

**A photon can now be injected at the base of the atmosphere, scatter its way through, and
either escape or be reabsorbed — the complete random walk.**

This is the core of the project. Photons start at the bottom moving in random upward
directions, take exponentially-distributed steps, and at each stop scatter off an electron via
Thomson scattering (which slightly prefers forward/backward over sideways). When a photon
reaches the top it escapes and we record its exit angle; if it drifts back down to the surface
it is absorbed. Everything is measured in *optical depth* rather than meters, which keeps the
physics general.

![Photon random walks through atmospheres of increasing thickness](docs/deep-dives/figures/07_random_walks.png)

📐 **Full derivation:** [Weeks 3–4 — Photon Transport](docs/deep-dives/weeks-3-4-photon-transport.md)

<details>
<summary>Technical details</summary>

- **Geometry:** 3D Cartesian in optical-depth coordinates; τ = 0 is the top (escape), τ =
  τ_total is the bottom (injection). Plane-parallel slab.
- **Injection:** isotropic over the upward hemisphere (uniform in cos θ, not in θ — see the
  [Weeks 3–4 deep dive](docs/deep-dives/weeks-3-4-photon-transport.md)).
- **Transport:** step size `Δτ = −ln(U)` (exponential free path); Thomson phase function
  `(3/4)(1 + μ²)` via rejection sampling; 3D direction update via `rotate_vector`.
- **Boundaries:** escape at τ ≤ 0 (record exit μ), absorb at τ ≥ τ_total.
- Code: `src/mcrt/monte_carlo.py` (`Photon`, `Simulation`).
</details>

**Next:** validate the random walk against analytic limits (Week 5).

---

### Weeks 1–2 — The building blocks are in place and tested
*2026-01-24 → 2026-03-14 · commits `ab83703`, `f50cd2a`*

**Every random-sampling primitive the simulation depends on is written and unit-tested.**

Before tracking a single photon we needed the small mathematical tools that the random walk is
built from: how far a photon travels before it scatters, which direction it scatters into, and
how to point it correctly in 3D afterward. Each of these is a short, independently-tested
function, so when the full engine was assembled in Weeks 3–4 we already trusted its parts.

![Sampling primitives: the Thomson phase function](docs/deep-dives/figures/04_thomson_phase.png)

📐 **Full derivation:** [Weeks 1–2 — Sampling Primitives](docs/deep-dives/weeks-1-2-sampling-primitives.md)

<details>
<summary>Technical details</summary>

- `sample_step_size()` — exponential free path via `−ln(U)`.
- `sample_thomson_angle()` — Thomson `(3/4)(1 + μ²)` by rejection sampling.
- `get_random_direction()`, `rotate_vector()` — isotropic directions and scatter rotation.
- pytest suite in `tests/test_physics.py` checks the exponential mean (≈ 1.0), angle bounds,
  and unit-norm preservation under repeated scatters.
- Code: `src/mcrt/utils.py`.
</details>

---

## Reference

### Physics model

A plane-parallel atmospheric slab in optical-depth coordinates:

```
τ = 0        ← TOP (escape surface)
   ↑
   │  photon scatters, propagates
   │
τ = τ_total  ← BOTTOM (injection point)
```

### Design decisions

| Choice | Decision | Rationale |
|--------|----------|-----------|
| **Scattering** | Thomson | Correct phase function P(μ) ∝ (1 + μ²) for electron-dominated atmospheres |
| **Bottom boundary** | Absorb | Standard approach — photons returning downward are "lost to the thermal source" |
| **Energy** | Monochromatic | Isolates the angular redistribution effect; avoids Compton/Klein-Nishina complexity |
| **Polarization** | Not tracked | Second-order effect with high implementation cost |

### What we defer (future work)

| Feature | Why deferred | Impact |
|---------|--------------|--------|
| **Compton scattering** | Requires Klein-Nishina cross-section; energy-dependent opacity | Would allow spectral analysis |
| **Polarization** | Requires Stokes vector tracking; Mueller matrix algebra | Would enable polarimetric predictions |
| **Magnetic effects** | O-mode/X-mode splitting in strong B-fields | Relevant for magnetars |
| **Curved geometry** | Full spherical atmosphere instead of plane-parallel slab | Needed for very extended atmospheres |

These are noted as limitations in the paper and provide clear directions for follow-up studies.

### Project structure

```
mc-radiative-transfer/
├── README.md                  # This file (overview + progress log)
├── pyproject.toml             # Package metadata (installable as `mcrt`)
├── requirements.txt           # numpy, matplotlib, pytest
├── src/
│   └── mcrt/                  # The simulation package
│       ├── __init__.py
│       ├── monte_carlo.py     # Photon + Simulation engine
│       └── utils.py           # Sampling & geometry primitives
├── scripts/                   # Runnable entry points
│   ├── validate_engine.py     # Validation + beaming-function plot
│   └── plot_paths.py          # 3D random-walk visualization
├── tests/
│   ├── conftest.py            # Makes `mcrt` importable without install
│   └── test_physics.py        # Unit tests for the primitives
├── docs/
│   ├── deep-dives/            # Per-entry math deep dives
│   │   ├── weeks-1-2-sampling-primitives.md
│   │   ├── weeks-3-4-photon-transport.md
│   │   ├── week-5-validation.md
│   │   ├── make_figures.py    # Regenerates the figures below
│   │   └── figures/           # Explanatory figures (01–08)
│   ├── monte_carlo_nicer.pdf  # Task list / research plan
│   ├── RNAA_draft.pdf         # Paper draft
│   └── proposal/              # Proposal + future directions
└── data/                      # Simulation outputs (plots + Week 6+ data)
```

### Setup & usage

```bash
pip install -e .              # makes `mcrt` importable everywhere
python scripts/validate_engine.py   # validation + beaming-function plot
python scripts/plot_paths.py        # random-walk visualization
pytest                              # run the unit tests
```

### Timeline

- [x] **Weeks 1-2**: Physics setup & environment (Thomson utilities, testing framework)
- [x] **Weeks 3-4**: Monte Carlo engine (photon transport, boundary handling)
- [x] **Week 5**: Validation & benchmarking (energy conservation, MFP, beaming function)
- [ ] **Weeks 6-7**: Beaming function extraction across τ_total values
- [ ] **Weeks 8-9**: Pulse profile synthesis (apply to NICER geometry)
- [ ] **Phase 4**: Analysis & paper completion

---

## How to update the progress log

Each week produces **two** linked pieces: a short entry at the **top** of the Progress Log, and
a companion **deep dive** in [`docs/deep-dives/`](docs/deep-dives/) holding the full math. Keep
the entry headline plain-physics (no code jargon); put code in `<details>` and the derivations
in the deep dive.

```markdown
### Week N — <plain-physics headline read in 30 seconds>
*<date>* · commits `abcd123`–`efgh456`

**<One bold sentence: what now works and why it matters for the science.>**

<2–4 sentences of physics context — what was added, what we can now do that we couldn't.>

![caption](docs/deep-dives/figures/NN_figure.png)

📐 **Full derivation:** [Week N — <title>](docs/deep-dives/week-N-<slug>.md)

<details>
<summary>Technical details</summary>

- Method / equation added (with the physics)
- What was validated + the numeric result
- Decisions and rationale
- Code: files touched
</details>

**Next:** <one line>
```

For the companion deep dive, copy an existing file in `docs/deep-dives/` as a template: open
with a link back to its progress-log entry, a **Builds on:** line pointing to the prior deep
dive, the derivation with figures (regenerate via `make_figures.py`), and a closing **Next:**.

Workflow: when a week's work is committed, both pieces are drafted from that week's commits and
figures, then the physics framing is reviewed before they land.
