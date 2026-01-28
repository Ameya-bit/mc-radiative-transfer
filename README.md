# Monte Carlo Radiative Transfer

Simulating photon transport through neutron star atmospheres to compute realistic beaming functions for NICER pulse profile modeling.

---

## Project Overview

**Goal**: Determine whether realistic angular emission patterns (beaming functions) from neutron star atmospheres significantly affect X-ray pulse profiles measured by NICER.

**Approach**: Monte Carlo simulation of photon scattering through an atmospheric slab → Extract beaming function I(μ) → Apply to synthetic pulse profile generation → Compare with standard isotropic assumption.

---

## Phase 1: Monte Carlo Engine (Current)

### Physics Model

This simulation tracks photon packets propagating through a plane-parallel atmospheric slab defined in optical depth coordinates.

```
τ = 0        ← TOP (escape surface)
   ↑
   │  photon scatters, propagates
   │
τ = τ_total  ← BOTTOM (injection point)
```

### Design Decisions

| Choice | Decision | Rationale |
|--------|----------|-----------|
| **Scattering** | Thomson | Correct phase function P(μ) ∝ (1 + μ²) for electron-dominated atmospheres |
| **Bottom boundary** | Absorb | Standard approach — photons returning downward are "lost to the thermal source" |
| **Energy** | Monochromatic | Isolates the angular redistribution effect; avoids Compton/Klein-Nishina complexity |
| **Polarization** | Not tracked | Second-order effect with high implementation cost |

### What We Gain

- **Physically meaningful**: Thomson scattering is the dominant process for X-ray photons in hot NS atmospheres
- **Tractable complexity**: ~8 functions, ~200 lines of Python
- **Validatable**: Known analytic limits exist (thin slab, semi-infinite slab)

### What We Defer (Future Work)

| Feature | Why Deferred | Impact |
|---------|--------------|--------|
| **Compton scattering** | Requires Klein-Nishina cross-section; energy-dependent opacity | Would allow spectral analysis |
| **Polarization** | Requires Stokes vector tracking; Mueller matrix algebra | Would enable polarimetric predictions |
| **Magnetic effects** | O-mode/X-mode splitting in strong B-fields | Relevant for magnetars |
| **Curved geometry** | Full spherical atmosphere instead of plane-parallel slab | Needed for very extended atmospheres |

These are noted as limitations in the paper and provide clear directions for follow-up studies.

---

## Project Structure

```
mc-radiative-transfer/
├── README.md               # This file
├── docs/
│   ├── monte_carlo_nicer.pdf   # Task list / research plan
│   └── RNAA_draft.pdf          # Paper draft
├── src/
│   └── monte_carlo.py          # Phase 1: MC engine
├── tests/
│   └── test_monte_carlo.py     # Unit tests
├── data/
│   └── (simulation outputs)
└── figures/
    └── (plots for paper)
```

---

## Phases

- [x] **Phase 1**: Monte Carlo engine (photon transport simulation)
- [ ] **Phase 2**: Beaming function extraction (histogram → I(μ) curve)
- [ ] **Phase 3**: Pulse profile synthesis (apply to NICER geometry)
- [ ] **Phase 4**: Analysis & paper completion

---

## Dependencies

```
numpy       # Random sampling, array operations
matplotlib  # Visualization
pytest      # Testing
```