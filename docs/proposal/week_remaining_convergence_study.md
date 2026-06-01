# This Week — Remaining Work: Convergence Study

> **Purpose.** A self-contained brief so this can be picked up cold in a new chat. It covers the
> one remaining item of the beaming-function phase (proposal Weeks 6–7) before pulse-profile
> synthesis (Weeks 8–9). The τ-sweep beaming-function library is **done** (v0.6.0) and its
> thin-τ injection defect is **fixed** (v0.6.1); what's left is the convergence study and the
> small seeding prerequisite it needs.

---

## Context — where the project is

- The Monte Carlo engine, the beaming function `I(μ)`, and the τ-swept **library** `I(μ; τ)`
  (`data/beaming_library.npz`) are all implemented and validated against Eddington / Chandrasekhar.
- Photon counts in the harness (5000, 1000, 200000) were chosen **by feel**, not by analysis.
  A convergence study justifies them and belongs in the paper's methods/appendix — see
  `docs/proposal/future_directions_after_completion.md §4` ("Convergence Analysis — How Many
  Photons Is Enough?"), the canonical description of this task.
- After this, the library feeds **pulse-profile synthesis** (the next phase, not this week).

The existing `scripts/convergence_study.py` is only a **thin** version: it tracks the fitted slope
`b` at a single τ=10 across photon counts. It is **not** the per-observable "error vs N on log-log,
find the knee" study the paper needs. The task is to replace/extend it.

---

## Task 1 — Reproducible seeding in the engine (prerequisite)

`src/mcrt/monte_carlo.py` `Simulation` currently calls global `np.random` with no seed, so runs
are not reproducible. The convergence study needs fixed, independent seed offsets per run.

- Add an optional `rng` (or `seed`) parameter to `Simulation.__init__`, defaulting to
  `np.random.default_rng()`, and route all sampling through it.
- This touches the sampling primitives in `src/mcrt/utils.py` (`sample_step_size`,
  `sample_thomson_angle`, `get_random_direction`) and the injection sampling in
  `monte_carlo.py` (`phi`, `costheta = sqrt(U)`). Thread the generator through, or have the
  primitives accept an `rng` argument. **Keep existing call sites working** (default → current
  behavior).
- **Test (TDD):** two `Simulation`s with the same seed produce identical `escaped_mu`; different
  seeds differ. Add to `tests/test_physics.py`.

> Note: `scripts/tau_sweep.py` currently uses a one-off global `np.random.seed(SEED)` for
> reproducibility. Once `Simulation` takes an `rng`, optionally migrate it to pass an explicit
> generator (low priority).

---

## Task 2 — The convergence study (the deliverable)

Extend `scripts/convergence_study.py` (keep the old slope-vs-N plot as one panel).

**Observables the paper reports:**
1. **Energy-conservation residual** — `(escaped + absorbed − injected)`; converges fast (every
   photon contributes).
2. **Mean-free-path estimate** — `Σ path / Σ scatters` vs the expected 1.0.
3. **Per-μ-bin beaming-function values** — with explicit attention to the **low-μ tail bins**
   (μ → 0), which converge slowest because few photons land there. (This is the same tail that
   makes τ = 30 noisy in the library — see `docs/deep-dives/v0.6.1-isotropic-injection.md`.)

**Method:**
- Sweep `N ∈ {1e3, 3e3, 1e4, 3e4, 1e5, 3e5, 1e6}` with fixed seed offsets (enabled by Task 1).
- There is no closed-form "truth" for the tail bins, so estimate the **statistical error** by
  running each N over several seeds and taking the spread (std across seeds).
- Plot **error vs N on log-log** per observable. While statistics-limited the slope is ≈ **−1/2**;
  it flattens once the systematic floor is hit. **The bend is the knee.**
- For each observable, recommend a production `N` just past its knee plus a small margin. Note
  explicitly that tail bins may demand larger N than bulk quantities (and that the τ ≥ 10 / low-μ
  corner is the binding case). Record these as the defensible values that replace the by-feel
  5000 / 1000 / 200000 in `scripts/validate_engine.py`.
- Put the knee-finding / per-bin-error logic in a small **pure helper** (e.g. `bin_error_vs_n`)
  with unit tests, not buried in `__main__`.

**Figures:** `data/convergence_*.png` (force-add like the other figures; `data/*` is gitignored
but result PNGs are tracked).

---

## Task 3 — Documentation

- New deep dive `docs/deep-dives/v0.7.0-convergence-study.md` (versioning is a suggestion — pick
  the next free version): the −1/2 slope argument, the per-observable knees, the recommended N
  values, and the tail-bin caveat. This is the "sample-size justification" appendix figure.
- Add a v0.7.0 progress-log entry at the top of `README.md` (newest first), update the project
  structure if new files are added, and bump `pyproject.toml`. Follow the "How to update the
  progress log" section at the bottom of `README.md`. **Include commit hashes** in the entry's
  date line (`· commits \`x\`–\`y\``), matching every other entry.

---

## Reuse (don't reimplement)

- `Simulation` — `src/mcrt/monte_carlo.py`
- `extract_intensity`, `fit_limb_darkening_slope` — `src/mcrt/beaming.py`
- Energy-conservation and mean-free-path logic — `scripts/validate_engine.py`
  (`validate_energy_conservation`, `validate_mean_free_path`)
- Library / sweep pattern and figure style — `scripts/tau_sweep.py`

## Verification

1. `PYTHONPATH=src python3 -m pytest -q` — existing 9 tests stay green + new seeding and
   convergence-helper tests pass.
2. `PYTHONPATH=src python3 scripts/convergence_study.py` — log-log error-vs-N panels show a ≈ −1/2
   slope that flattens (a visible knee) for the bulk observables; tail bins may not have hit their
   floor yet — note this rather than hide it.
3. Sanity: a τ = 10 point still reproduces the validated v0.5.1 / v0.6.1 slope (~1.7–1.8).

## Explicitly NOT this week

- Pulse-profile synthesis (rotating NS + hot spot + ray-tracing) — next phase; consumes
  `data/beaming_library.npz`.
- Real NICER-data comparison — documented later stretch (needs GR ray-tracing tuned to a source).
  The core comparison remains realistic-`I(μ)` vs. isotropic emission.
