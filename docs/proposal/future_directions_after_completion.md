# Future Directions After Core Project Completion

This document captures three optional follow-on directions discussed during planning.
None of these are part of the core REU deliverable (which is the Monte Carlo engine,
beaming function, pulse profiles, and the comparison against isotropic emission).
They are listed in the order I would attempt them if time permits after the paper
draft is written.

> **Scope rule:** none of these should be started until the core paper figures
> (validated beaming function + synthetic pulse profile + comparison) are done.
> Each one is a self-contained side project, not a dependency of the main result.

---

## 1. Speed Work — Vectorization / Numba / GPU

### Motivation

The current scalar NumPy engine runs ~10k photons/sec. That is fast enough for the
core paper (a single beaming-function curve at 200k–1M photons finishes in well
under two minutes). It becomes painful only when the project widens into:

- multi-dimensional parameter sweeps over `(τ, T_eff, composition, ...)`
- very high statistics (10⁶+ photons) for tail behavior at extreme angles
- iterative experimentation where each test costs minutes

If any of those happen post-paper, a speedup is the natural next step.

### Three tiers, in order of cost

**Tier A — Numba (recommended first attempt).** Add `@njit` to the hot functions
in `utils.py` and the per-photon inner loop in `monte_carlo.py`. Requires
refactoring the `Photon` class into plain functions over arrays (Numba does not
love Python classes). Expected speedup: 50–100×. Effort: an afternoon. Algorithm
is unchanged, so the Eddington benchmark still validates correctness directly.

**Tier B — Vectorized NumPy or PyTorch on CPU.** Treat all N photons as columns
of a tensor: positions `(N, 3)`, directions `(N, 3)`, `alive` mask `(N,)`. Replace
the per-photon `while` loop with a global `while any(alive)` loop that advances
all live photons together using boolean masks. Expected speedup: 10–50× over
scalar NumPy. Effort: a few days. Adds debugging difficulty (cannot single-step
one photon).

**Tier C — GPU (PyTorch / JAX / CuPy).** Once Tier B is done, moving tensors to
a GPU is mostly mechanical. Expected speedup: another 10–100× on top of Tier B,
depending on photon count and hardware. Worth it only at ≥10⁶ photons; below
that, GPU launch overhead dominates.

### Validation strategy

Always keep the original scalar engine as `monte_carlo_reference.py`. Any new
implementation must reproduce the Eddington `I(μ) ∝ 1 + 1.5μ` benchmark and
match the reference engine's beaming-function histogram to within Poisson noise
on a fixed seed.

### Expected scientific payoff

Modest. The core paper does not need this. The win is enabling broader parameter
sweeps for a follow-up paper or for inclusion as supplementary material.

---

## 2. Mechanistic Interpretability of a Beaming-Function Surrogate

### Motivation

This is the direction with the highest personal-development value and the most
unusual angle for a physics paper. Standard astrophysics ML uses neural networks
as black-box surrogates. Here we would instead **train a tiny network on a known
physics problem and dissect what it learned**, using the closed-form Eddington
result as ground truth.

### What it looks like

1. Run the (already-validated) Monte Carlo engine across a grid of atmospheric
   parameters: τ ∈ [0.1, 30], T_eff, possibly composition. Save `I(μ; params)`.
2. Train a small MLP (2–4 layers, width ≤ 64) to predict `I(μ)` from
   `(μ, τ, T_eff, ...)`.
3. Verify the surrogate matches MC output to within statistical noise.
4. **Interpret the network.** Concrete questions worth probing:
   - In the small-τ regime, does the network internally encode the linear
     `1 + 1.5μ` form? Can you find a direction in activation space that
     corresponds to that linearity?
   - Are there individual neurons that activate selectively at specific μ values
     or specific τ regimes?
   - Train a sparse autoencoder on the hidden activations. Do its features
     correspond to interpretable physical quantities (limb-darkening strength,
     escape probability, optical-depth regime)?
   - When the network smoothly interpolates between the optically-thin and
     optically-thick regimes, what does the transition look like in
     activation space?

### Why this is good as a mech-interp exercise

The usual mech-interp problem with LLMs is that you do not know the ground
truth — you cannot tell whether the "feature" you found is real or a
post-hoc story. Here the ground truth is a closed-form physics equation. If
you claim "neuron 7 represents limb darkening," you can falsify that claim
by checking against the Eddington formula. That is a much cleaner research
target than poking at GPT-2.

### Effort estimate

A focused two-to-four-week side project. Most of the time is in the interp,
not the training (the network is tiny). Best done in PyTorch because the
mech-interp tooling (TransformerLens-style hooks, SAE libraries) is
PyTorch-native.

### Expected scientific payoff

Probably not part of the main paper. Could become a short methods note or a
workshop paper at an ML-for-physics venue. Strong portfolio piece for graduate
applications in either ML or computational astrophysics.

---

## 3. Simulation-Based Inference (SBI) for Pulse Profiles

### Motivation

The actual inference problem NICER teams solve is the inverse of what this
project does: given an observed pulse profile, infer the underlying neutron
star parameters `(M, R, geometry, atmosphere)`. They use MCMC. Modern SBI
methods — neural posterior estimation, neural ratio estimation, normalizing
flows — solve the same problem in a way that treats the simulator as a black
box and is typically much faster at inference time once trained.

This is the closest direction to what real research groups are publishing now.

### What it looks like

1. Wrap the Monte Carlo + pulse-profile pipeline as a single function
   `simulate(params) -> pulse_profile`.
2. Generate a training set: thousands of `(params, pulse_profile)` pairs by
   sampling from a prior over parameters.
3. Train a normalizing flow (e.g., with the `sbi` Python package) to learn
   `p(params | pulse_profile)`.
4. Given a synthetic "observed" pulse profile, recover the posterior over
   parameters and compare against the truth.
5. Optional: apply to a real NICER pulse profile and compare against the
   official NICER team's published posterior.

### Effort estimate

A four-to-eight-week project. Most of the cost is in (a) the parameter sweep
to build a training set (this is where the speed work from §1 actually
matters), and (b) getting the SBI plumbing to converge cleanly.

### Expected scientific payoff

Potentially publishable as a methods paper if it works well: "neural SBI for
NICER pulse-profile inference." Real research bite. Heavier ML engineering
than the mech-interp direction; less interpretive content.

---

## 4. Convergence Analysis — How Many Photons Is Enough?

> **Scope caveat:** this is the one item on the list that probably should
> *not* wait until after the paper. The current N values in the validation
> harness (5000 for energy conservation, 1000 for mean free path, etc.) were
> picked by what felt fast, not by analysis. A convergence study justifies
> those numbers and belongs in the methods section or appendix of the core
> paper itself. I am listing it here only because the other three directions
> live in this document and it would feel odd to put a methodology check
> somewhere else.

### Motivation

Monte Carlo statistical error scales as 1/√N. Each doubling of photons cuts
noise by only ~30% while doubling runtime, so the cost-to-precision tradeoff
gets steadily worse as N grows. Eventually statistical noise drops below the
*systematic* error floor — finite μ-bin width in the beaming-function
histogram, the Eddington approximation itself, the validation tolerance you
chose — and additional photons stop buying real precision. The "right" N is
not one number; it is per observable. A total-energy conservation check
converges fast (every photon contributes). A tail μ-bin of the beaming
function (μ near 0 or near 1) converges much more slowly, because few
photons land there.

Without doing this analysis I cannot honestly answer the natural reviewer
question "why 5000 and not 50000, or 500?"

### What it looks like

1. Pick the metrics that the paper actually reports: energy conservation
   residual, mean-free-path estimate, and per-bin beaming-function values
   with special attention to the lowest-statistics bins.
2. Sweep N across decades — e.g. 1e3, 3e3, 1e4, 3e4, 1e5, 3e5, 1e6 — with
   fixed seed offsets so the runs are reproducible.
3. For each metric, plot error vs N on log-log axes. While statistics-limited
   you should see a clean -1/2 slope. The slope flattens once you hit the
   systematic floor. The bend is the knee.
4. For each observable, set production N just past its knee plus a small
   safety margin. Tail bins may demand a larger N than bulk quantities; that
   is fine and worth saying out loud.
5. Optional: if the tail bins demand absurdly large N to converge, that is
   the natural opening for **variance reduction** — importance sampling
   biased toward the under-sampled angles, with weights to keep the estimator
   unbiased. That is a small extension, not a separate project.

### Effort estimate

A day or two. The engine, the validation harness, and the plotting utilities
already exist. Most of the work is the sweep script and one or two figures.

### Expected scientific payoff

Modest but real. Concretely:
- An appendix figure ("sample-size justification") for the paper. Reviewers
  like this and it is cheap insurance against "did you converge?" comments.
- Honest, defensible N values in the methods section instead of round
  numbers chosen by feel.
- A natural lead-in to a short discussion of variance reduction if the tail
  bins motivate it.

The personal-development value is low compared to §2 and §3, but the
*scientific hygiene* value is high enough that this is the only item on this
page I would actually do before the paper is locked.

---

## Comparison

| Direction | Effort | ML content | Paper payoff | Personal-development value |
|---|---|---|---|---|
| Speed work | low (afternoon for Numba; week+ for full vectorization) | none | enables bigger sweeps, supplementary figures | low–medium |
| Mech interp | medium (2–4 weeks) | high; mech-interp-specific | short methods note / workshop paper | high (mech interp portfolio) |
| SBI | high (4–8 weeks) | high; standard ML engineering | possible methods paper | high (modern Bayesian ML) |
| Convergence study | very low (1–2 days) | none | sample-size justification in paper; cheap | low–medium (good scientific hygiene) |

## Decision rule

The convergence study is the odd one out: it should be folded into the core
paper rather than treated as post-completion work. Do it before the final
figures are locked.

After the core paper is drafted, for the other three:

- If the bottleneck has been "I cannot iterate fast enough," do **speed work** first.
- If the goal is **a portfolio piece for grad school in interpretability**, do
  the **mech interp** project. The Eddington ground truth is the unique
  advantage here that does not exist in LLM mech interp.
- If the goal is **publishable ML-for-astrophysics methods work**, do **SBI**.
  This is the direction with the most real-research traction in the NICER
  community right now.

These are not mutually exclusive. Speed work strictly helps both of the
others, so if you do more than one, do speed work first.
