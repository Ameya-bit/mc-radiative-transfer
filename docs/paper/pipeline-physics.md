# Pipeline & Physics — the equation chain

> The whole codebase, written as **one continuous equation and unfolded one
> substitution at a time** — from the paper's number at the top down to the raw
> random numbers at the bottom. Every numbered equation names the function that
> implements it, links the reference it comes from
> ([`references.md`](references.md)), and points at the check that verifies it.
> The next link — which paper claim each result supports — is
> [`claims-evidence.md`](claims-evidence.md). Derivations and result discussion
> live in `docs/deep-dives/`.
>
> *Supersedes the pipeline-diagram version of this file (2026-07-07). All
> numerical results quoted here are the v0.9.8–v0.9.12 production values.*

## The whole thing, one line

The paper's central number is

$$\Delta \mathrm{PF} \;=\; \mathrm{PF}\big[F_\star^{\text{realistic}}\big] \;-\; \mathrm{PF}\big[F_\star^{\text{isotropic}}\big],$$

and the single-spot flux that everything is built from, fully substituted, reads

$$F(\varphi)\;\propto\;\underbrace{\gamma\,\delta^{\,n}}_{\texttt{rotating.py}}\;\cdot\;\underbrace{D(\psi)}_{\texttt{bending.py}}\;\cdot\;\underbrace{I\big(\,\delta\cos\alpha\,;\;\tau\big)}_{\texttt{beaming.py}\;\leftarrow\;\texttt{monte\_carlo.py}}\;\cdot\;\underbrace{\delta\cos\alpha}_{\text{projection}}\;\cdot\;\mathbf{1}\!\left[\cos\alpha\ge 0\right],\qquad \cos\alpha=\mathrm{bend}\!\big(\cos\psi(\varphi)\big)$$

The rest of this document unfolds that expression symbol by symbol, equations
**(1)–(23)**. The isotropic-vs-realistic comparison changes exactly **one
factor**: $I \equiv 1$ versus the Monte-Carlo-built $I(\mu;\tau)$ — every other
symbol is held identical between the two runs.

## How the files chain

- **Building $I(\mu;\tau)$ (run once, ~40 min, Monte Carlo):**
  `src/mcrt/utils.py` (random samplers) → `src/mcrt/monte_carlo.py` (photon walk)
  → `scripts/tau_sweep.py` (sweep + pooling) → `data/beaming_library.npz`
  → `src/mcrt/beaming.py` (`beaming_lookup` turns a saved row into the callable $I(\mu)$).
- **Using it (milliseconds, deterministic):**
  `src/mcrt/pulse.py` (one spot's $F(\varphi)$, importing `bending.py` and
  `rotating.py`) → `scripts/anchor_lib.py` (sum spots into a star) →
  the experiment scripts (each one runs the chain twice and reports $\Delta\mathrm{PF}$).

## Notation

| symbol | meaning | defined by |
|---|---|---|
| $\varphi$ | rotation phase over one spin, $[0, 2\pi)$ | uniform grid of $n_\text{phase} = 1024$ points, `pulse.compute_profile` |
| $i,\ \theta_s$ | observer inclination (angle between spin axis and line of sight); spot colatitude (angle between spin axis and spot center) | inputs (published fits) |
| $u$ | compactness $2GM/(Rc^2)$ = Schwarzschild radius / stellar radius — how deep the star sits in its own gravity (some papers, incl. the NICER fits in [`references.md`](references.md), call $GM/Rc^2 = u/2$ the "compactness") | input; the neutron-star anchors have $u \approx 0.3\text{–}0.5$ |
| $\psi$ | angle between spot normal and line of sight | eq. (5) |
| $\alpha$ | emission angle at the surface (from the local vertical) | eqs. (6)/(7) |
| $\mu,\ \mu'$ | $\cos\alpha$; its aberrated value in the spot's own rest frame | eq. (12) |
| $D$ | lensing Jacobian $d\cos\alpha/d\cos\psi$ — solid-angle stretch of the bending map | eqs. (6)/(7) |
| $F,\ F_\star$ | one spot's flux; the whole star's — both in arbitrary units, only ratios matter | eqs. (4), (3) |
| $w_k,\ \varphi_{0,k},\ A_k,\ T_k$ | spot $k$'s brightness weight, longitude, area, temperature | eq. (3) |
| $\nu,\ R,\ c$ | spin frequency (Hz, as measured far from the star); stellar radius; speed of light | inputs to `Rotation(ν, R)` |
| $\beta,\ \gamma,\ \xi,\ \delta$ | spot speed as a fraction of $c$; Lorentz factor; ray–velocity angle; Doppler factor | eqs. (8)–(11) |
| $B,\ n$ | the Doppler flux boost; its exponent (4 = energy flux, 3 = photon counts) | eqs. (13)/(14) |
| $\tau$ | optical depth of the scattering slab: its thickness in units of photon mean free paths | input to the Monte Carlo (Part IV); downstream reads the $\tau = 10$ row |
| $I(\mu;\tau)$ | the beaming function: surface brightness vs. viewing angle | eqs. (15)–(22) |
| $N_j,\ \bar\mu_j$ | escaped-photon count in $\mu$-bin $j$; that bin's mean escaped $\mu$ | eqs. (16)/(17) |
| $a,\ b$ | fitted brightness scale (intercept) and limb-darkening slope of $I = a(1+b\mu)$ | eq. (18) |
| $U$ | a uniform random number on $(0,1)$ | seeded `numpy` Generator |
| $\mathbf{1}[\,\cdot\,]$ | indicator: 1 where the bracketed condition holds, 0 elsewhere | — |
| $\mathrm{PF}$ | pulsed fraction | eq. (2) |
| $A_0, A_1, A_2$ | Fourier amplitudes of the pulse $F(\varphi)$: mean, fundamental, second harmonic | see [Quantities that live between equations](#quantities-that-live-between-equations) |

## Project shorthand

Recurring names and jargon, defined once:

- **Anchors** — the two real pulsars, PSR J0030+0451 ("J0030") and PSR J0740+6620
  ("J0740"), whose published spot geometries are hard-coded as test cases in
  `scripts/anchor_lib.py`.
- **Riley / Miller** — the two independent analysis teams that each published a
  fit to the same data for each anchor; the pipeline runs both fits.
- **NICER** — the X-ray telescope (on the International Space Station) whose data
  those fits come from; its band is roughly 0.2–2 keV.
- **SD1a / SD1c** — standard test problems from the
  [Bogdanov et al. (2019)](references.md#bogdanov-et-al-2019--paper-ii-l26)
  code-comparison suite, with reference waveforms computed by independent groups:
  the same star and spot spun at 1 Hz (SD1a — velocity effects negligible) and
  200 Hz (SD1c). Matching them certifies the pulse machinery.
- **PB06** — shorthand for
  [Poutanen & Beloborodov (2006)](references.md#poutanen--beloborodov-2006).
- **Gate G1** — this project's pre-registered decision rule: switch the anchors
  to exact bending (7) if the measured linear-map bias mattered. It fired (v0.9.9).
- **Gate G2** — the companion pre-registered rule for the rotation layer (Part III):
  report the anchors with spin (Doppler + aberration) included if it moved the
  result. It fired at 37σ (v0.9.10), which is why the headline is quoted both static
  and spinning.
- **Grey** — treating all photon energies the same (no color information); the
  Monte Carlo and its library are grey.
- **Bolometric** — summed over all photon energies, as opposed to a telescope band.
- **Limb darkening** — a glowing surface looking dimmer toward its edge (grazing
  viewing angles), like the Sun's dimmer rim; its strength here is the slope $b$.
- **Eddington law / $H(\mu)$** — two classical predictions for the emergent
  brightness of an infinitely deep scattering atmosphere:
  [Eddington's](references.md#eddington-limb-darkening) approximate
  $I \propto 1 + 1.5\mu$, and
  [Chandrasekhar's](references.md#chandrasekhar-1960) $H$-function, which is
  exact when the scattering is isotropic. They bracket our thick-slab result
  from below and above.
- **Seed** — one independent rerun of the Monte Carlo with a different
  random-number stream; the spread across 5 seeds gives the error bars.
- **Reduced $\chi^2$** — goodness-of-fit statistic; ≈ 1 means the curves agree
  to within the error bars.
- **Tiling** — whether a star's spots between them keep some flux visible at
  every phase; see [Quantities that live between equations](#quantities-that-live-between-equations).

---

## Part I — the paper's numbers

**(1) The headline difference** — every experiment script

$$\Delta \mathrm{PF} = \mathrm{PF}\big[F_\star^{\text{realistic}}\big] - \mathrm{PF}\big[F_\star^{\text{isotropic}}\big]$$

Run the identical chain twice — once with `beaming=None` ($I \equiv 1$), once
with the library row — and subtract; positive means realistic limb-darkened
beaming makes the star look *more* pulsed than the isotropic assumption predicts.
*Headline:* J0740 static $+0.137 \pm 0.003$ (Riley) / $+0.195 \pm 0.005$ (Miller); with spin and the NICER band included, $+0.02\ldots+0.06$.

*where the pulsed fraction of any light curve is:*

**(2) Pulsed fraction** — `pulsed_fraction()` · `src/mcrt/pulse.py`

$$\mathrm{PF} = \frac{F_{\max} - F_{\min}}{F_{\max} + F_{\min}}$$

The standard pulse-strength summary; note it **saturates at 1** whenever the flux
touches zero at any phase, which is what makes multi-spot *tiling* (below) decisive.
*Check:* closed form (23); `tests/test_pulse.py`.

*where the star's flux is the weighted sum of its spots:*

**(3) A star = weighted spots** — `multi_spot_flux()` · `scripts/anchor_lib.py`

$$F_\star(\varphi) = \sum_k w_k\, F\!\big(\varphi - \varphi_{0,k}\big), \qquad w_k \propto A_k\, T_k^4$$

Light is additive, so each spot is one full run of the chain below, shifted to its
longitude $\varphi_{0,k}$ (`np.roll` by $\mathrm{round}(\varphi_{0,k} \cdot n_\text{phase})$
grid points) and weighted by how much light it puts out: its area $A_k$ times the
fourth power of its temperature $T_k$ — a hot surface radiates as $T^4$ (the
Stefan–Boltzmann law), so a slightly hotter spot outshines a cooler one by a lot.
Only relative weights within one star matter — never compare $w_k$ across anchors
(Miller quotes temperatures in keV, Riley in kelvin).
*Check:* PB06's visibility classes for antipodal (exactly opposite) spot pairs, `validate_phase_diagram.py`.

---

## Part II — one spot's flux

**(4) The master flux equation** — `point_spot_flux()` · `src/mcrt/pulse.py`

$$F(\varphi) \;\propto\; \underbrace{\gamma}_{(9)}\;\underbrace{B(\delta)}_{(13)\,\text{or}\,(14)}\;\underbrace{D}_{(6)\,\text{or}\,(7)}\;\underbrace{I(\mu')}_{(15)}\;\underbrace{\mu'}_{(12)}\;\cdot\;\mathbf{1}\!\left[\cos\alpha \ge 0\right]$$

Each factor is one physical effect — comoving area ($\gamma$), Doppler boost
($B$), gravitational solid-angle stretch ($D$), surface brightness ($I$), and
projection ($\mu'$) — and the flux is zero wherever the spot has set. Constants
(redshift powers, area, distance) drop out of shapes and PF, so only ratios are
physical. On a **frozen star** (`rotation=None`, the verified baseline)
$\gamma = B = 1$ and $\mu' = \cos\alpha$, leaving $F \propto D\, I(\cos\alpha)\cos\alpha$.
*Ref:* [Poutanen & Beloborodov (2006)](references.md#poutanen--beloborodov-2006);
rotation terms [Bogdanov et al. (2019) §2](references.md#bogdanov-et-al-2019--paper-ii-l26).
*Check:* SD1a reference waveform, residual 0.11% with exact bending (`code_comparison.py`).

*where the spot–observer angle at each phase is pure spherical geometry:*

**(5) Rotation geometry** — `cos_psi()` · `src/mcrt/pulse.py`

$$\cos\psi(\varphi) = \cos i \cos\theta_s + \sin i \sin\theta_s \cos\varphi$$

As the star turns, the spot swings between its closest approach
$\cos(i-\theta_s)$ at $\varphi = 0$ and farthest $\cos(i+\theta_s)$ at $\varphi = \pi$.
*Ref:* standard; [Poutanen & Beloborodov (2006)](references.md#poutanen--beloborodov-2006).

*where gravity maps that geometric angle to the emission angle — two interchangeable maps:*

**(6) Light bending, linear (default)** — `bend()`, `visibility_threshold()` · `src/mcrt/pulse.py`

$$\cos\alpha = u + (1-u)\cos\psi, \qquad D = \frac{d\cos\alpha}{d\cos\psi} = 1-u, \qquad \text{visible} \iff \cos\psi \ge \frac{-u}{1-u}$$

Beloborodov's one-line approximation (accurate to ~1% for $u \lesssim 0.5$); the
visibility threshold is *negative*, so you literally see part of the star's far side.
*Ref:* [Beloborodov (2002)](references.md#beloborodov-2002).
*Check:* `tests/test_bending.py`; $u \to 0$ gives $\cos\alpha = \cos\psi$ exactly.

**(7) Light bending, exact** — `deflection_angle()`, `ExactBending` · `src/mcrt/bending.py`

$$\psi(\alpha) = \int_0^1 \frac{\tilde\beta\; dx}{\sqrt{1 - \tilde\beta^2 x^2 (1 - u x)}}, \qquad \tilde\beta = \frac{\sin\alpha}{\sqrt{1-u}}, \qquad D = \frac{d\cos\alpha}{d\cos\psi}\ \text{(numerical)}$$

The true general-relativistic ray path — Schwarzschild means the exact spacetime
around a non-spinning mass — integrated by Gauss–Legendre quadrature (a standard
numerical-integration rule), tabulated once per $u$, then inverted. Valid below
the photon sphere, $u < 2/3$ (the compactness beyond which grazing light orbits
the star and the map breaks down); at $u = 0$ it
returns $\psi = \alpha$ exactly. Used at the anchors since Gate G1: at J0740's
$u = 0.494$ the linear map biased $\Delta\mathrm{PF}$ by $-0.0066$ and hid a
~14% per-spot eclipse.
*Ref:* [Pechenick, Ftaclas & Cohen (1983)](references.md#pechenick-ftaclas--cohen-1983);
context [La Placa et al. (2019)](references.md#la-placa-et-al-2019).
*Check:* flat-space limit test; SD1a residual 0.80% → 0.11% (`b2_exact_bending.py`).

---

## Part III — the rotation terms (off by default)

Passing `rotation=Rotation(ν, R)` — spin frequency and stellar radius — switches
on the special-relativistic layer;
$\nu \to 0$ recovers the frozen chain bit-for-bit. **Aberration (12) is the term
that couples to the beaming swap** — it moves the $\mu$ at which $I(\mu)$ is
sampled, so rotation and the swap don't factor apart.
*Ref for (8)–(14):* [Bogdanov et al. (2019) §2](references.md#bogdanov-et-al-2019--paper-ii-l26)
(their eq. numbers in brackets). *Check:* SD1c waveform reproduction, ~1% (`c1_doppler_validate.py`).

**(8) Spot speed** [eq. 11] — `spot_speed()` · `src/mcrt/rotating.py`

$$\beta = \frac{2\pi\nu R \sin\theta_s}{c\,\sqrt{1-u}}$$

The $\sqrt{1-u}$ is the gravitational-redshift enhancement of the locally
measured velocity; J0740 at 346.5 Hz gives $\beta \approx 0.12$.

**(9) Lorentz factor** — `lorentz_gamma()` · `src/mcrt/rotating.py`

$$\gamma = \frac{1}{\sqrt{1-\beta^2}}$$

Enters the flux once, as the comoving-area factor in (4).

**(10) Ray–velocity angle** — `cos_xi()` · `src/mcrt/rotating.py`

$$\cos\xi = \frac{-\sin\alpha\, \sin i\, \sin\varphi}{\sin\psi}$$

The projection of the photon direction onto the spot's (azimuthal) velocity; at
$\varphi = 0, \pi$ the motion is transverse and $\cos\xi = 0$.

**(11) Doppler factor** [eq. 12] — `doppler_factor()` · `src/mcrt/rotating.py`

$$\delta = \frac{1}{\gamma\,(1 - \beta\cos\xi)}$$

$\delta > 1$ approaching (brighter, bluer), $\delta < 1$ receding; transverse
phases sit at $1/\gamma$.

**(12) Aberration** [eq. 13] — `point_spot_flux()` · `src/mcrt/pulse.py`

$$\mu' = \delta\cos\alpha \quad (\text{clipped to } [-1, 1])$$

The same motion that boosts the flux also tilts the apparent emission angle —
both the beaming lookup $I(\mu')$ and the projection cosine in (4) use $\mu'$,
which is consistent because a surface element's *projected* area comes out the
same in either frame ($dS\cos\alpha = dS'\cos\alpha'$, with $dS$ the emitting area).

**(13) Bolometric boost** [eq. 20] — `Rotation.flux_exponent` · `src/mcrt/rotating.py`

$$B = \delta^{\,n}, \qquad n = 4\ (\text{energy flux})\ \text{or}\ 3\ (\text{photon flux})$$

The grey library $I(\mu)$ is an energy intensity, so $\delta^4$ is the consistent
bolometric default.

**(14) Band-limited boost** — `band_boost()` · `src/mcrt/rotating.py`

$$B = \frac{\Phi_k(\delta)}{\Phi_k(1)}, \qquad \Phi_k(\delta) = \int_{E_1}^{E_2} \frac{E^k\, dE}{\exp\!\big[E / \big(\delta\sqrt{1-u}\; kT\big)\big] - 1}, \qquad k = 2\ (\text{counts})\ \text{or}\ 3\ (\text{energy})$$

Here $E$ is photon energy, $[E_1, E_2]$ the instrument band (NICER's calibrated
0.3–1.5 keV for J0740), and $kT$ the spot's temperature in energy units. For a
blackbody (an ideal thermal emitter — the stand-in for the fitted atmosphere
spectra) seen through a fixed band, every Doppler and redshift factor collapses
into one effective temperature $\delta\sqrt{1-u}\,kT$; a very wide band recovers
$\delta^3/\delta^4$, but when the band sits on the exponentially falling
high-energy (Wien) tail of the spectrum — the J0740 regime — the weight is
*steeper than any power* of $\delta$.
*Check:* the $\delta^3/\delta^4$/band spread is the $+0.02\ldots+0.06$ convention range (`c3_band_doppler.py`).

*Out of production scope, both beaming-independent: light-travel delay
(`travel_time_delay`, Bogdanov's eq. 18 — driver-level, used only by `c1_doppler_validate.py`
and `c4_caveat_audit.py`) and oblateness (the spinning star's slight equatorial
bulge, ~3% at 346.5 Hz,
[AlGendy & Morsink 2014](references.md#algendy--morsink-2014)).*

---

## Part IV — where $I(\mu;\tau)$ comes from (the Monte Carlo)

This is the only stochastic part, run once by `scripts/tau_sweep.py` for
$\tau \in \{0.1, 0.3, 1, 3, 10, 30\}$ × 5 seeds (~$4 \times 10^5$ *escaped*
photons each) and saved to `data/beaming_library.npz`; everything above reads
the file. Downstream uses the $\tau = 10$ row — $\tau = 10$ and $30$ agree
within error (the thick-slab plateau).

**(15) The lookup** — `beaming_lookup()` · `src/mcrt/beaming.py`

$$I(\mu) = \text{linear interpolation of } \big(\bar\mu_j,\, I_j\big), \quad \text{held flat outside the tabulated range}$$

Turns one saved library row into the plain callable that (4) samples; the flat
clamp is deliberate — the grazing $\mu \to 0$ tail is the noisiest data, and
extrapolating it would amplify Monte Carlo noise.
*Check:* tail not load-bearing — perturbing/splicing it moves $\Delta\mathrm{PF}$ by $\le 0.006$ (`e2_tail_sensitivity.py`).

*where each tabulated intensity comes from binned escape counts:*

**(16) Counts → intensity** — `curve_from_counts()` · `scripts/tau_sweep.py`

$$I_j = \frac{N_j / \bar\mu_j}{a}$$

$N_j$ is the number of photons that escaped into $\mu$-bin $j$ (20 bins across
$0 \le \mu \le 1$). A photon escaping at angle $\theta$ contributes
$\mu = \cos\theta$ to the outward flux, so binned counts measure flux
$\propto I(\mu)\,\mu$ — dividing by the bin's **mean escaped $\mu$** (fix 8b; the
geometric bin center would bias the grazing tail) recovers intensity, and
dividing by $a$ — the overall brightness scale fitted in (18), averaged over ~18
bins (fix 8a) — normalizes the curve without inheriting one noisy edge bin's error.
*Note:* `beaming.extract_intensity` is the simpler bin-center version kept for
unit tests; this pooled path is what builds the production library.
*Check:* validation against the exact $H(\mu)$, compared in flux space ($I\mu$ — the quantity the counts directly measure): reduced $\chi^2 = 0.70$ (`d_isotropic_validate.py`).

**(17) The shared $\mu$ grid** — `pooled_mu_grid()` · `scripts/tau_sweep.py`

$$\bar\mu_j = \frac{\sum_{\text{runs}} \sum_{\text{escapes in bin } j} \mu}{\sum_{\text{runs}} N_j}$$

Escapes crowd the high edge of each bin, so the pooled (over all runs) mean is
the low-variance estimate of the $\mu$ each count actually represents.

**(18) The limb-darkening fit** — `fit_limb_darkening()` · `src/mcrt/beaming.py`

$$I(\mu) = a\,(1 + b\mu), \qquad \text{least squares over bins with } \mu > 0.1$$

$a$ is the overall brightness scale (the intercept, reused as the normalizer in
(16)); $b$ is the one-number summary of the limb darkening — how steeply the
surface dims from face-on toward grazing: thin slab
$b \to 0$ (nearly isotropic), thick slab $b \approx 1.79$, bracketed by
[Eddington's $1 + 1.5\mu$](references.md#eddington-limb-darkening) below and
[Chandrasekhar's $H(\mu)$](references.md#chandrasekhar-1960) above.
*Check:* the $\Delta\mathrm{PF}$ magnitude tracks $b$ across all three laws — the effect is a property of limb darkening, not of our slab (`e3_atmosphere_laws.py`).

*where the counts $N_j$ are the histogram of escape angles from the photon random walk
(`Simulation.run()` · `src/mcrt/monte_carlo.py`), whose sampling rules are:*

**(19) Injection** — `Simulation.run()` · `src/mcrt/monte_carlo.py`

$$\cos\theta_0 = \sqrt{U}$$

A thermal source is equally bright in every direction, so the photon *count*
crossing the base goes as $N(\mu) \propto \mu$; sampling $\sqrt{U}$ reproduces
exactly that (Lambertian) law, where uniform sampling would over-produce grazing
photons. The base is a perfect absorber, so a photon that random-walks back down
is lost: $I(\mu;\tau)$ is the *transmission* beaming of a passive scattering
layer, conditioned on escape.
*Ref:* deep-dive v0.6.1; method [Cashwell & Everett (1959)](references.md#cashwell--everett-1959), [Whitney (2011)](references.md#whitney-2011).

**(20) Step length** — `sample_step_size()` · `src/mcrt/utils.py`

$$\Delta\tau = -\ln U$$

Free paths between scatterings are exponentially distributed,
$P(\Delta\tau) = e^{-\Delta\tau}$ — one mean free path on average, which is the
definition of optical depth.
*Check:* mean free path ≈ 1 (`validate_engine.py`); sampler distributions (`tests/test_physics.py`).

**(21) Scattering angle** — `sample_thomson_angle()` · `src/mcrt/utils.py`

$$P(\mu_{\text{sc}}) = \tfrac{3}{4}\big(1 + \mu_{\text{sc}}^2\big)$$

$\mu_{\text{sc}}$ is the cosine of the angle between the photon's directions
before and after the scatter, and $P$ is the **phase function** — the probability
distribution of scattering directions — here Thomson's, the physical law for a
photon bouncing off a free electron. It is sampled by rejection (draw a
candidate, keep it with probability $\propto P$), with a uniform azimuth around
the old direction (`rotate_vector`). The opt-in
`phase_function="isotropic"` ($P = \tfrac12$, `sample_isotropic_angle`) exists
solely so the machinery can be checked against Chandrasekhar's exact $H(\mu)$
solution — see [Two meanings of "isotropic"](#two-meanings-of-isotropic).
*Check:* Thomson vs. isotropic emergent $I(\mu)$ differ ≤ 2.8% — the measured phase-function sensitivity.

**(22) Boundaries** — `Simulation.run()` · `src/mcrt/monte_carlo.py`

$$\tau_{\text{pos}} \le 0 \;\Rightarrow\; \text{escape, record } \mu = -d_z; \qquad \tau_{\text{pos}} \ge \tau_{\text{total}} \;\Rightarrow\; \text{absorbed (lost)}$$

$\tau_{\text{pos}}$ is the photon's current vertical depth in the slab (it is
injected at the base, $\tau_{\text{pos}} = \tau_{\text{total}}$, the slab's full
optical thickness; the surface is $\tau_{\text{pos}} = 0$), and $d_z$ is the
vertical component of its direction vector. The recorded escape cosines are the
entire hand-off from the Monte Carlo to everything above: histogram them (16),
fit them (18), save them — and the chain never touches randomness again.
*Check:* energy conservation (`validate_engine.py`); error $\propto N^{-1/2}$ with no systematic floor (`convergence.py`).

---

## The pen-and-paper cross-check

**(23) Closed-form PF** — `analytic_isotropic_pf()` · `src/mcrt/pulse.py`

$$\mathrm{PF} = \frac{(1-u)\,(\cos\psi_{\max} - \cos\psi_{\min})}{2u + (1-u)\,(\cos\psi_{\max} + \cos\psi_{\min})}, \qquad \cos\psi_{\max} = \cos(i-\theta_s),\quad \cos\psi_{\min} = \cos(i+\theta_s)$$

Collapse the whole chain by hand for an always-visible isotropic spot on a
frozen star — substituting (6) into (4) with $I \equiv 1$ makes $F$ monotonic in
$\cos\psi$, so the extremes sit at $\varphi = 0, \pi$ and (2) evaluates in closed
form; this is the benchmark the numerical chain must reproduce (it raises if the
spot ever sets, where the assumption breaks).
*Check:* `tests/test_pulse.py`.

---

## The experiment scripts — which equations each one exercises

Every script runs the same chain; the columns say what it swaps and what it
reports. **(V)** marks a validation script (checks the machinery against an
external truth), **(E)** a robustness experiment (checks the result doesn't
hinge on a modeling choice).

| script | chain used | what is swapped | result |
|---|---|---|---|
| `beaming_pulse_sweep.py` | (2)–(6), (15)–(18) | $I \equiv 1$ vs. each $\tau$ row | $\Delta\mathrm{PF}(\tau)$ for 3 always-visible geometries |
| `j0030_anchor.py`, `j0740_anchor.py` | (1)–(7), (15) | beaming, per published fit; exact bending | J0030: spots fail to tile → PF pins at 1, $\Delta\mathrm{PF} \approx 0$, systematic moves into waveform shape; J0740: $+0.137$/$+0.195$ static |
| `phase_diagram.py` | (1)–(6) on a grid | two-spot geometry $(\Delta\varphi, \theta)$ | map of live vs. saturated $\Delta\mathrm{PF}$; analytic boundary from single-spot eclipsed fraction |
| `c2_doppler_coupling.py` | + (8)–(13) | rotation off/on × beaming | static $+0.137/+0.195$ → spinning $+0.037/+0.061$ ($\delta^4$); the peak-normalized shape metric barely moves — the systematic is routed, not erased |
| `c3_band_doppler.py` | + (14) | $\delta^3$ / $\delta^4$ / NICER band | diluted headline is convention-dependent: $+0.02\ldots+0.06$ |
| `c4_caveat_audit.py` | (1) + delay/oblateness bounds | — | caveat numbers re-sized to the diluted claim |
| `code_comparison.py` (V) | (2)–(7) | linear vs. exact bending at SD1a | residual vs. independent reference waveform: 0.80% → 0.11% |
| `c1_doppler_validate.py` (V) | + (8)–(13) + delay | SD1c − SD1a | isolates pure Doppler + aberration; ~1% (the linear-bending floor) |
| `b2_exact_bending.py` (V) | (6) vs. (7) | bending map only | $\Delta\mathrm{PF}$ shift $-0.0066$ at $u = 0.494$ → Gate G1 enacted |
| `d_isotropic_validate.py` (V) | (16)–(22), isotropic mode | phase function | flux-space $\chi^2 = 0.70$ vs. exact $H(\mu)$; Thomson control ≤ 2.8% |
| `e2_tail_sensitivity.py` (E) | (15)–(16) | grazing-bin perturb/splice | $\lvert\delta(\Delta\mathrm{PF})\rvert \le 0.006$ — tail not load-bearing |
| `e3_atmosphere_laws.py` (E) | (15) | Eddington ($b{=}1.50$) / $H(\mu)$ ($b{=}1.68$) / slab ($b{=}1.79$) | same sign & geometry-routing; magnitude tracks $b$ |
| `finite_cap.py` (E) | (3)–(4), ring-tiled caps | point spot vs. finite cap | $\Delta\mathrm{PF}$ moves ~$-0.003$ — a measured small bias |
| `validate_phase_diagram.py` (V) | (3), (5)–(6) | — | engine eclipse onset ≡ analytic condition; PB06 antipodal classes exact |

---

## Quantities that live between equations

- **The tiling criterion** — PF saturates at 1 exactly when the *summed* flux
  (3) touches zero, i.e. when the spots between them fail to keep something
  visible at every phase; it is *not* single-spot eclipse (J0740-Miller's spots
  each hide ~21% of the cycle, yet the anti-phased pair tiles the rotation and
  PF stays live). Novelty claim, pending advisor verification.
- **Fourier harmonics ($A_0, A_1, A_2, \ldots$)** — the pulse $F(\varphi)$ is
  periodic in phase, so it decomposes into a Fourier series; $A_0$ is the mean
  level, $A_1$ the fundamental (once-per-rotation) amplitude, $A_2$ the second
  harmonic (twice-per-rotation), and so on. These are standard X-ray-pulsar
  observables. The ratio $A_1/A_0$ is an alternative pulsed-amplitude measure to
  PF (eq 2) — the two agree only for a pure sinusoid. This matters because spin
  (Part III) pumps a large *common* $A_2$ into both the isotropic and realistic
  profiles, saturating the min/max PF while leaving $A_1/A_0$ and the waveform
  shape largely intact (see [`claims-evidence.md`](claims-evidence.md) A4).
- **Error bars** — Monte Carlo quantities carry seed-to-seed standard deviations
  (5 seeds): $\pm 0.003$ (Riley) / $\pm 0.005$ (Miller) on the headline.
  Deterministic quantities are exact to grid resolution.

### Two meanings of "isotropic"

- **Isotropic *scattering*** — a phase-function choice inside the Monte Carlo,
  eq. (21): each scatter goes in a fully random direction. Used only by the
  transport validation `d_isotropic_validate.py`.
- **Isotropic *beaming*** — the finished surface glowing equally in all
  directions, $I \equiv 1$ in eq. (4). This is the baseline of the paper's comparison.

The headline experiment is **not** "isotropic vs. Thomson scattering" — it is one
Thomson-built brightness curve $I(\mu;\tau)$ versus a flat $I \equiv 1$, pushed
through identical geometry.

---

## What validates what

| Equation(s) | Claim | Check |
|---|---|---|
| (19)–(21) | free-path & phase-function sampling correct | `tests/test_physics.py` |
| (20), (22) | energy conservation, mean free path ≈ 1 | `scripts/validate_engine.py` |
| (16)–(22) | transport machinery exact | isotropic mode vs. $H(\mu)$, flux-space $\chi^2 = 0.70$ (`tests/test_transport.py`, v0.9.11) |
| (15)–(18) | emergent $I(\mu)$ at $\tau = 10$ physical | Eddington/$H(\mu)$ bracket, ≤ 2.8% from $H$ (measured) |
| (2)–(6) | pulse machinery (geometry + bending + flux) | closed form (23); SD1a waveform (`tests/test_pulse.py`) |
| (7) | exact bending map | flat-space limit; SD1a residual 0.80% → 0.11% (`tests/test_bending.py`, v0.9.9) |
| (8)–(13) | Doppler + aberration layer | SD1c waveform, ~1% (`c1_doppler_validate.py`, v0.9.10) |
| (3), (6) | eclipse/visibility logic | analytic onset; PB06 antipodal classes (`test_pulse.py`) |
| (3)–(4) | point-spot reduction harmless | finite-cap tiling, ~$-0.003$ (`finite_cap.py`) |
| (15)–(16) | grazing tail not load-bearing | E2 perturb/splice, ≤ 0.006 (`tests/test_sensitivity.py`, v0.9.12) |
| (15), (18) | not-just-our-slab robustness | E3 Eddington/$H$ swap rows (v0.9.12) |
| (16)–(22) | photon-count budget sufficient | convergence study, error $\propto N^{-1/2}$ (v0.9.7) |
