"""Track E2/E3 — tail-sensitivity and atmosphere-law robustness driver helpers.

E2 asks whether the noisy, clamped grazing tail of the beaming library destabilizes
ΔPF; E3 asks whether the systematic survives swapping the Thomson slab for independent
analytic limb-darkening laws. These tests certify the driver helpers behave — the
perturbation is a no-op at zero σ, the H-splice touches only the tail, ΔPF is invariant
to a law's overall normalization, and analytic laws reproduce the positive, spin-diluted
signal — plus a data-gated check that the drivers reproduce the published headline.

The driver modules live in ``scripts/`` (alongside ``anchor_lib``), so the scripts
directory is put on ``sys.path`` here the way the shell invocation's PYTHONPATH does.
"""

import sys
from pathlib import Path

import numpy as np
import pytest

SCRIPTS = Path(__file__).resolve().parents[1] / "scripts"
if str(SCRIPTS) not in sys.path:
    sys.path.insert(0, str(SCRIPTS))

from mcrt import ExactBending, Rotation, chandrasekhar_h, eddington_limb_darkening

import e2_tail_sensitivity as e2
import e3_atmosphere_laws as e3
import j0740_anchor as j0740

LIBRARY = Path(__file__).resolve().parents[1] / "data" / "beaming_library.npz"

# A synthetic, monotone I(μ) on a 20-bin μ-grid — enough to exercise the helpers
# without depending on the production library.
MU = np.linspace(0.03, 0.98, 20)
CURVE = 1.0 + 1.5 * MU
SIGMA = 0.05 * np.ones_like(CURVE)
TAIL = MU < e2.MU_TAIL


# --- E2 helper unit tests ---------------------------------------------------

def test_zero_perturbation_reproduces_baseline_bitwise():
    """σ = 0 ⇒ every resampled draw equals the unperturbed ΔPF and the spread is 0."""
    anchor = j0740.RILEY
    dpf = e2.make_delta_pf(anchor, MU, ExactBending(anchor.compactness), None)
    rng = np.random.default_rng(0)
    baseline, sigma, draws = e2.resample_sigma(
        dpf, CURVE, np.zeros_like(SIGMA), TAIL, rng, n_draws=64)
    assert sigma == 0.0
    assert np.all(draws == baseline)


def test_h_splice_touches_only_the_tail():
    """The H-tail splice modifies μ < 0.1 bins and leaves the trusted curve untouched."""
    spliced = e2.h_spliced_curve(CURVE, MU, TAIL)
    assert np.array_equal(spliced[~TAIL], CURVE[~TAIL])       # trusted part preserved
    assert np.all(spliced[TAIL] != CURVE[TAIL])               # tail replaced
    # Chandrasekhar H is monotone increasing, so the spliced tail rises with μ and
    # stays below the first trusted bin it was anchored to (a continuous join).
    assert np.all(np.diff(spliced[TAIL]) > 0)
    first_kept = int(np.argmax(~TAIL))
    assert np.all(spliced[TAIL] < CURVE[first_kept])


def test_resample_sigma_is_small_and_positive_with_real_scatter():
    """A nonzero tail σ induces a nonzero but modest ΔPF spread (no blow-up)."""
    anchor = j0740.RILEY
    dpf = e2.make_delta_pf(anchor, MU, ExactBending(anchor.compactness), None)
    rng = np.random.default_rng(1)
    _base, sigma, _draws = e2.resample_sigma(dpf, CURVE, SIGMA, TAIL, rng, n_draws=200)
    assert 0.0 < sigma < 0.05


# --- E3 helper unit tests ---------------------------------------------------

def test_delta_pf_invariant_to_law_normalization():
    """ΔPF depends only on I(μ)'s shape, not its overall scale (PF is a ratio)."""
    anchor = j0740.RILEY
    bending = ExactBending(anchor.compactness)
    dpf_unit, _ = e3.anchor_metrics(anchor, eddington_limb_darkening, bending, None)
    scaled = lambda mu: 3.7 * eddington_limb_darkening(mu)  # noqa: E731
    dpf_scaled, _ = e3.anchor_metrics(anchor, scaled, bending, None)
    assert dpf_unit == pytest.approx(dpf_scaled, abs=1e-12)


@pytest.mark.parametrize("law", [eddington_limb_darkening, chandrasekhar_h])
def test_analytic_law_is_pf_live_and_spin_dilutes(law):
    """Each independent law gives a positive static ΔPF that the real spin dilutes."""
    anchor = j0740.RILEY
    bending = ExactBending(anchor.compactness)
    rot = Rotation(spin_hz=e3.J0740_SPIN_HZ, radius_km=e3.RADIUS_KM[anchor.label])
    dpf0, rms0 = e3.anchor_metrics(anchor, law, bending, None)
    dpf_rot, _ = e3.anchor_metrics(anchor, law, bending, rot)
    assert dpf0 > 0.10                    # PF-live under an independent law
    assert 0.0 < dpf_rot < dpf0           # spin dilutes but does not flip the sign
    assert rms0 > 0.05                    # a real waveform-shape difference exists


# --- Data-gated cross-checks against the published headline ------------------

@pytest.mark.skipif(not LIBRARY.exists(),
                    reason="production beaming library not present on a fresh checkout")
def test_e2_baseline_reproduces_published_headline():
    """E2's clamp baseline equals the exact-bending J0740 headline (Riley/Miller)."""
    d = np.load(LIBRARY)
    from anchor_lib import shape_tau_index
    ti10 = shape_tau_index(d["tau_values"])
    mu = d["mu_centers"]
    curve = d["intensity_by_tau"][ti10]
    expected = {"Riley 2021 (X-PSI ST-U)": 0.137,
                "Miller 2021 (Illinois–Maryland)": 0.195}
    for anchor in j0740.ANCHORS:
        dpf = e2.make_delta_pf(anchor, mu, ExactBending(anchor.compactness), None)
        assert dpf(curve) == pytest.approx(expected[anchor.label], abs=0.003)


@pytest.mark.skipif(not LIBRARY.exists(),
                    reason="production beaming library not present on a fresh checkout")
def test_e2_tail_effects_below_seed_error_bar():
    """The tail's induced σ and H-splice systematic are both under the ±0.003 seed bar."""
    d = np.load(LIBRARY)
    from anchor_lib import shape_tau_index
    ti10 = shape_tau_index(d["tau_values"])
    mu = d["mu_centers"]
    curve = d["intensity_by_tau"][ti10]
    sigma = d["intensity_std_by_tau"][ti10]
    tail = mu < e2.MU_TAIL
    anchor = j0740.RILEY
    dpf = e2.make_delta_pf(anchor, mu, ExactBending(anchor.compactness), None)
    rng = np.random.default_rng(7)
    base, sig, _ = e2.resample_sigma(dpf, curve, sigma, tail, rng, n_draws=300)
    spliced = e2.h_spliced_curve(curve, mu, tail)
    assert sig < 0.003
    assert abs(dpf(spliced) - base) < 0.003
