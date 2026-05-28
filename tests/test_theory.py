import numpy as np

from mcrt import chandrasekhar_h, eddington_limb_darkening


def test_eddington_endpoints():
    """Eddington law is 1 at μ=0 and 2.5 at μ=1."""
    assert np.isclose(eddington_limb_darkening(0.0), 1.0)
    assert np.isclose(eddington_limb_darkening(1.0), 2.5)


def test_h_function_at_zero_is_one():
    """H(0) = 1 exactly for any albedo (definitional)."""
    assert np.isclose(chandrasekhar_h(0.0, albedo=1.0)[0], 1.0)
    assert np.isclose(chandrasekhar_h(0.0, albedo=0.5)[0], 1.0)


def test_h_function_monotonic_increasing():
    """H(μ) increases monotonically across [0, 1]."""
    mu = np.linspace(0, 1, 25)
    H = chandrasekhar_h(mu, albedo=1.0)
    assert np.all(np.diff(H) > 0)


def test_h_function_conservative_endpoint():
    """Conservative isotropic H(1) ≈ 2.908 (Chandrasekhar's tabulated value)."""
    H1 = chandrasekhar_h(1.0, albedo=1.0)[0]
    assert 2.85 < H1 < 2.95


def test_h_function_grows_with_albedo():
    """A more reflective atmosphere (higher ω) yields a larger H."""
    h_low = chandrasekhar_h(1.0, albedo=0.5)[0]
    h_high = chandrasekhar_h(1.0, albedo=1.0)[0]
    assert h_high > h_low > 1.0
