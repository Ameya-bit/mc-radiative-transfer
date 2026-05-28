import numpy as np
import pytest

from mcrt import sample_step_size, sample_thomson_angle, get_random_direction, rotate_vector

def test_sample_step_size():
    """Verify that step sizes are positive and follow an exponential distribution roughly."""
    samples = [sample_step_size() for _ in range(10000)]
    assert all(s > 0 for s in samples)
    assert np.isclose(np.mean(samples), 1.0, rtol=0.1)

def test_sample_thomson_angle():
    """Verify Thomson angles are within [-1, 1]."""
    samples = [sample_thomson_angle() for _ in range(1000)]
    assert all(-1 <= s <= 1 for s in samples)

def test_get_random_direction():
    """Verify unit vector length."""
    for _ in range(100):
        v = get_random_direction()
        assert np.isclose(np.linalg.norm(v), 1.0)

def test_rotate_vector_conservation():
    """Verify that rotation preserves the norm of the vector."""
    v_orig = np.array([0, 0, 1])
    for _ in range(100):
        costheta = sample_thomson_angle()
        v_new = rotate_vector(v_orig, costheta)
        assert np.isclose(np.linalg.norm(v_new), 1.0)
