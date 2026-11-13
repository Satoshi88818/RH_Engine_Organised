"""Basic smoke tests for the energy module."""
import numpy as np
import pytest

from rh_engine.energy import energy, compute_energy_profile


def test_energy_symmetry():
    """E(σ) should be symmetric about 1/2 when zeros are on the line."""
    # Use a tiny set of artificial zeros for speed
    gammas = np.array([14.1347, 21.0220, 25.0109])
    t = np.linspace(30, 50, 200)
    e0 = energy(0.0, gammas, t)
    e1 = energy(1.0, gammas, t)
    e05 = energy(0.5, gammas, t)
    assert abs(e0 - e1) < 1e-10
    assert e05 < e0


def test_compute_profile_runs():
    """Smoke test that the full profile computation does not crash on tiny N."""
    # This will try to load/compute zeros; for unit tests we keep N small
    # and rely on the cache or a short run.
    sig, E, Epp = compute_energy_profile(n_zeros=5, n_t=100, n_sigma=21)
    assert len(sig) == 21
    assert len(E) == 21
    assert np.isfinite(E).all()
