"""Candidate energy E(σ) computed from known zeta zeros.

E(σ) = mean over t of log |ξ(σ + it) / ξ(0)| via the Hadamard product.
This is a sanity-check / teaching tool, not evidence for RH.
"""

from __future__ import annotations

import os
from pathlib import Path

import numpy as np

def load_zeros(n: int, cache_dir: Path | str | None = None) -> np.ndarray:
    """Load or compute the first n nontrivial zeros (imaginary parts)."""
    cache_dir = Path(cache_dir or ".")
    cache = cache_dir / f"zeros_{n}.npy"
    if cache.exists():
        return np.load(cache)
    import mpmath as mp

    g = np.array([float(mp.zetazero(k).imag) for k in range(1, n + 1)])
    cache_dir.mkdir(parents=True, exist_ok=True)
    np.save(cache, g)
    return g


def energy(sigma: float, gammas: np.ndarray, t: np.ndarray) -> float:
    """Mean over t of log|ξ(σ+it)/ξ(0)| using the Hadamard product.

    Pairing ρ = ½ ± iγ gives the factor
        ((x + it)² + γ²) / (¼ + γ²),  x = σ − ½.
    """
    x = sigma - 0.5
    z = (x + 1j * t)[:, None] ** 2 + gammas[None, :] ** 2
    return float(np.log(np.abs(z) / (0.25 + gammas**2)[None, :]).sum(axis=1).mean())


def compute_energy_profile(
    n_zeros: int = 300,
    t0: float = 30.0,
    n_sigma: int = 101,
    n_t: int = 4000,
    cache_dir: Path | str | None = None,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return (sigma, E, E'') on a grid."""
    gam = load_zeros(n_zeros, cache_dir)
    t1 = 0.5 * gam[-1]
    t = np.linspace(t0, t1, n_t)
    sig = np.linspace(0.0, 1.0, n_sigma)
    E = np.array([energy(s, gam, t) for s in sig])
    h = sig[1] - sig[0]
    Epp = np.gradient(np.gradient(E, h), h)
    return sig, E, Epp


def main(n_zeros: int = 300, out_png: str = "rh_energy.png") -> None:
    sig, E, Epp = compute_energy_profile(n_zeros)
    print(f"{n_zeros} zeros")
    print(f"min of E at sigma = {sig[E.argmin()]:.2f}; E(1/2)={E[50]:.3f}")
    print(f"E(0)-E(1/2) = {E[0]-E[50]:.3f}   symmetry error = {np.abs(E-E[::-1]).max():.2e}")
    away = np.abs(sig - 0.5) > 0.15
    print(f"mean |E''| away from 1/2: {np.abs(Epp[away]).mean():.3f}; peak at 1/2: {Epp[50]:.1f}")

    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 3.5))
        ax1.plot(sig, E)
        ax1.set(title="E(sigma)", xlabel="sigma")
        ax2.plot(sig[2:-2], Epp[2:-2])
        ax2.set(title="E''(sigma)", xlabel="sigma")
        plt.tight_layout()
        plt.savefig(out_png, dpi=130)
        print(f"saved {out_png}")
    except ImportError:
        pass


if __name__ == "__main__":
    import sys

    n = int(sys.argv[1]) if len(sys.argv) > 1 else 300
    main(n)
