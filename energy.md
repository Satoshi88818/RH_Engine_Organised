# `rh_energy.py`: a candidate E(σ) computed from known zeta zeros

## What this is for

Some proposed approaches to the Riemann Hypothesis (e.g. the "Dynamic Spectral Magnetic Field" repo) describe an "energy" E(σ) that is claimed to be strictly convex with a unique minimum at σ = ½. This script makes that idea concrete and testable by defining one precise E(σ) and computing it numerically from the first N nontrivial zeros.

It is a **sanity-check and teaching tool, not a proof or evidence for RH.**

## The definition

For s = σ + it, write ξ(s) for the completed zeta function. Over a window t ∈ [T₀, T₁]:

> **E(σ) = mean over t of log |ξ(σ + it) / ξ(0)|**

Computed from zeros via the Hadamard product. Pairing ρ = ½ + iγ with its conjugate gives, with x = σ − ½:

> ξ(s)/ξ(0) = ∏ₙ [ ((x + it)² + γₙ²) / (¼ + γₙ²) ]

This product converges absolutely. The zeros are fetched with `mpmath.zetazero` and cached in a `.npy` file (~50 s for 300 zeros the first time).

## Usage

```bash
pip install mpmath numpy matplotlib
python rh_energy.py 300     # N zeros; saves rh_energy.png
```

## What it shows (300 zeros, t ∈ [30, 271])

- E(σ) is exactly symmetric about σ = ½ (error 0.00).
- Its minimum is at σ = 0.50.
- E″(σ) ≈ 0 away from ½ (mean |E″| ≈ 0.008) and has a sharp spike at ½.

So E looks like a **V shape (a cusp), not a smooth strictly convex bowl.**

## How to interpret this honestly

1. **The minimum at ½ is an input, not an output.** The script feeds in zeros that are *already on the line* (RH is numerically verified to enormous height, so `zetazero` returns line zeros). The V-shape just reflects that. If you could feed in a hypothetical off-line zero, E would develop its kink at that zero's real part instead.
2. **This is a known fact in disguise (Littlewood's lemma / Jensen's formula).** The second derivative of the averaged log|ξ| is, asymptotically, a measure supported on the real parts of the zeros. Hence RH ⇔ E is piecewise linear with a single kink at ½. This is a correct *reformulation* of RH, not a proof.
3. **"Strictly convex everywhere" is the wrong target.** Under RH, E has zero curvature for σ ≠ ½. And convexity alone is automatic: for any analytic function, such averages are convex in σ (Hardy's convexity theorem), so convexity cannot discriminate between zero locations.
4. **Where a real proof would have to work:** showing the curvature measure has no mass off σ = ½ without assuming the zero locations. That is the open problem.

## Limitations

- Finite window and finite N introduce edge effects; keep T₁ well below the largest γ used (the script uses T₁ = γ_N / 2).
- Numerical results can never establish RH.
- Not independent of the zero data; to test the framework you'd need a definition of E that doesn't take zero locations as input.

## Suggested extensions

- Compare against direct evaluation of ξ(σ + it) with `mpmath` to check the product truncation.
- Insert a *fake* zero pair at ½ ± δ + iγ₀ and watch E″ develop mass at σ = ½ ± δ.
- Vary N and the window to see the cusp slope track ½·log(t/2π).

---

## Script

```python
"""Candidate energy E(sigma) built from known zeta zeros (assumed on the line)."""
import os, sys
import numpy as np

N_ZEROS = int(sys.argv[1]) if len(sys.argv) > 1 else 300
T0 = 30.0                                  # start of the t-window
CACHE = f"zeros_{N_ZEROS}.npy"

def load_zeros(n):
    if os.path.exists(CACHE):
        return np.load(CACHE)
    import mpmath as mp                      # slow: ~0.2 s per zero
    g = np.array([float(mp.zetazero(k).imag) for k in range(1, n + 1)])
    np.save(CACHE, g)
    return g

def energy(sigma, gammas, t):
    """E(sigma) = mean over t in window of log|xi(sigma+it)/xi(0)|,
    using the Hadamard product over zeros rho = 1/2 +- i*gamma:
    pair factor = ((x+it)^2 + gamma^2) / (1/4 + gamma^2), x = sigma - 1/2."""
    x = sigma - 0.5
    z = (x + 1j * t)[:, None] ** 2 + gammas[None, :] ** 2
    return np.log(np.abs(z) / (0.25 + gammas ** 2)[None, :]).sum(axis=1).mean()

if __name__ == "__main__":
    gam = load_zeros(N_ZEROS)
    T1 = 0.5 * gam[-1]                        # stay well below the truncation edge
    t = np.linspace(T0, T1, 4000)
    sig = np.linspace(0.0, 1.0, 101)
    E = np.array([energy(s, gam, t) for s in sig])
    h = sig[1] - sig[0]
    Epp = np.gradient(np.gradient(E, h), h)   # numerical E''(sigma)
    print(f"{len(gam)} zeros, window t in [{T0}, {T1:.0f}]")
    print(f"min of E at sigma = {sig[E.argmin()]:.2f}; E(1/2)={E[50]:.3f}")
    print(f"E(0)-E(1/2) = {E[0]-E[50]:.3f}   symmetry error = {np.abs(E-E[::-1]).max():.2e}")
    away = np.abs(sig - 0.5) > 0.15
    print(f"mean |E''| away from 1/2: {np.abs(Epp[away]).mean():.3f}; peak at 1/2: {Epp[50]:.1f}")
    try:
        import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
        f, (a, b) = plt.subplots(1, 2, figsize=(10, 3.5))
        a.plot(sig, E); a.set(title="E(sigma)", xlabel="sigma")
        b.plot(sig[2:-2], Epp[2:-2]); b.set(title="E''(sigma)", xlabel="sigma")
        plt.tight_layout(); plt.savefig("rh_energy.png", dpi=130)
    except ImportError:
        pass
```
