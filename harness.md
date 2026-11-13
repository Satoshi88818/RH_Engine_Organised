# Empirical harness: local collision times of the zeros of H_t on windows up to γ ≈ 10⁶

## Purpose

Test the picture from `README_no_collision_bound.md` on real zeta zeros: the **gap-lifetime lemma** (a real pair stays apart backward in time for at least r²/8) and how the available lifetime shrinks with height.

This is **evidence, not proof.** It does not touch the hard direction (t ≥ 0, where an off-line zero would live); see Limitations.

## What I built, and one change from the plan

The plan was to evaluate the Polymath15 A + B approximation at small t. I did something different, for a good reason: for **t ≤ 0 the flow is exact and simple**. Multiplying by e^{−su²} in Fourier space means Gaussian smoothing of H₀ along the real axis:

> H_{−s}(z) = (k_s ∗ H₀)(z), k_s = N(0, 2s), with H₀(z) ∝ (¼ + t²)·|Γ(¼ + it/2)|·Z(t), t = z/2 (so z = 2γ).

So no approximation is needed. Z(t) is Hardy's function (`mpmath.siegelz`, 20 digits). The smoothing is done spectrally on a window, after multiplying H₀ by e^{πz/8} to cancel its decay (a tilt that only shifts the Gaussian's mean by 2·(π/8)·s).

**Per window** (length 40 in z, grid step 0.1, 4× spectral up-sampling for sign-change counting):
1. Find the real zeros at s = 0 (spline + root-finding).
2. Increase s and record the first s_c at which the number of sign changes drops, i.e. two real zeros merge. The local collision time is t = −s_c. Counting is done between midpoints of the outermost zeros so that edge drift is not mistaken for a collision.
3. Compare s_c with the lemma bound min r²/8 over all consecutive pairs in the window.

**Windows:** 6 per height at γ ≈ 10², 10³, 10⁴, 10⁵, 10⁶, each 20 wide in γ. That gives 30 windows and 655 zeros. Runtime is about 4 minutes on one core.

## Validation

- At γ ≈ 100 the 9 zeros found match the known zeros 2γₙ to 1.6 × 10⁻⁷.
- **Bug caught during development:** an early version flagged a "collision" at 0.57 × the lemma bound. It was a zero drifting across a fixed window edge (count dropped by 1, not 2). Using midpoint-bounded counting fixed it, and the same window then gave ratio 1.14.

## Results

| γ | windows | zeros | mean gap (z) | min gap | min s_c | median s_c | min r²/8 | s_c / (r²/8) | windows with no collision by s = 0.30 |
|---|---|---|---|---|---|---|---|---|---|
| 10² | 6 | 53 | 3.62 | 1.432 | 0.2821 | 0.2850 | 0.2562 | 1.10 | 4 |
| 10³ | 6 | 85 | 2.44 | 0.442 | 0.02485 | 0.1134 | 0.02444 | 1.02–1.14 | 0 |
| 10⁴ | 6 | 130 | 1.67 | 0.371 | 0.01766 | 0.0394 | 0.01724 | 1.02–1.15 | 0 |
| 10⁵ | 6 | 172 | 1.29 | 0.407 | 0.02259 | 0.0300 | 0.02068 | 1.08–1.15 | 0 |
| 10⁶ | 6 | 215 | 1.05 | 0.052 | 0.00033 | 0.0065 | 0.00034 | 0.98–1.15 | 0 |

![Lemma vs measurement, and collision time vs height](rh_harness.png)

**Findings**

1. **The lemma holds, and is nearly sharp.** In 26 windows with a collision, s_c / (r²/8) ranged from 0.984 to 1.151 (median 1.08). 25 of 26 are ≥ 1. The one at 0.984 involves a gap of 0.052, close to the position-resolution limit (root error of order 10⁻³ gives a few percent error in r²), so I count it as consistent with the bound.
2. **Spacing scaling confirmed.** mean gap × log(γ/2π) = 12.5 at every height (12.33–12.59), matching 4π. So the natural lifetime scale (4π/log)²/8 falls from 1.6 at γ = 100 to 0.14 at 10⁶.
3. **Medians fall, minima are erratic.** Median s_c drops from 0.285 to 0.0065 over four decades. Minima are dominated by single close pairs: at γ ≈ 10⁶ one pair with gap 0.052 gives s_c = 0.00033.
4. **The 1/log²x guess is not confirmed.** Median s_c × log²(γ/2π) = 3.4, 3.0, 2.1, 2.8, 0.93. With 6 windows per height, noise swamps any clean law. The 10⁶ value is lower, as expected if close pairs become more extreme as more zeros are sampled.
5. **Consequence for the conjectured t(x) in `README_no_collision_bound.md`:** any workable t(x) must lie below r_min(x)²/8, and r_min over a stretch of N zeros shrinks with N (random-matrix heuristic: roughly N^(−1/3) in normalised units). So t(x) probably has to decay **faster** than 1/log²x. This is a heuristic, not tested here.
6. **Connection to the literature:** a close pair at height γ gives a lower bound Λ > −s_c. This harness reproduces that mechanism at small scale (e.g. Λ > −0.00033 from 215 zeros near 10⁶). The published bounds (about −10⁻¹¹) use far higher, carefully chosen pairs.

## Limitations

- **Only t ≤ 0 is probed.** Any off-line zero for the true ζ would show up at t > 0 (Λ > 0), which this harness cannot see. For t > 0 the flow is not a convolution; one would need either the A + B approximation or the identity H_t(z) = E[H₀(z − i√(2t)·W)], W ~ N(0,1), which uses ξ at complex arguments plus an argument-principle zero count. That is the natural next build.
- Windows are narrow (40 in z) and few (6 per height): small samples, so statistics are noisy.
- "First collision" is resolved to about 10⁻⁵ in s, but detection can miss pairs closer than the up-sampled grid (0.025 in z).
- The lemma's rigorous proof assumes all zeros are real (finite-system case; infinite extension assumed). Zeros outside a window are not controlled.
- Nothing here bears on whether any zero is off the critical line.

## Reproduce

```bash
pip install mpmath numpy scipy matplotlib
python rh_harness.py 100,1000,10000,100000,1000000 6 results.json
```

`results_final.json` holds the per-window output, including the zero positions.

## Script

```python
"""Empirical harness: backward heat flow H_{-s} of the Riemann Xi function on windows of the critical line.

H_{-s}(z) = (k_s * H_0)(z), k_s = N(0, 2s)   [exact; smoothing, no A+B approximation needed for t<=0].
H_0(z) ∝ (1/4 + t^2) |Gamma(1/4 + i t/2)| Z(t), t = z/2, so z = 2*gamma.
In each window we find the first s_c at which two real zeros merge (sign-change count drops),
i.e. the local collision time t = -s_c, and compare with the Gap-lifetime Lemma (s_c >= min r^2/8).
"""
import sys, json, time
import numpy as np, mpmath as mp
from scipy.interpolate import CubicSpline
from scipy.optimize import brentq

mp.mp.dps = 20
A = np.pi / 8                      # tilt that cancels the e^{-pi z/8} decay of H_0
DZ, MARGIN, LW = 0.1, 3.6, 40.0    # grid step (z), margin, window length (z)
UP = 4                             # spectral up-sampling for sign-change counting

def g_values(zgrid):
    out = np.empty(len(zgrid))
    for i, z in enumerate(zgrid):
        t = mp.mpf(z) / 2
        lg = mp.re(mp.loggamma(mp.mpc(0.25, t / 2)))
        out[i] = float(mp.siegelz(t) * (mp.mpf(1) / 4 + t * t) * mp.exp(lg + A * z))
    return out

def smooth(g, s, dz, up=1):
    """F_s = (N(2*A*s, 2s) * g): same zeros as H_{-s} (positive tilt removed). Spectral, optional up-sampling."""
    n = len(g); N = 2 * n
    G = np.fft.rfft(g, N)
    k = 2 * np.pi * np.fft.rfftfreq(N, dz)
    G = G * np.exp(-s * k**2 - 1j * k * (2 * A * s))
    if up > 1:
        Gp = np.zeros(up * N // 2 + 1, dtype=complex); Gp[:len(G)] = G
        f = np.fft.irfft(Gp, up * N) * up
        return f[: up * n], dz / up
    return np.fft.irfft(G, N)[:n], dz

def count_interior(f, dz, lo, hi):
    x = np.arange(len(f)) * dz
    m = (x >= lo) & (x <= hi)
    fs = np.sign(f[m]); return int((fs[1:] * fs[:-1] < 0).sum())

def analyse(gamma0, tag=""):
    z0 = 2 * gamma0
    zg = z0 - MARGIN + DZ * np.arange(int((LW + 2 * MARGIN) / DZ) + 1)
    g = g_values(zg)
    lo, hi = MARGIN, MARGIN + LW
    f0, d0 = smooth(g, 0.0, DZ, UP)
    # zero positions at s=0 (spline roots on the original grid)
    cs = CubicSpline(zg, g); x = np.arange(len(zg)) * DZ
    idx = [i for i in range(len(zg) - 1) if lo <= x[i] <= hi and g[i] * g[i + 1] < 0]
    roots = np.array([brentq(cs, zg[i], zg[i + 1]) for i in idx])
    gaps = np.diff(roots)
    # collision region bounded by midpoints, so edge zeros drifting across a fixed edge are not counted
    lo = 0.5 * (roots[0] + roots[1]) - zg[0]; hi = 0.5 * (roots[-2] + roots[-1]) - zg[0]
    c0 = count_interior(f0, d0, lo, hi)
    gaps_in = np.diff(roots[1:-1])
    lemma = (gaps_in**2).min() / 8
    # first collision: scan s, then bisect
    ss = np.arange(0.0005, 0.30, 0.0005); sc = np.nan
    for s in ss:
        f, d = smooth(g, s, DZ, UP)
        if count_interior(f, d, lo, hi) < c0:
            a, b = s - 0.0005, s
            for _ in range(25):
                m = 0.5 * (a + b); f, d = smooth(g, m, DZ, UP)
                if count_interior(f, d, lo, hi) < c0: b = m
                else: a = m
            sc = b; break
    ncoll = []
    for s in (0.05, 0.10, 0.20):
        f, d = smooth(g, s, DZ, UP); ncoll.append(c0 - count_interior(f, d, lo, hi))
    return dict(gamma0=gamma0, n_zeros=c0, n_roots=len(roots), mean_gap=float(gaps.mean()),
                min_gap=float(gaps_in.min()), lemma=float(lemma), s_c=float(sc),
                lost_by=dict(zip(["0.05", "0.10", "0.20"], ncoll)), roots=roots.tolist())

if __name__ == "__main__":
    heights = [float(h) for h in sys.argv[1].split(",")]
    nwin = int(sys.argv[2]); out = sys.argv[3]
    res = []
    for h in heights:
        for w in range(nwin):
            t0 = time.time(); r = analyse(h + w * 40.0)
            r["time"] = time.time() - t0; res.append(r)
            print(f"gamma0={r['gamma0']:.0f} zeros={r['n_zeros']} min_gap={r['min_gap']:.3f} lemma={r['lemma']:.5f} s_c={r['s_c']:.5f} "
                  f"ratio={r['s_c']/r['lemma']:.2f} lost@0.05/0.1/0.2={list(r['lost_by'].values())} ({r['time']:.0f}s)", flush=True)
            json.dump(res, open(out, "w"))
```
