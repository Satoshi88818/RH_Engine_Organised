# Riemann–Siegel build: lifetime scaling of zero pairs up to γ ≈ 10⁹

## Purpose

`README_harness.md` was limited to γ ≲ 10⁶ and about 650 zeros by the cost of `mpmath` (25 ms per Z evaluation). This build replaces it with a vectorised Riemann–Siegel Z(t), so that I can (1) survey **millions** of zeros, (2) test the gap-lifetime lemma on the *most extreme* close pairs, and (3) test the scaling of the collision time with height against the random-matrix (GUE) prediction.

It is evidence about how "barely true" RH would be, **not a proof**, and nothing here bears on whether a zero is off the line.

## 1. The Riemann–Siegel evaluator

> Z(t) = 2 Σ_{n ≤ N} n^{−1/2} cos(θ(t) − t ln n) + (−1)^{N−1}(2π/t)^{1/4}[C₀(p) + C₁(p)(2π/t)^{1/2} + C₂(p)(2π/t)], N = ⌊√(t/2π)⌋, p = √(t/2π) − N.

- **C₀, C₁, C₂** come from Φ(p) = cos(2π(p² − p − 1/16))/cos(2πp) and its derivatives, computed by Taylor series in 60-digit arithmetic and stored as Chebyshev tables (`rs_tables.py`, table error 5 × 10⁻¹⁵).
- **Phases** are formed in 80-bit long double, reduced mod 2π, and then expanded around the window centre to 4th order. Z on a uniform grid then costs one complex multiply per term (a recurrence), so a 2500-point window takes 5 ms at γ = 10⁵ and 27 ms at 10⁹.
- Keeping t as (centre, small offset) matters: at γ ≈ 10⁹, double precision has spacing 1.2 × 10⁻⁷, which alone would cost 10⁻⁶ accuracy.

**Accuracy vs `mpmath.siegelz` at exact t:**

| γ | max error |
|---|---|
| 2 × 10³ | 1 × 10⁻⁶ |
| 10⁵ | 6 × 10⁻¹² |
| 10⁷ | 2 × 10⁻¹¹ |
| 10⁹ | 5 × 10⁻⁹ |
| 10¹¹ | 1 × 10⁻⁶ |

## 2. Zero survey

`rh_survey.py` samples Z on a grid of 0.04 and finds roots by bracketed bisection on a cubic spline. Wherever two adjacent intervals both contain a sign change (zeros closer than 2h), or a local extremum is suspiciously close to zero, it re-samples at 0.002 and replaces the coarse roots there.

**Checks.** (a) The number of zeros found matches the Riemann–von Mangoldt count (θ(T₁) − θ(T₀))/π to within ±1.8 in every one of 7 runs (2.3 million zeros), so no zeros are missed or duplicated beyond the known small fluctuation. (b) For the closest pairs, `mpmath` shows two sign changes and |Z| ≈ 10⁻¹¹ at both ends.

**Two bugs found and fixed during the build:**
- The first root solver (Newton) clipped to an interval end for pairs straddling a grid point. It produced fake "pairs" with gaps near 10⁻⁵, far below anything GUE allows. `mpmath` showed only one sign change there. Replacing it with bisection removed them.
- Flow windows with zero-padded edges showed fake early collisions: the edge discontinuity creates Gibbs ripples comparable to the shallow dip between a very close pair. Tapering the window edges (Tukey) fixed it.

**Results** (N consecutive zeros starting at the listed height; ν = gap × local mean density; P = probability under GUE that the minimum would be at most as small as observed):

| start γ | zeros N | span ends | ⟨ln(γ/2π)⟩ | count − RvM | min ν | GUE median min ν | P | min gap (γ) | lifetime min = gap²/2 |
|---|---|---|---|---|---|---|---|---|---|
| 10³ | 15,064 | 1.5 × 10⁴ | 7.0 | +0.3 | 0.0421 | 0.0347 | 0.71 | 0.0377 | 7.1 × 10⁻⁴ |
| 10⁴ | 60,043 | 5.5 × 10⁴ | 8.5 | +0.9 | 0.0407 | 0.0219 | 0.99 | 0.0295 | 4.4 × 10⁻⁴ |
| 10⁵ | 200,088 | 2.2 × 10⁵ | 10.1 | −0.8 | 0.0265 | 0.0147 | 0.98 | 0.0163 | 1.3 × 10⁻⁴ |
| 10⁶ | 500,098 | 1.3 × 10⁶ | 12.1 | −0.4 | 0.0130 | 0.0108 | 0.70 | 0.00679 | 2.3 × 10⁻⁵ |
| 10⁷ | 500,027 | 1.0 × 10⁷ | 14.3 | 0.0 | 0.0146 | 0.0108 | 0.82 | 0.00640 | 2.0 × 10⁻⁵ |
| 10⁸ | 500,109 | 1.0 × 10⁸ | 16.6 | −1.8 | 0.0097 | 0.0108 | 0.39 | 0.00366 | 6.7 × 10⁻⁶ |
| 10⁹ | 500,092 | 1.0 × 10⁹ | 18.9 | +0.9 | 0.0115 | 0.0108 | 0.56 | 0.00381 | 7.3 × 10⁻⁶ |

The GUE prediction used here is the small-gap law: the fraction of gaps with ν < x is F(x) ≈ (π²/9)x³ − (2π⁴/225)x⁵, so the median minimum over N gaps is ν ≈ (9 ln 2/(π² N))^{1/3}. This is the standard small-gap expansion (my recollection of the second coefficient is not independently verified, but it affects only the second digit).

**Small-gap counts, observed / GUE expected for ν < 0.1:** 8/16, 44/65, 178/218, 450/544, 515/544, 481/544, 499/544 (ratios 0.50, 0.68, 0.82, 0.83, 0.95, 0.88, 0.92). So at lower heights there are fewer very close pairs than GUE predicts (stronger repulsion), and the deficit shrinks as γ grows. This is the known slow convergence of zeta zeros to GUE, which I did not model. At 10⁷–10⁹ the counts are within about 10% of GUE.

## 3. The gap-lifetime lemma on extreme pairs and typical windows

`rh_flow2.py` runs the exact backward flow (Gaussian smoothing, as in `README_harness.md`) around the 3 closest pairs at each of the 7 heights (21 pairs, ν as small as 0.01), and on 30 random windows of width 40 (in z) at each height from 10⁵ to 10⁹ (150 windows).

- **Closest pairs:** s_c / (r²/8) ranges from **1.00008 to 1.0059** (median 1.0002). For an isolated close pair the lemma is essentially an equality. The collision time is r²/8 to within 0.6%, at every height up to 10⁹.
- **Typical windows:** s_c / (min r²/8) has median 1.03–1.06 and range 0.98–1.49. The slack is the tidal correction when the closest pair is not very close. One window (10⁷) had ratio 0.02 at grid 0.05 because of an unresolved close pair; with grid 0.01 it gave 1.001 and is reported that way. The 0.98 at 10⁹ is within the root-resolution error.
- **Typical-window scaling:** median s_c × ln²(γ/2π) = 2.1, 1.65, 1.57, 1.60, 1.09 at γ = 10⁵ … 10⁹. This is roughly constant, consistent with a 1/log²γ law for typical windows, though the 10⁹ value is lower and 30 windows per height is noisy.

![Closest pair vs GUE prediction; lemma vs measurement up to 10⁹](rh_rs_scaling.png)

## 4. What this says about the t(x) problem

1. **For the lifetime of the closest pair in a stretch of N zeros,** lifetime ≈ r²/8 = (gap_γ)²/2 with gap_γ = ν · 2π/ln(γ/2π). With the GUE small-gap law the median is

> τ_min(N, γ) ≈ (2π²/ln²(γ/2π)) · (9 ln 2 / (π² N))^{2/3}.

It matches the observed minima at γ ≥ 10⁶ (e.g., 7 × 10⁻⁶ observed vs 6 × 10⁻⁶ predicted at 10⁹).

2. **This corrects the conjecture in `README_no_collision_bound.md`.** The 1/log²x guess is right for *typical* windows, but if t(x) must hold for every pair among the N ≈ xL/2π zeros near x, the binding constraint is the closest pair, and τ_min falls like (log x)^{−2} N^{−2/3}, roughly x^{−2/3} times powers of log. A no-collision window criterion must shrink at that rate, much faster than 1/log²x.

3. **Extrapolation (heuristic, assumes GUE small-gap statistics persist):** using all zeros up to X, the median lifetime of the closest pair is about 2 × 10⁻⁸ at 10⁹, 9 × 10⁻¹¹ at 10¹², 1.5 × 10⁻¹¹ at 10¹³, 5 × 10⁻¹³ at 10¹⁵, and 10⁻¹⁶ at 10²⁰. Since the closest pair at height X gives Λ ≳ −τ_min, this suggests Λ lower bounds of order −10⁻¹¹ at γ ≈ 10¹³. That is the same order as the best published lower bound (about −10⁻¹¹, from memory; I did not re-verify it). The lower-height deficit of close pairs means the true values at moderate heights are somewhat smaller in magnitude than this extrapolation.

4. **Meaning for RH:** RH is "barely true" in this quantitative sense. If RH holds then Λ = 0, and the first collision time in the backward flow tends to 0 like a power of the height. No fixed margin below t = 0 exists.

## 5. What this does not show

- No claim about any zero being off the line, or about t > 0 beyond `README_forward_flow.md`.
- The extrapolation rests on GUE statistics persisting; it is not proved, and zeta zeros at finite height differ measurably (the small-gap deficit above).
- The survey covers 2.3 million zeros out of an unbounded set, at 7 start heights. Minima are extreme-value statistics, so each height is one sample.
- The lemma is exact for finite real systems; the infinite extension is still assumed.

## 6. Next steps

- Push to γ ≈ 10¹¹ (evaluator error 10⁻⁶ there, which needs care for close pairs).
- Compare the closest-pair list with Odlyzko's published zero tables up to 10²².
- Test the (WNC) idea directly: choose t(x) ~ x^{−2/3} and check, window by window, that backward flow keeps all zeros real until t(x).

## Reproduce

```bash
pip install mpmath numpy scipy matplotlib
python rs_tables.py                       # Chebyshev tables
python rh_survey.py 1e3,1e4,1e5,1e6,1e7,1e8,1e9 1.5e4,6e4,2e5,5e5,5e5,5e5,5e5 survey_results.json   # ~3 min
python rh_flow2.py                        # ~3 min
```

Files: `rs_z.py`, `rs_tables.py`, `rs_cheb_coef.npy`, `rh_survey.py`, `rh_flow2.py`, `survey_results.json`, `flow2_results.json`, `rh_rs_scaling.png`. `rh_flow2.py` imports the smoothing routines from `rh_harness.py`.

## Core evaluator (`rs_z.py`)

```python
"""Vectorised Riemann-Siegel Z(t): main sum + C0,C1,C2 remainder.  Phases in long double."""
import numpy as np
LD = np.longdouble
_COEF = np.load("rs_cheb_coef.npy"); _A, _B = -0.15, 1.15
PI_L = LD(np.pi) if np.finfo(LD).eps < 1e-18 else LD(np.pi)
TWO_PI_L = LD("6.283185307179586476925286766559005768394")

def theta_ld(t):
    t = LD(t)
    return (t / 2) * np.log(t / TWO_PI_L) - t / 2 - TWO_PI_L / 16 + 1 / (48 * t) + 7 / (5760 * t**3)

def corr(t):
    """RS remainder (without the main sum)."""
    t = np.asarray(t, float); r = np.sqrt(t / (2 * np.pi)); N = np.floor(r); p = r - N
    x = (2 * p - (_A + _B)) / (_B - _A)
    C = [np.polynomial.chebyshev.chebval(x, _COEF[j]) for j in range(3)]
    q = (2 * np.pi / t)
    sgn = np.where(N % 2 == 0, -1.0, 1.0)             # (-1)^(N-1)
    return sgn * q**0.25 * (C[0] + C[1] * q**0.5 + C[2] * q)

def Z(t):
    """Direct evaluation at arbitrary (array) t."""
    t = np.atleast_1d(np.asarray(t, dtype=float)); out = np.empty_like(t)
    for i, ti in enumerate(t):
        N = int(np.floor(np.sqrt(ti / (2 * np.pi)))); n = np.arange(1, N + 1, dtype=LD)
        ph = theta_ld(ti) - LD(ti) * np.log(n)
        ph = np.mod(ph, TWO_PI_L).astype(float)
        out[i] = 2 * np.sum(np.cos(ph) / np.sqrt(n.astype(float))) + corr(ti)
    return out


def corr_N(t, N):
    """RS remainder with a FIXED N (p may leave [0,1] slightly near an N-boundary; tables cover [-0.15, 1.15])."""
    t = np.asarray(t, float); p = np.sqrt(t / (2 * np.pi)) - N
    x = (2 * p - (_A + _B)) / (_B - _A)
    C = [np.polynomial.chebyshev.chebval(x, _COEF[j]) for j in range(3)]
    q = 2 * np.pi / t
    return (-1.0 if N % 2 == 0 else 1.0) * q**0.25 * (C[0] + C[1] * q**0.5 + C[2] * q)

def Z_uniform(tc, h, M, N=None):
    """Z at t_k = tc + h*(k - (M-1)/2), k < M  (tc exactly representable; offsets handled relative to tc).
    Taylor-expanded phases + complex recurrence (fast)."""
    if N is None:
        N = int(np.floor(np.sqrt(tc / (2 * np.pi))))
    n = np.arange(1, N + 1, dtype=LD)
    b = np.mod(theta_ld(tc) - LD(tc) * np.log(n), TWO_PI_L).astype(float)
    lnn = np.log(np.arange(1, N + 1, dtype=float))
    thp = 0.5 * np.log(tc / (2 * np.pi)) - 1 / (48 * tc**2)
    w = thp - lnn
    d0 = -h * (M - 1) / 2
    u = np.arange(1, N + 1, dtype=float) ** -0.5 * np.exp(1j * (b + d0 * w))
    step = np.exp(1j * h * w)
    S = np.empty(M, complex)
    for k in range(M):
        S[k] = u.sum(); u *= step
    d = d0 + h * np.arange(M)
    thpp = 1 / (2 * tc) + 1 / (24 * tc**3); th3 = -1 / (2 * tc**2); th4 = 1 / tc**3
    c = 0.5 * thpp * d**2 + th3 * d**3 / 6 + th4 * d**4 / 24
    t = tc + d
    return 2 * np.real(np.exp(1j * c) * S) + corr_N(t, N)
```
