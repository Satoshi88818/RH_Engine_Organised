# Forward flow (t > 0): zero counting for H_t with the argument principle

## Purpose

`README_harness.md` could only probe t ≤ 0, where the flow is a Gaussian smoothing along the real axis. An off-line zero of the true ζ would show up on the other side, at t > 0. This build is the instrument for that side. It counts **all** zeros of H_t in a window of the strip, real or not, and compares with the count of real ones.

It is a working tool and a consistency check, **not a proof**, and in the range tested it can only confirm what is already known.

## Method

**Exact representation.** From H_t(z) = ½∫Φ(u)e^{tu²}e^{izu}du and e^{tu²} = E[e^{√(2t)·u·W}] with W ~ N(0,1):

> H_t(x + iy) ∝ E_δ[ ξ(½ − y/2 + δ + i·x/2) ], δ ~ N(0, t/2), with z = 2γ.

So H_t is a Gaussian blur of ξ across σ, and a complex point z = x + iy corresponds to the point σ = ½ − y/2 in the strip. I evaluate the average by 20-point Gauss–Hermite quadrature, with `mpmath` ξ at 25 digits. H_t is real on the real axis; I use conjugate-paired nodes there and conjugate symmetry for y < 0.

**Per window** (width 16 in z, about 1–7 zeros):
1. Find the real zeros at t = 0 and put the left/right edges at midpoints of the outer gaps. This gives c₀ zeros inside.
2. For each t: count real sign changes on the real axis (N_real), and compute the **winding number** of H_t around the rectangle [x₀, x₁] × [−1, 1] (N_tot), with adaptive phase tracking (every phase step ≤ 0.9 rad; worst observed 0.45). The strip |Im z| < 1 is the full range 0 < σ < 1.
3. **Non-real zeros in the window = N_tot − N_real.**

**Positive control.** I multiply H_t by ((x − a) + iy)² + b² with b = 0.3, which inserts a hand-made pair a ± ib inside the window. The winding number must rise by exactly 2 and N_real must not change.

## Validation

- **Against the ODE flow of the zeros.** At γ ≈ 100, zero positions at t = 0.05 and 0.2 differ from my 300-zero ODE (`rh_flow_sketch.py`) by 0.017 and 0.068. That entire gap is the missing tail zeros beyond the truncation: adding the leading-order tail force (−4x·ρ-integral, about −0.0016·x per unit time) reduces the mismatch to **2 × 10⁻⁴ and 6 × 10⁻⁴**.
- **Quadrature.** At γ ≈ 1000, t = 0.2: switching from 20 to 32 nodes shifts three zeros by less than 10⁻⁸.
- **Positive control:** in all 32 (window, t) runs the control raised the winding number by exactly 2.00.

## Results

8 windows, 36 zeros in total, t ∈ {0.02, 0.05, 0.1, 0.2}, strip 0 < σ < 1.

| γ | windows | zeros | t values | non-real zeros found | control detected | N_real = c₀ at all t | worst phase step | CPU time |
|---|---|---|---|---|---|---|---|---|
| ~60–90 | 2 | 3 | 4 | 0 | +2.00 in 8/8 | yes | 0.23 | 29 s |
| ~300–330 | 2 | 8 | 4 | 0 | +2.00 in 8/8 | yes | 0.29 | 234 s |
| ~1000–1030 | 2 | 11 | 4 | 0 | +2.00 in 8/8 | yes | 0.37 | 438 s |
| ~3000–3030 | 2 | 14 | 4 | 0 | +2.00 in 8/8 | yes | 0.45 | 829 s |

- **No non-real zeros** in any window at any t tested; winding = N_real = c₀ everywhere, so no zero crossed the window boundary either.
- **Zeros drift as the flow predicts.** At γ ≈ 3000 the zeros move by up to 0.44 in z between t = 0 and t = 0.2, with neighbouring zeros moving in opposite directions (repulsion).
- Cost is about 25 s per (window, t) at γ ≈ 300 and 100 s at γ ≈ 3000, dominated by off-line ζ evaluations.

## What this does and does not show

- It confirms the instrument works: it sees a planted pair every time and reports none otherwise.
- **The result is expected, not new.** All zeros of H_t are known to be real for t ≥ 0.22 (Polymath15), RH is numerically verified far beyond these heights, and zeros cannot go non-real at t > 0 once real. So finding none here tells us nothing about RH.
- **Where it could matter:** at heights beyond direct verification, using this same winding test with a fast ξ approximation. That is the Polymath barrier idea in miniature.
- **Sensitivity.** A genuine off-line pair at offset b₀ (in z) would, by the two-body estimate of `README_phase_velocity_map.md`, collide near t ≈ b₀²/2. Windows tested at t = 0.02 would therefore see any pair with b₀ ≳ 0.2, i.e. σ-offset ≳ 0.1.

## Limitations

- Only 8 windows, 36 zeros, heights ≤ 3000, widths 16. This is a demonstration, not a survey.
- Floating-point, finite-sampling numerics. It is not an interval-arithmetic or rigorous verification.
- The strip is |Im z| < 1; zeros further out would be missed (none are expected).
- Window edges are fixed at t = 0 zero gaps; I check N_real = c₀ at each t to catch zeros crossing them. All runs passed.
- Scaling to large heights needs a fast ξ approximation (A + B terms), not mpmath ζ.

## Next steps

1. Replace the mpmath ξ with the Polymath A + B approximation so heights of 10⁶ and beyond are feasible.
2. Add small t (down to 10⁻³) and track how the winding count behaves as t → 0 at larger heights.
3. Use the lifetime bounds from `README_no_collision_bound.md` to choose t(x) and window sizes automatically.

## Reproduce

```bash
pip install mpmath numpy scipy
python rh_forward.py 60,300,1000,3000 2 0.02,0.05,0.1,0.2 forward_results.json 16
```

About 25 minutes on one core.

## Script

```python
"""Forward heat flow (t > 0) of the Riemann Xi function: zero counting with the argument principle.

Exact representation (derived from H_t(z) = (1/2)∫Φ(u)e^{tu²}e^{izu}du and e^{tu²} = E[e^{√(2t)·u·W}]):
    H_t(x+iy) ∝ E_δ[ ξ(1/2 - y/2 + δ + i x/2) ],  δ ~ N(0, t/2),   z = 2γ.
i.e. a Gaussian blur of ξ across σ.  Evaluated by Gauss-Hermite quadrature.
Per window [x0,x1] x [-Y,Y]: N_real = sign changes on the real axis, N_tot = winding number of H_t on the boundary.
N_tot - N_real = number of NON-REAL zeros of H_t in the window (strip |Im z| < Y, i.e. 0 < σ < 1 for Y = 1).
"""
import sys, json, time
import numpy as np, mpmath as mp
from scipy.optimize import brentq

mp.mp.dps = 25
NODES = 20
Y = 1.0

def xi(s):
    return 0.5 * s * (s - 1) * mp.exp(-(s / 2) * mp.log(mp.pi) + mp.loggamma(s / 2)) * mp.zeta(s)

class Flow:
    def __init__(self, t, n=NODES):
        self.t = t; self.n = n
        u, w = np.polynomial.hermite.hermgauss(n)
        self.delta = np.sqrt(t) * u; self.w = w / np.sqrt(np.pi); self._memo = {}
    def H(self, x, y=0.0):
        T = mp.mpf(x) / 2; scale = mp.exp(mp.pi * T / 4)
        if y == 0.0 and self.t > 0:                      # symmetric nodes: use conjugate pairs
            tot = mp.mpf(0)
            for d, w in zip(self.delta, self.w):
                if d > 0: tot += 2 * w * mp.re(xi(mp.mpc(0.5 + d, T)))
                elif d == 0: tot += w * mp.re(xi(mp.mpc(0.5, T)))
            return float(tot * scale)
        tot = mp.mpc(0)
        for d, w in zip(self.delta, self.w):
            tot += w * xi(mp.mpc(0.5 - y / 2 + d, T))
        v = tot * scale
        return complex(v)
    def Hc(self, x, y):                                    # any sign of y (H real on the axis)
        key = (round(float(x), 9), round(float(y), 9))
        if key not in self._memo:
            self._memo[key] = self.H(x, y) if y >= 0 else np.conj(self.H(x, -y))
        return self._memo[key]

def real_zeros(F, x0, x1, dx=0.1):
    xs = np.arange(x0, x1 + dx / 2, dx); v = np.array([F(x) for x in xs])
    out = []
    for i in range(len(xs) - 1):
        if v[i] * v[i + 1] < 0:
            out.append(brentq(F, xs[i], xs[i + 1], xtol=1e-7))
    return np.array(out), xs, v

def winding(Hfun, pts, maxstep=0.9):
    """pts: list of (x,y) around a closed contour (first == last). Adaptive phase tracking."""
    vals = [Hfun(*p) for p in pts]
    path = [(pts[0], vals[0])]; total = 0.0; worst = 0.0
    stack = []
    for (p1, v1), (p2, v2) in zip(zip(pts[:-1], vals[:-1]), zip(pts[1:], vals[1:])):
        stack.append((p1, v1, p2, v2, 0))
        while stack:
            a, va, b, vb, depth = stack.pop()
            dphi = np.angle(vb / va)
            if abs(dphi) > maxstep and depth < 12:
                m = (0.5 * (a[0] + b[0]), 0.5 * (a[1] + b[1])); vm = Hfun(*m)
                stack.append((m, vm, b, vb, depth + 1)); stack.append((a, va, m, vm, depth + 1))
            else:
                total += dphi; worst = max(worst, abs(dphi))
    return total / (2 * np.pi), worst

def contour(x0, x1, dxc=0.2, dyc=0.1):
    nx = int(np.ceil((x1 - x0) / dxc)); ny = int(np.ceil(2 * Y / dyc))
    bottom = [(x, -Y) for x in np.linspace(x0, x1, nx + 1)]
    right = [(x1, y) for y in np.linspace(-Y, Y, ny + 1)][1:]
    top = [(x, Y) for x in np.linspace(x1, x0, nx + 1)][1:]
    left = [(x0, y) for y in np.linspace(Y, -Y, ny + 1)][1:]
    return bottom + right + top + left

def analyse(gamma0, ts, width=16.0, control=True):
    z0 = 2 * gamma0
    F0 = Flow(0.0, 1).H if False else (lambda x: float(mp.re(xi(mp.mpc(0.5, mp.mpf(x) / 2))) * mp.exp(mp.pi * mp.mpf(x) / 8)))
    zr, _, _ = real_zeros(F0, z0 - 1.5, z0 + width + 1.5)
    gaps_mid = 0.5 * (zr[1:] + zr[:-1])
    x0 = gaps_mid[0]; x1 = gaps_mid[-1]; c0 = int(((zr > x0) & (zr < x1)).sum())
    out = dict(gamma0=gamma0, x0=float(x0), x1=float(x1), c0=c0, zeros0=zr[(zr > x0) & (zr < x1)].tolist(), per_t=[])
    cont = contour(x0, x1)
    for t in ts:
        t_start = time.time(); fl = Flow(t)
        zt, _, _ = real_zeros(lambda x: fl.H(x, 0.0), x0, x1)
        w, worst = winding(fl.Hc, cont)
        rec = dict(t=t, n_real=len(zt), winding=float(w), worst_step=float(worst), zeros=zt.tolist())
        if control:                                          # positive control: insert a non-real pair a ± i b
            a = 0.5 * (x0 + x1); b = 0.3
            wc, _ = winding(lambda x, y: fl.Hc(x, y) * (((x - a) + 1j * y) ** 2 + b * b), cont)
            rec["control_winding"] = float(wc)
        rec["secs"] = time.time() - t_start; out["per_t"].append(rec)
    return out

if __name__ == "__main__":
    heights = [float(h) for h in sys.argv[1].split(",")]; nwin = int(sys.argv[2])
    ts = [float(t) for t in sys.argv[3].split(",")]; outp = sys.argv[4]
    width = float(sys.argv[5]) if len(sys.argv) > 5 else 16.0
    res = []
    for h in heights:
        for k in range(nwin):
            r = analyse(h + k * (width / 2 + 10), ts, width); res.append(r)
            line = " | ".join(f"t={p['t']}: real={p['n_real']} wind={p['winding']:.2f} ctrl={p.get('control_winding', float('nan')):.2f} ({p['secs']:.0f}s)" for p in r["per_t"])
            print(f"gamma0={r['gamma0']:.0f} c0={r['c0']} :: {line}", flush=True)
            json.dump(res, open(outp, "w"))
```
