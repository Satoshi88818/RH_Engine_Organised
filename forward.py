"""Forward heat flow (t > 0) of the Riemann Xi function: zero counting with the argument principle.

H_t is evaluated as a Gaussian blur of ξ across σ via Gauss–Hermite quadrature.
Per window: N_real (sign changes on the real axis) and N_tot (winding number).
N_tot − N_real = number of non-real zeros in the strip.
"""

from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import mpmath as mp
import numpy as np
from scipy.optimize import brentq

mp.mp.dps = 25
NODES = 20
Y = 1.0


def xi(s):
    return 0.5 * s * (s - 1) * mp.exp(-(s / 2) * mp.log(mp.pi) + mp.loggamma(s / 2)) * mp.zeta(s)


class Flow:
    def __init__(self, t: float, n: int = NODES):
        self.t = t
        self.n = n
        u, w = np.polynomial.hermite.hermgauss(n)
        self.delta = np.sqrt(t) * u
        self.w = w / np.sqrt(np.pi)
        self._memo: dict = {}

    def H(self, x, y=0.0):
        T = mp.mpf(x) / 2
        scale = mp.exp(mp.pi * T / 4)
        if y == 0.0 and self.t > 0:
            tot = mp.mpf(0)
            for d, w in zip(self.delta, self.w):
                if d > 0:
                    tot += 2 * w * mp.re(xi(mp.mpc(0.5 + d, T)))
                elif d == 0:
                    tot += w * mp.re(xi(mp.mpc(0.5, T)))
            return float(tot * scale)
        tot = mp.mpc(0)
        for d, w in zip(self.delta, self.w):
            tot += w * xi(mp.mpc(0.5 - y / 2 + d, T))
        v = tot * scale
        return complex(v)

    def Hc(self, x, y):
        key = (round(float(x), 9), round(float(y), 9))
        if key not in self._memo:
            self._memo[key] = self.H(x, y) if y >= 0 else np.conj(self.H(x, -y))
        return self._memo[key]


def real_zeros(F, x0, x1, dx=0.1):
    xs = np.arange(x0, x1 + dx / 2, dx)
    v = np.array([F(x) for x in xs])
    out = []
    for i in range(len(xs) - 1):
        if v[i] * v[i + 1] < 0:
            out.append(brentq(F, xs[i], xs[i + 1], xtol=1e-7))
    return np.array(out), xs, v


def winding(Hfun, pts, maxstep=0.9):
    """Adaptive phase tracking around a closed contour."""
    vals = [Hfun(*p) for p in pts]
    total = 0.0
    worst = 0.0
    stack = []
    for (p1, v1), (p2, v2) in zip(zip(pts[:-1], vals[:-1]), zip(pts[1:], vals[1:])):
        stack.append((p1, v1, p2, v2, 0))
        while stack:
            a, va, b, vb, depth = stack.pop()
            dphi = np.angle(vb / va)
            if abs(dphi) > maxstep and depth < 12:
                m = (0.5 * (a[0] + b[0]), 0.5 * (a[1] + b[1]))
                vm = Hfun(*m)
                stack.append((m, vm, b, vb, depth + 1))
                stack.append((a, va, m, vm, depth + 1))
            else:
                total += dphi
                worst = max(worst, abs(dphi))
    return total / (2 * np.pi), worst


def contour(x0, x1, dxc=0.2, dyc=0.1):
    nx = int(np.ceil((x1 - x0) / dxc))
    ny = int(np.ceil(2 * Y / dyc))
    bottom = [(x, -Y) for x in np.linspace(x0, x1, nx + 1)]
    right = [(x1, y) for y in np.linspace(-Y, Y, ny + 1)][1:]
    top = [(x, Y) for x in np.linspace(x1, x0, nx + 1)][1:]
    left = [(x0, y) for y in np.linspace(Y, -Y, ny + 1)][1:]
    return bottom + right + top + left


def analyse(gamma0: float, ts: list[float], width: float = 16.0, control: bool = True) -> dict:
    z0 = 2 * gamma0
    F0 = lambda x: float(mp.re(xi(mp.mpc(0.5, mp.mpf(x) / 2))) * mp.exp(mp.pi * mp.mpf(x) / 8))
    zr, _, _ = real_zeros(F0, z0 - 1.5, z0 + width + 1.5)
    gaps_mid = 0.5 * (zr[1:] + zr[:-1])
    x0 = float(gaps_mid[0])
    x1 = float(gaps_mid[-1])
    c0 = int(((zr > x0) & (zr < x1)).sum())
    out = dict(
        gamma0=gamma0,
        x0=x0,
        x1=x1,
        c0=c0,
        zeros0=zr[(zr > x0) & (zr < x1)].tolist(),
        per_t=[],
    )
    cont = contour(x0, x1)
    for t in ts:
        t_start = time.time()
        fl = Flow(t)
        zt, _, _ = real_zeros(lambda x: fl.H(x, 0.0), x0, x1)
        w, worst = winding(fl.Hc, cont)
        rec = dict(
            t=t,
            n_real=len(zt),
            winding=float(w),
            worst_step=float(worst),
            zeros=zt.tolist(),
        )
        if control:
            a = 0.5 * (x0 + x1)
            b = 0.3
            wc, _ = winding(
                lambda x, y: fl.Hc(x, y) * (((x - a) + 1j * y) ** 2 + b * b), cont
            )
            rec["control_winding"] = float(wc)
        rec["secs"] = time.time() - t_start
        out["per_t"].append(rec)
    return out


def main() -> None:
    heights = [float(h) for h in sys.argv[1].split(",")]
    nwin = int(sys.argv[2])
    ts = [float(t) for t in sys.argv[3].split(",")]
    outp = sys.argv[4]
    width = float(sys.argv[5]) if len(sys.argv) > 5 else 16.0
    res = []
    for h in heights:
        for k in range(nwin):
            r = analyse(h + k * (width / 2 + 10), ts, width)
            res.append(r)
            line = " | ".join(
                f"t={p['t']}: real={p['n_real']} wind={p['winding']:.2f} "
                f"ctrl={p.get('control_winding', float('nan')):.2f} ({p['secs']:.0f}s)"
                for p in r["per_t"]
            )
            print(f"gamma0={r['gamma0']:.0f} c0={r['c0']} :: {line}", flush=True)
            Path(outp).write_text(json.dumps(res, indent=2))


if __name__ == "__main__":
    main()
