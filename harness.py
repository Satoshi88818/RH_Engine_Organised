"""Empirical harness: local collision times under backward heat flow of Xi.

Tests the gap-lifetime lemma (a real pair stays apart for at least r²/8)
on windows of real zeros up to moderate height.
"""

from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import mpmath as mp
import numpy as np
from scipy.interpolate import CubicSpline
from scipy.optimize import brentq

mp.mp.dps = 20
A = np.pi / 8  # tilt that cancels the e^{-π z/8} decay of H₀
DZ, MARGIN, LW = 0.1, 3.6, 40.0
UP = 4


def g_values(zgrid: np.ndarray) -> np.ndarray:
    out = np.empty(len(zgrid))
    for i, z in enumerate(zgrid):
        t = mp.mpf(z) / 2
        lg = mp.re(mp.loggamma(mp.mpc(0.25, t / 2)))
        out[i] = float(mp.siegelz(t) * (mp.mpf(1) / 4 + t * t) * mp.exp(lg + A * z))
    return out


def smooth(g: np.ndarray, s: float, dz: float, up: int = 1):
    """Gaussian smoothing corresponding to H_{-s}."""
    n = len(g)
    N = 2 * n
    G = np.fft.rfft(g, N)
    k = 2 * np.pi * np.fft.rfftfreq(N, dz)
    G = G * np.exp(-s * k**2 - 1j * k * (2 * A * s))
    if up > 1:
        Gp = np.zeros(up * N // 2 + 1, dtype=complex)
        Gp[: len(G)] = G
        f = np.fft.irfft(Gp, up * N) * up
        return f[: up * n], dz / up
    return np.fft.irfft(G, N)[:n], dz


def count_interior(f: np.ndarray, dz: float, lo: float, hi: float) -> int:
    x = np.arange(len(f)) * dz
    m = (x >= lo) & (x <= hi)
    fs = np.sign(f[m])
    return int((fs[1:] * fs[:-1] < 0).sum())


def analyse(gamma0: float) -> dict:
    z0 = 2 * gamma0
    zg = z0 - MARGIN + DZ * np.arange(int((LW + 2 * MARGIN) / DZ) + 1)
    g = g_values(zg)
    lo, hi = MARGIN, MARGIN + LW
    f0, d0 = smooth(g, 0.0, DZ, UP)

    cs = CubicSpline(zg, g)
    x = np.arange(len(zg)) * DZ
    idx = [i for i in range(len(zg) - 1) if lo <= x[i] <= hi and g[i] * g[i + 1] < 0]
    roots = np.array([brentq(cs, zg[i], zg[i + 1]) for i in idx])
    gaps = np.diff(roots)

    # Collision region bounded by midpoints to avoid edge-drift artefacts
    lo = 0.5 * (roots[0] + roots[1]) - zg[0]
    hi = 0.5 * (roots[-2] + roots[-1]) - zg[0]
    c0 = count_interior(f0, d0, lo, hi)
    gaps_in = np.diff(roots[1:-1])
    lemma = float((gaps_in**2).min() / 8)

    ss = np.arange(0.0005, 0.30, 0.0005)
    sc = np.nan
    for s in ss:
        f, d = smooth(g, s, DZ, UP)
        if count_interior(f, d, lo, hi) < c0:
            a, b = s - 0.0005, s
            for _ in range(25):
                m = 0.5 * (a + b)
                f, d = smooth(g, m, DZ, UP)
                if count_interior(f, d, lo, hi) < c0:
                    b = m
                else:
                    a = m
            sc = b
            break

    ncoll = []
    for s in (0.05, 0.10, 0.20):
        f, d = smooth(g, s, DZ, UP)
        ncoll.append(c0 - count_interior(f, d, lo, hi))

    return {
        "gamma0": gamma0,
        "n_zeros": c0,
        "n_roots": len(roots),
        "mean_gap": float(gaps.mean()),
        "min_gap": float(gaps_in.min()),
        "lemma": lemma,
        "s_c": float(sc),
        "lost_by": dict(zip(["0.05", "0.10", "0.20"], ncoll)),
        "roots": roots.tolist(),
    }


def run_harness(
    heights: list[float],
    n_windows: int = 6,
    out_path: str | Path = "results.json",
) -> list[dict]:
    res = []
    for h in heights:
        for w in range(n_windows):
            t0 = time.time()
            r = analyse(h + w * 40.0)
            r["time"] = time.time() - t0
            res.append(r)
            ratio = r["s_c"] / r["lemma"] if r["lemma"] > 0 else float("nan")
            print(
                f"gamma0={r['gamma0']:.0f} zeros={r['n_zeros']} "
                f"min_gap={r['min_gap']:.3f} lemma={r['lemma']:.5f} "
                f"s_c={r['s_c']:.5f} ratio={ratio:.2f} "
                f"({r['time']:.0f}s)",
                flush=True,
            )
            Path(out_path).write_text(json.dumps(res, indent=2))
    return res


def main() -> None:
    heights = [float(h) for h in sys.argv[1].split(",")]
    nwin = int(sys.argv[2]) if len(sys.argv) > 2 else 6
    out = sys.argv[3] if len(sys.argv) > 3 else "results.json"
    run_harness(heights, nwin, out)


if __name__ == "__main__":
    main()
