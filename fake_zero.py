"""Fake off-line zero test for the candidate energy E(σ).

Inserts a hypothetical off-line zero quadruple and measures the effect on E and E''.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np

from .energy import load_zeros

T0, N_DEFAULT = 30.0, 300


def F(a: np.ndarray, y: np.ndarray) -> np.ndarray:
    """Antiderivative in a of 0.5 * log(y² + a²)."""
    r2 = y * y + a * a
    lg = np.log(np.where(r2 > 0, r2, 1.0))
    return 0.5 * (a * lg - 2 * a + 2 * np.abs(y) * np.arctan2(a, np.abs(y)))


def E(sig: np.ndarray, beta: np.ndarray, gam: np.ndarray, T1: float) -> np.ndarray:
    y = (sig - 0.5)[:, None] - beta[None, :]
    L = T1 - T0
    mean_log = ((F(T1 - gam, y) - F(T0 - gam, y)) / L).sum(axis=1)
    return mean_log - np.log(np.hypot(0.5 + beta, gam)).sum()


def dE(sig: np.ndarray, beta: np.ndarray, gam: np.ndarray, T1: float) -> np.ndarray:
    """Exact E'(sigma)."""
    y = (sig - 0.5)[:, None] - beta[None, :]
    L = T1 - T0
    at = lambda a: np.sign(y) * np.arctan2(a, np.abs(y))
    return ((at(T1 - gam) - at(T0 - gam)) / L).sum(axis=1)


def zero_set(gam: np.ndarray, fake: tuple[float, float] | None = None):
    """Real zeros (±gamma, beta=0) plus an optional fake quadruple (delta, gamma0)."""
    beta = np.zeros(2 * len(gam))
    g = np.concatenate([gam, -gam])
    if fake:
        d, g0 = fake
        beta = np.concatenate([beta, [d, -d, d, -d]])
        g = np.concatenate([g, [g0, g0, -g0, -g0]])
    return beta, g


def analyse(gam: np.ndarray, T1: float, fake: tuple[float, float] | None, h: float = 0.005) -> dict:
    sig = np.round(np.arange(0, 1 + h / 2, h), 10)
    b, g = zero_set(gam, fake)
    e = E(sig, b, g, T1)
    epp = np.gradient(np.gradient(e, h), h)
    far = np.abs(sig - 0.5) > 0.012
    thr = 20 * np.median(np.abs(epp[far]))
    peaks = [
        s
        for i, s in enumerate(sig[1:-1], 1)
        if far[i] and epp[i] > thr and epp[i] >= epp[i - 1] and epp[i] >= epp[i + 1]
    ]
    jumps = []
    if fake:
        b0, g0_ = zero_set(gam)
        for c in (0.5 + fake[0], 0.5 - fake[0]):
            w = 2 * h
            s2 = np.array([c - w, c + w])
            jumps.append(float(np.diff(dE(s2, b, g, T1) - dE(s2, b0, g0_, T1))[0]))
    return dict(
        sig=sig,
        E=e,
        Epp=epp,
        peaks=peaks,
        jumps=jumps,
        argmin=float(sig[e.argmin()]),
        minEpp=float(epp[3:-3].min()),
        symerr=float(np.abs(e - e[::-1]).max()),
    )


def main(n_zeros: int = N_DEFAULT, out_png: str = "rh_fake_zero.png") -> None:
    gam = load_zeros(n_zeros)
    T1 = 0.5 * gam[-1]
    L = T1 - T0
    print(f"window [{T0}, {T1:.1f}], L={L:.1f}, predicted mass per cusp 2*pi/L = {2 * np.pi / L:.4f}")
    base = analyse(gam, T1, None)
    print(f"baseline: argmin={base['argmin']}, off-axis peaks={base['peaks']}, min E''={base['minEpp']:.3f}")

    print("\nTable 1: delta sweep, gamma0=150.3")
    for d in (0.02, 0.05, 0.10, 0.20, 0.30, 0.40):
        r = analyse(gam, T1, (d, 150.3))
        print(
            f"delta={d:.2f} argmin={r['argmin']:.3f} blind peaks={[round(p, 3) for p in r['peaks']]} "
            f"jumps={[round(j, 4) for j in r['jumps']]} minE''={r['minEpp']:.2f}"
        )

    print("\nTable 2: gamma0 sweep, delta=0.20")
    for g0 in (40.0, 100.0, 150.3, 250.0, 400.0):
        r = analyse(gam, T1, (0.20, g0))
        print(
            f"gamma0={g0:6.1f} blind peaks={[round(p, 3) for p in r['peaks']]} "
            f"jumps={[round(j, 4) for j in r['jumps']]}"
        )

    sig = base["sig"]
    r = analyse(gam, T1, (0.20, 150.3))
    dE_ = r["E"] - base["E"]
    print(
        f"\nDelta E (fake - base): at 1/2 = {dE_[100]:.3f}, min={dE_.min():.3f}, "
        f"max={dE_.max():.3f}, all finite={np.isfinite(r['E']).all()}"
    )

    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        fig, ax = plt.subplots(1, 2, figsize=(11, 3.8))
        m = np.abs(sig - 0.5) > 0.012
        ax[0].plot(sig[m], base["Epp"][m], label="no fake zero", lw=1)
        ax[0].plot(sig[m], r["Epp"][m], label="fake zeros at 1/2±0.2", lw=1)
        ax[0].set(xlabel="sigma", title="E''(sigma) off the axis", ylim=(-1, 7))
        ax[0].legend()
        ax[1].plot(sig, base["E"], label="baseline")
        ax[1].plot(sig, r["E"], "--", label="with fake zeros")
        ax[1].set(xlabel="sigma", title="E(sigma): minimum stays at 1/2")
        ax[1].legend()
        plt.tight_layout()
        plt.savefig(out_png, dpi=130)
        print(f"saved {out_png}")
    except ImportError:
        pass


if __name__ == "__main__":
    import sys

    n = int(sys.argv[1]) if len(sys.argv) > 1 else N_DEFAULT
    main(n)
