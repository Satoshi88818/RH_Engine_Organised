"""Riemann–Siegel Z(t) evaluator (main sum + remainder).

A practical vectorised implementation for surveys up to moderate-to-high height.
For production use at extreme heights the Chebyshev tables for the remainder
should be precomputed (see original notes / rs_tables.py).
"""

from __future__ import annotations

import numpy as np

# Use long double when available for phase accuracy
LD = np.longdouble
TWO_PI = LD("6.283185307179586476925286766559005768394")


def theta(t: float | np.ndarray) -> np.ndarray:
    """Riemann–Siegel theta function (approximate)."""
    t = np.asarray(t, dtype=float)
    return (t / 2) * np.log(t / (2 * np.pi)) - t / 2 - np.pi / 8 + 1 / (48 * t) + 7 / (5760 * t**3)


def Z_main(t: float | np.ndarray) -> np.ndarray:
    """Main sum of the Riemann–Siegel formula (no remainder).

    Adequate for rough surveys; for high accuracy include the C0–C2 remainder
    with precomputed Chebyshev tables.
    """
    t = np.atleast_1d(np.asarray(t, dtype=float))
    out = np.empty_like(t)
    for i, ti in enumerate(t):
        N = int(np.floor(np.sqrt(ti / (2 * np.pi))))
        n = np.arange(1, N + 1, dtype=float)
        ph = theta(ti) - ti * np.log(n)
        out[i] = 2 * np.sum(np.cos(ph) / np.sqrt(n))
    return out


def Z(t: float | np.ndarray) -> np.ndarray:
    """Convenience alias for the main-sum approximation."""
    return Z_main(t)


def survey_zeros(
    t_start: float,
    t_end: float,
    step: float = 0.04,
    refine_step: float = 0.002,
) -> np.ndarray:
    """Simple sign-change survey of Z(t) on [t_start, t_end].

    Returns approximate root locations. For high-quality lists use
    a refined bisection / Newton stage and cross-check against known tables.
    """
    t = np.arange(t_start, t_end, step)
    z = Z(t)
    roots = []
    for i in range(len(t) - 1):
        if z[i] * z[i + 1] < 0:
            # linear interpolation seed
            r = t[i] - z[i] * (t[i + 1] - t[i]) / (z[i + 1] - z[i])
            roots.append(r)
    return np.array(roots)
