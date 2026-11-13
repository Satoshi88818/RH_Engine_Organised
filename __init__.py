"""RH-Engine: numerical experiments around the Riemann Hypothesis.

Modules
-------
energy
    Candidate energy functional E(σ) from the Hadamard product over zeros.
fake_zero
    Insert hypothetical off-line zeros and measure the effect on E / E''.
harness
    Backward heat-flow collision-time harness (gap-lifetime lemma).
forward
    Forward heat flow (t > 0) with winding-number zero counting.
riemann_siegel
    Practical Riemann–Siegel Z(t) evaluator and simple zero survey.

None of the code constitutes a proof of the Riemann Hypothesis.
All numerics are floating-point and non-rigorous.
"""

__version__ = "0.1.0"

from .energy import compute_energy_profile, energy, load_zeros
from .fake_zero import analyse as analyse_fake_zero
from .harness import analyse as analyse_harness, run_harness

__all__ = [
    "__version__",
    "load_zeros",
    "energy",
    "compute_energy_profile",
    "analyse_fake_zero",
    "analyse_harness",
    "run_harness",
]
