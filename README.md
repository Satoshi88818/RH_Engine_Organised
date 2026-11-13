# RH-Engine

**Riemann Hypothesis Engine** — a cleaned, organised Python package of numerical experiments and theoretical notes exploring objects connected to the Riemann Hypothesis (RH).

> **Not a proof and not a claim — just an engine.**

This repository reorganises the original flat collection of research notes and embedded scripts (November 2026) into a maintainable Python package. All numerics are floating-point and non-rigorous. Nothing here constitutes evidence that the Riemann Hypothesis is true or false.

---

## What this package does

The material began as an exploration (and later critical assessment) of speculative variational approaches to RH, then shifted to precise, well-studied mathematical objects:

| Theme | Description |
|-------|-------------|
| **Energy functionals** | Concrete \(E(\sigma)\) built from zeros via the Hadamard product; tests of “minimum only at \(\sigma=1/2\)” claims |
| **Fake-zero insertion** | What happens to \(E\) and \(E''\) when a hypothetical off-line zero is added |
| **Backward heat flow** | Gap-lifetime lemma and empirical collision times of zeros under the de Bruijn–Newman flow (\(t\le 0\)) |
| **Forward heat flow** | Winding-number zero counting for \(H_t\) at \(t>0\) (argument principle) |
| **Riemann–Siegel** | Practical \(Z(t)\) evaluator and simple zero surveys |
| **Form factors / BK** | Comparisons with the Bogomolny–Keating formula and function-field analogues (documented) |
| **Scaling** | How closest-pair lifetimes shrink with height (GUE heuristics) |

A detailed self-assessment of an earlier speculative framework is included in `docs/assessment.md`.

---

## Installation

```bash
git clone https://github.com/Satoshi88818/RH-Engine.git   # or unpack the zip
cd RH-Engine

python -m venv .venv
source .venv/bin/activate          # Windows: .venv\Scripts\activate

pip install -e .                   # recommended
# or
pip install -r requirements.txt
```

**Requirements:** Python ≥ 3.9, `numpy`, `scipy`, `mpmath`, `matplotlib`.

---

## Quick start

### 1. Candidate energy functional

```bash
python -m rh_engine.energy 300
# or
python scripts/run_energy.py 300
```

Computes \(E(\sigma)\) from the first 300 zeros and saves `rh_energy.png`. The minimum sits at \(\sigma=1/2\) because the input zeros already lie on the critical line.

### 2. Fake off-line zero test

```bash
python -m rh_engine.fake_zero 300
# or
python scripts/run_fake_zero.py 300
```

Inserts a hypothetical off-line quadruple and measures the resulting kinks in \(E''\). Shows that the minimum of this particular \(E\) remains at \(1/2\) and that the change in \(E\) is finite.

### 3. Backward-flow collision harness

```bash
python -m rh_engine.harness 100,1000,10000 3 results.json
# or
python scripts/run_harness.py 100,1000,10000 3 results.json
```

Tests the gap-lifetime lemma (a real pair stays apart for at least \(r^2/8\)) on a few windows. Full runs (six windows per height up to \(10^6\)) take several minutes on one core.

### 4. Forward heat flow (winding numbers)

```bash
python -m rh_engine.forward 60,300 1 0.02,0.05,0.1 forward_results.json 16
# or
python scripts/run_forward.py 60,300 1 0.02,0.05,0.1 forward_results.json 16
```

Counts real and non-real zeros of \(H_t\) via the argument principle. Includes a positive-control planted pair. Runtime grows quickly with height (mpmath \(\xi\) evaluations).

### 5. Riemann–Siegel survey (basic)

```python
from rh_engine.riemann_siegel import survey_zeros, Z
roots = survey_zeros(14, 50, step=0.05)
print(roots[:5])
```

---

## Package layout

```
rh-engine/
├── pyproject.toml
├── requirements.txt
├── LICENSE                          # MIT
├── README.md                        # this file
├── .gitignore
│
├── src/rh_engine/
│   ├── __init__.py
│   ├── energy.py                    # E(σ) from Hadamard product
│   ├── fake_zero.py                 # off-line zero insertion test
│   ├── harness.py                   # gap-lifetime / collision harness
│   ├── forward.py                   # forward flow + winding numbers
│   ├── riemann_siegel.py            # Z(t) evaluator + simple survey
│   └── py.typed
│
├── scripts/                         # thin CLI wrappers
│   ├── run_energy.py
│   ├── run_fake_zero.py
│   ├── run_harness.py
│   └── run_forward.py
│
├── docs/                            # theoretical notes (cleaned)
│   ├── README.md
│   ├── assessment.md                # critical review of earlier approach
│   ├── energy.md
│   ├── harness.md
│   ├── lifetime_scaling.md
│   ├── function_field.md
│   └── forward_flow.md
│
├── tests/
│   └── test_energy.py
│
├── data/                            # caches, JSON results (gitignored)
├── figures/                         # generated plots
└── notebooks/
```

---

## Documentation map

| Document | Content |
|----------|---------|
| [docs/assessment.md](docs/assessment.md) | Honest review of the earlier “Dynamic Spectral Magnetic Field” variational framework; what holds, what does not |
| [docs/energy.md](docs/energy.md) | Definition and interpretation of the candidate energy \(E(\sigma)\) |
| [docs/harness.md](docs/harness.md) | Gap-lifetime lemma and empirical collision-time results |
| [docs/lifetime_scaling.md](docs/lifetime_scaling.md) | Scaling of closest-pair lifetimes with height (GUE) |
| [docs/function_field.md](docs/function_field.md) | Function-field analogues and Bogomolny–Keating formula |
| [docs/forward_flow.md](docs/forward_flow.md) | Forward heat flow and winding-number methodology |

Additional experimental notes from the original repository (Li criterion, Dirichlet L-functions, primes-and-repulsion statistics, full Riemann–Siegel tables, etc.) can be restored as further modules if needed.

---

## Important caveats

- **No proof.** The code and notes never claim to prove the Riemann Hypothesis.
- **Floating-point only.** Results are not interval-arithmetic verified.
- **Known zeros as input.** Most experiments take zeros already known to lie on the critical line.
- **Heuristic extrapolations.** Scaling laws that reach heights beyond published tables assume GUE small-gap statistics continue to hold.
- **Performance.** Forward-flow and high-precision zero computations are intentionally slow (high-precision `mpmath`); they are demonstration tools, not production scanners.

---

## Scientific context (brief)

The de Bruijn–Newman constant \(\Lambda\) is defined so that RH \(\Leftrightarrow \Lambda\le 0\). It is known that \(\Lambda\ge 0\) (Rodgers–Tao) and \(\Lambda\le 0.22\) (Polymath15). The heat-flow picture of the zeros of \(\Xi\) (or of \(H_t\)) supplies a quantitative language in which “how barely true” RH would be can be discussed: closest-pair lifetimes under backward flow shrink with height, consistent with \(\Lambda=0\) if RH holds.

Littlewood’s lemma supplies an equivalent formulation of RH in terms of the support of a certain measure derived from \(\log|\xi|\). The energy experiments in this package make that statement concrete for finite windows.

None of the numerical work here reaches the analytic difficulty of controlling the flow at \(t=0\) with arithmetic input—the hard direction of any barrier argument.

---

## Licence

MIT — see `LICENSE`.

## Provenance

Scientific content originates from [Satoshi88818/RH-Engine](https://github.com/Satoshi88818/RH-Engine) (November 2026). This reorganisation extracts the embedded scripts, imposes a conventional package layout, adds packaging metadata, tests, and a coherent README while preserving the original scientific claims and caveats.
