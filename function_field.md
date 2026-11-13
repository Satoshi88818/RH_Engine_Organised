# Function-field analogue and the full Bogomolny–Keating formula

## Summary

1. **Function fields (where RH is a theorem) show the same ramp-then-saturation form factor, and the saturation point is exactly the number of zeros.** For hyperelliptic curves y² = f(x) of genus g over F_q, the family-averaged form factor K = E[t_n²]/2g ramps up and equals 1 from n = 2g onward, as the symplectic-matrix theory predicts (measured for q = 17, g = 4: 1.02 ± 0.02 at τ = n/2g = 1, and 0.98–1.02 beyond; at q = 3 the values at τ = 1 are about 1.07 because of the small-q excess below). The finite-field correction to the matrix prediction scales like 1/q (q·deviation ≈ 1.6–2.0 for q = 3 … 17). All point counts pass Weil's check |α| = √q.
2. **The full Bogomolny–Keating (BK) two-point function, evaluated numerically from the primes with no free parameters, matches the data from γ ≈ 10³ to 10²¹.** The rms z-score of close-pair counts over all 12 datasets and 4 thresholds is 0.80, against 9.8 for plain random-matrix theory and 2.8 for the first-order correction.
3. **It fixes the failure of the first-order correction at low height.** For ln(γ/2π) < 13 (N_eff < 3) the rms z-scores are 15.5 (random matrix), 4.7 (first-order), **0.7 (full BK)**.
4. **It also resolves my earlier line-weight puzzle.** The excess lines in K(τ) near τ = 1 (at 1 + ln q / L) are reproduced by the full formula (measured/BK ratio 0.97 on average for q = 2 and 6). My earlier comparison with the simple weight g(q)/L² was the wrong comparison, because it ignored the smooth background that the same formula puts around the lines.
5. **The BK formula itself predicts the saturation K = 1 for τ > 1** (numeric K(1.3) = 1.0007, K(1.5) = 0.9994 at γ ≈ 10⁹). In that formula it comes from the off-diagonal term, which is built from prime-pair correlations (Hardy–Littlewood). So the "why does it stop at τ = 1" question has a conjectural arithmetic answer, which agrees with the data; it is not a proof.

None of this bears on whether RH is true.

## 1. Function-field analogue

**Setup.** For y² = f(x), f squarefree and monic of degree 2g + 1 over F_q (q prime), the L-function is a polynomial of degree 2g whose roots α_i = √q e^{iθ_i} play the role of the zeros. I counted points by brute force over the fields F_{q^n}, n ≤ g (log tables), got the L-polynomial from Newton's identities and the functional equation, and computed the normalised traces t_n = Σ_i e^{inθ_i}. Families: g = 3 over q = 3 … 53, g = 4 over q = 3 … 17, plus g = 6, 8, 10 over q = 5, 3, 3. Between 4,000 and 30,000 random curves each.

**Checks.** For every curve tested, |α_i| = √q to 10⁻⁸ or better (Weil's theorem, which would fail if a count were wrong); the traces are integers after rescaling; E[t_n] = −1 for even n ≤ 2g and 0 for odd n, as expected (q = 17, g = 4: n = 2, 4, 6, 8 give -0.97, -0.93, -1.00, -1.04).

**The matrix prediction.** As q → ∞ the Frobenius classes of this family equidistribute in USp(2g) (Katz–Sarnak), for which E[t_n²] = n + η_n for n ≤ g, n − 1 + η_n for g < n ≤ 2g, and 2g for n > 2g (η_n = 1 for even n, else 0). So K(τ) = E[t_n²]/2g with τ = n/2g ramps to 1 at n = 2g and stays there. The saturation is a theorem here: the traces for n > 2g are determined by the first 2g.

**Deviation from USp(2g) versus field size** (rms of E[t_n²] − prediction over n = 1 … 2g):

| g | q | curves | rms deviation | statistical noise | q × deviation |
|---|---|---|---|---|---|
| 3 | 3 | 30,000 | 0.650 | 0.027 | 1.95 |
| 3 | 5 | 30,000 | 0.355 | 0.027 | 1.78 |
| 3 | 7 | 30,000 | 0.252 | 0.028 | 1.76 |
| 3 | 11 | 30,000 | 0.185 | 0.028 | 2.04 |
| 3 | 17 | 15,000 | 0.118 | 0.039 | 2.00 |
| 3 | 29 | 8,000 | 0.087 | 0.054 | 2.53 |
| 3 | 53 | 4,000 | 0.095 | 0.077 | 5.02 |
| 4 | 3 | 30,000 | 0.627 | 0.038 | 1.88 |
| 4 | 5 | 30,000 | 0.364 | 0.037 | 1.82 |
| 4 | 7 | 30,000 | 0.224 | 0.038 | 1.57 |
| 4 | 11 | 15,000 | 0.153 | 0.053 | 1.68 |
| 4 | 13 | 10,000 | 0.136 | 0.064 | 1.76 |
| 4 | 17 | 6,000 | 0.107 | 0.083 | 1.81 |

The deviation falls like 1/q (the last column is roughly constant at 1.6–2.0 until it hits the noise floor).

**Form factor for g = 4, q = 17:**

| n | τ = n/2g | measured K | USp |
|---|---|---|---|
| 1 | 0.125 | 0.121 ± 0.002 | 0.125 |
| 2 | 0.250 | 0.350 ± 0.006 | 0.375 |
| 3 | 0.375 | 0.382 ± 0.007 | 0.375 |
| 4 | 0.500 | 0.608 ± 0.011 | 0.625 |
| 5 | 0.625 | 0.499 ± 0.009 | 0.500 |
| 6 | 0.750 | 0.749 ± 0.012 | 0.750 |
| 7 | 0.875 | 0.742 ± 0.013 | 0.750 |
| 8 | 1.000 | 1.020 ± 0.017 | 1.000 |
| 9 | 1.125 | 1.022 ± 0.017 | 1.000 |
| 10 | 1.250 | 0.977 ± 0.016 | 1.000 |
| 12 | 1.500 | 0.994 ± 0.016 | 1.000 |
| 16 | 2.000 | 0.998 ± 0.016 | 1.000 |

**Large genus at small q** (K at τ = 0.25, 0.5, 0.75, 1, 1.5, 2):

| family | 0.25 | 0.5 | 0.75 | 1.0 | 1.5 | 2.0 |
|---|---|---|---|---|---|---|
| q=3, g=8 | 0.270 ± 0.003 | 0.587 ± 0.006 | 0.753 ± 0.007 | 1.068 ± 0.010 | 0.993 ± 0.009 | 1.022 ± 0.010 |
| q=3, g=10 | 0.244 ± 0.005 | 0.571 ± 0.011 | 0.712 ± 0.014 | 1.068 ± 0.020 | 1.008 ± 0.020 | 1.016 ± 0.020 |
| q=5, g=6 | 0.242 ± 0.002 | 0.593 ± 0.006 | 0.673 ± 0.006 | 1.050 ± 0.010 | 1.004 ± 0.009 | 1.009 ± 0.010 |

Just above the Heisenberg point there is a small excess that decays: E[t_n²] − 2g for n − 2g = 1, 2, 3, 4 is +0.95, +0.37, +0.10, +0.05 for q = 3, g = 8 (noise about ±0.15) and +0.38, -0.01, -0.11, -0.03 for q = 5, g = 6 (noise about ±0.11). This is qualitatively like the excess above τ = 1 for the zeta zeros, but I cannot map it quantitatively.

**Dictionary.** A prime power of the zeta function with log n = u corresponds to a prime polynomial of degree d with u = d ln q. The point τ = 1, u = ln(T/2π), corresponds to d = 2g, i.e. u = 2g ln q: **"Heisenberg time" is the degree at which the number of independent traces exceeds the number of zeros.** In function fields that forces K = 1 exactly. For the zeta zeros the same picture is the natural heuristic, but there it is not a proof.

## 2. The full Bogomolny–Keating formula, numerically

**Formula.** r₂(ε) = ρ² + r_diag + r_off with

> r_diag = −(1/4π²) ∂²_ε [ ln|ζ(1+iε)|² + 2 Σ_p Σ_{r≥2} (1−r)/(r² p^r) cos(εr ln p) ],
> r_off = (1/4π²) |ζ(1+iε)|² e^{i2πρε} ∏_p [1 − (1−p^{iε})²/(p−1)²] + c.c.

with ρ = ln(γ/2π)/2π and R₂(s) = r₂(s/ρ)/ρ². I evaluated it with mpmath (ζ and its first two derivatives at 1 + iε), the Euler product over primes ≤ 2×10⁵ and the prime-power sum over p ≤ 10⁴, r ≤ 14, and Fourier-transformed it (Hann window, S = 400 spacings, resolution about 0.005 in τ) to get K(τ).

**Validation.** (i) The Fourier-transform code gives the exact ramp K = τ for the random-matrix R₂ = 1 − sinc². (ii) The 1/ε² poles of r_diag and r_off cancel, so R₂ → 0 like s², with small-s strength 0.94–0.97 of the random-matrix value at heights 10⁸–10¹⁰.

**Form factor, measured versus BK** (χ² per bin, BK / plain min(τ,1); the number of bins is 32, 18, 14 and the statistical errors are about 0.03–0.1):

| height | ln(γ/2π) | 0.3 < τ < 0.95 | 0.95 < τ < 1.3 | 1.3 < τ < 1.58 | correlation of fine structure, τ = 1–1.3 |
|---|---|---|---|---|---|
| 1e8 | 16.6 | 0.24 / 0.28 | 0.81 / 1.93 | 1.31 / 1.41 | 0.64 |
| 1e9 | 18.9 | 0.09 / 0.08 | 0.19 / 0.83 | 0.80 / 0.81 | 0.66 |
| 1e10 | 21.2 | 0.38 / 0.38 | 0.20 / 0.33 | 1.06 / 1.05 | 0.37 |

BK is as good as the ramp where the ramp is already good, and better in the region 0.95 < τ < 1.3 where the lines live (e.g. 0.81 against 1.93 at γ = 10⁸; figure, lower left).

**Line weights in a ±0.009 window around τ = 1 + ln q/L** (excess over K = 1):

| q | height | measured | full BK | simple analytic g(q)/L² |
|---|---|---|---|---|
| 2 | 1e8 | 0.00425 | 0.00392 | 0.00728 |
| 2 | 1e9 | 0.00269 | 0.00319 | 0.00560 |
| 2 | 1e10 | 0.00296 | 0.00238 | 0.00445 |
| 6 | 1e8 | 0.00236 | 0.00335 | 0.00546 |
| 6 | 1e9 | 0.00247 | 0.00275 | 0.00420 |
| 6 | 1e10 | 0.00211 | 0.00201 | 0.00334 |
| 3 | 1e8 | 0.00051 | 0.00060 | 0.00273 |
| 3 | 1e9 | 0.00039 | 0.00025 | 0.00210 |
| 3 | 1e10 | -0.00015 | 0.00034 | 0.00167 |
| 10 | 1e8 | 0.00070 | 0.00049 | 0.00227 |
| 10 | 1e9 | 0.00004 | 0.00008 | 0.00175 |
| 10 | 1e10 | 0.00017 | 0.00024 | 0.00139 |

For the two strong lines (q = 2, 6) the measured/BK ratio is 1.08, 0.71, 0.84, 0.90, 1.24, 1.05 (mean 0.97); the simple analytic weight is 1.6–2 times too large for q = 2, 6 and 5–20 times too large for the weak lines, so the shortfall I reported last time came from comparing with the wrong quantity. The weak lines (q = 3, 10) are at the noise level (±0.0002).

**Close-pair counts, all heights** (x = 0.5 shown; observed, random-matrix expected, first-order BBLM expected, full BK expected; z-scores in that order; the other thresholds are in `bk_lowL.json`):

| dataset | ln(γ/2π) | observed | RM | BBLM | full BK | z (RM / BBLM / BK) |
|---|---|---|---|---|---|---|
| zeros6 [1e+03,1e+04 | 6.7 | 850 | 1,073 | 680 | 844 | -6.8 / +6.5 / +0.2 |
| zeros6 [1e+04,1e+05 | 9.0 | 12,504 | 14,463 | 11,626 | 12,420 | -16.3 / +8.1 / +0.8 |
| zeros6 [1e+05,1e+06 | 11.3 | 164,471 | 181,915 | 159,804 | 164,203 | -40.9 / +11.7 / +0.7 |
| zeros6 [1e+06,1.1e+06 | 12.0 | 26,401 | 28,705 | 25,694 | 26,230 | -13.6 / +4.4 / +1.1 |
| zeros3 (2.7e11 | 24.5 | 1,112 | 1,130 | 1,103 | 1,104 | -0.5 / +0.3 / +0.2 |
| zeros4 (1.4e20 | 44.6 | 1,173 | 1,130 | 1,122 | 1,105 | +1.3 / +1.5 / +2.0 |
| zeros5 (1.4e21 | 46.8 | 1,124 | 1,130 | 1,123 | 1,105 | -0.2 / +0.0 / +0.6 |
| mine 1e7 | 14.3 | 52,857 | 56,533 | 52,374 | 52,922 | -15.5 / +2.1 / -0.3 |
| mine 1e8 | 16.6 | 53,948 | 56,531 | 53,472 | 53,780 | -10.9 / +2.1 / +0.7 |
| mine 1e9 | 18.9 | 54,422 | 56,528 | 54,187 | 54,374 | -8.9 / +1.0 / +0.2 |
| mine 1e10 | 21.2 | 32,844 | 33,930 | 32,821 | 32,892 | -5.9 / +0.1 / -0.3 |
| mine 1e11 | 23.5 | 16,492 | 16,970 | 16,521 | 16,544 | -3.7 / -0.2 / -0.4 |

Summary of rms z-scores over x = 0.1, 0.2, 0.3, 0.5: low heights (L < 13): 15.5 / 4.7 / 0.70; L ≥ 13: 5.05 / 1.03 / 0.85. For the nearest-neighbour density I used p(s) ≈ R₂(s), which is accurate to a negligible correction for s < 0.5, normalised so that the random-matrix limit is the exact determinant value.

![Function-field and BK results](rh_ff_bk.png)

## 3. What this adds to the mechanism

| step | status after this work |
|---|---|
| zeros ↔ primes | exact (explicit formula), measured (`README_primes_and_repulsion.md`) |
| ramp τ < 1 | primes; verified |
| saturation at τ = 1 | **theorem in function fields** (finite dimension 2g); **reproduced by the BK formula** for the zeta zeros, where it comes from the off-diagonal prime-pair term (conjectural) |
| finite-height weakening of repulsion | **quantitatively reproduced by the full BK formula over 10³ … 10²¹ with no free parameters** |
| lines near τ = 1 | reproduced by the full BK formula |

So the picture "primes make the ramp; prime-pair correlations stop it at the Heisenberg time and shape the finite-height corrections" is now supported by a parameter-free numerical test across 18 orders of magnitude in height. It remains conditional on the Hardy–Littlewood prime-pair conjecture and on the further assumptions behind the BK formula.

## 4. Corrections to the previous report

`README_primes_and_repulsion.md`, §6 reported measured line weights about 25–30% below the Hardy–Littlewood-based prediction and left the cause open. The cause is that I compared with the simple analytic weight and not with the full formula; with the full formula the measured weights agree within noise (§2 above). §5 there showed the first-order correction failing for N_eff < 3; the full formula removes that failure. I added a note to that file.

## 5. Limitations

- Function fields: odd-degree hyperelliptic family over prime fields only; q·deviation ≈ constant is an empirical scaling from q = 3 … 17; no attempt to fit constants.
- The numeric BK function uses truncated Euler products and prime sums (truncation errors that I estimate at about 10⁻⁶ in the product and 10⁻³ in the prime-power sum), a window of 400 spacings, and the approximation p(s) ≈ R₂(s) for s < 0.5.
- The z-scores treat counts as Poisson. Neighbouring gaps are correlated and the four thresholds are nested, so the rms z-scores are optimistic and not independent.
- The BK formula depends on the Hardy–Littlewood conjecture and on further assumptions in its derivation; this test supports it numerically, it does not prove it.
- The zeros are from my computation or from a GitHub mirror of Odlyzko's tables (sanitised, numbers unchanged per its author).
- Nothing here bears on RH.

## 6. Next steps

- Test the BK prediction for higher correlations (three-point), which the same prime-pair structure predicts.
- Push the function-field analogue to the regime that corresponds to the zeta zeros at fixed q and growing g, with a quantitative theory of the small-q excess above n = 2g.
- Examine the weak lines (q = 3, 10) with larger zero samples.

## Files

`ff_hyperelliptic.py`, `run_ff.py`, `run_rest.py`, `ff_analyze.py` (function fields); `bk_pred.py`, `bk_compare.py`, `bk_lowL.py` (BK numerics); `build_ff_bk.py`, `write_ffbk_readme.py`; results `ff_results.json`, `ff_analysis.json`, `bk_compare.json`, `bk_lowL.json`, `bk_R2_table.npz`, `bk_L*.npz`; figure `rh_ff_bk.png`. They need `primes_2e8.npy`, the position files from `run_all.sh`, and `odlyzko/zeros{3,4,5,6}`.

```bash
pip install mpmath numpy scipy matplotlib
python run_ff.py && python ff_analyze.py              # function fields
python bk_pred.py 18.89 400 0.03 bk_L18.89.npz        # K(tau) at gamma ~ 1e9 (about 6 min)
python bk_compare.py && python bk_lowL.py             # comparison and close-pair counts
```
