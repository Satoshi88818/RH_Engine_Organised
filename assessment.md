Assessing the “Reimann-Hypothesis-” repository: what holds, what doesn’t, and what the numerics show
Consolidated report · 1 November 2026
1. Bottom line
• The repository does not prove the Riemann Hypothesis. Its central object, the energy E(σ), is fixed entirely by the gamma factor in the functional equation and contains no information about the zeros, so it cannot distinguish RH from its negation.
• Other steps rest on undefined terms (“endomorphic wave”, “magnetic gravity tensor”) or on unconditional theorems that are equally consistent with RH being false.
• What is sound in its picture can be made precise with existing mathematics: the heat flow of the zeros of Ξ (de Bruijn–Newman), Littlewood’s lemma, and a simple gap-lifetime lemma. I built and tested that framework numerically on millions of computed and published zeros, up to γ ≈ 2.7×10¹¹, plus published zeros near γ ≈ 1.4×10²⁰ and 1.4×10²¹.
• The numerics agree with known theory everywhere and found nothing off the critical line, as expected, since RH is verified far beyond these heights. Their value is quantitative: they show how “barely true” RH would be. The first backward-collision time of the zeros shrinks roughly like x^(−2/3) with height (a heuristic, supported by tests on windows of up to 3×10⁴ zeros).
• None of this bears on whether RH is true.
2. What was examined, and its limits
• The repository README (fetched from GitHub) and the file RHTheorem.txt (uploaded by you). I read the opening part of RHTheorem.txt and searched it for the definition of the energy. I did not read the whole file, and I did not read RHProposal or any other file in the repository.
• Everything numerical is my own code, run on one CPU core, in floating-point arithmetic (not interval arithmetic). No result is a rigorous verification.
• Odlyzko’s published zero tables were obtained through a GitHub mirror (nightbay/odzeros), whose author says the files were whitespace-sanitised with the numbers unchanged. I did not compare them byte-for-byte with Odlyzko’s originals.
3. Problems with the repository’s argument
1. The energy is vacuous. RHTheorem.txt defines E(σ) = ∫ [log|ζ(σ+it)| − log|ζ(1−σ+it)|]² dt. By ζ(s) = χ(s)ζ(1−s) and conjugation symmetry, the integrand is exactly (log|χ(σ+it)|)², which depends only on the gamma factor (I checked this numerically at three points to 10 digits). Since |χ| = 1 on the critical line and log|χ| ≈ (½−σ)·ln(t/2π) elsewhere, E is 0 at σ = ½ and +∞ for every other σ, whatever the zeros do. Inserting an off-line zero changes nothing in this E, so the claimed “infinite action” and “ground state” are properties of χ, not of ζ’s zeros.
2. RH itself is stated precisely in Stage 1 of RHTheorem.txt (all non-trivial zeros have real part ½). The README’s final statement says only “the zeros of ζ(s)”, which is false as written because of the trivial zeros.
3. Undefined terms and reinterpretation. “Endomorphic”, “magnetic gravity tensor”, “load balancer” and similar terms are never defined in a way that supports a deduction. Recasting RH as a ground-state, stability or conservation statement is asserted, not shown to be equivalent.
4. Unconditional theorems presented as obstructions. The claim that an off-line zero would force the AM-GM inequality, the zero-density formula or Jensen’s formula to fail is wrong: these are theorems that hold whether or not RH is true.
5. The “infinite action” conclusion is unsupported (see §5: with a concrete energy, an off-line zero changes E by a tiny finite amount and the minimum stays at ½).
6. The Yang–Mills analogy is odd, since confinement and the mass gap are themselves open Millennium Problems.
7. Impression, not a finding: the long prose, the Stage 1–21 structure and the absence of any computation resemble an AI-assisted brainstorming session more than a verified argument.
4. A precise framework
Normalisation (Polymath15 / Rodgers–Tao). H_t(z) = ∫ e^{tu²} Φ(u) cos(zu) du, with H₀(z) = ξ(½ + iz/2)/8. So zeros sit at z = 2γ and s = ½ + iz/2. Then ∂_t H = −∂_zz H, and simple zeros move by
ż_k = 2 Σ_{j≠k} 1/(z_k − z_j).
This is the law Newman–Wu cite from Csordas–Smith–Varga 1994. The de Bruijn–Newman constant Λ is the smallest t for which all zeros are real. Known: RH ⇔ Λ ≤ 0; Λ ≥ 0 (Rodgers–Tao 2018); Λ ≤ 0.22 (Polymath15); real zeros stay real forward in t.
Littlewood’s lemma. The second derivative in σ of the window-averaged log|ξ| is, asymptotically, a measure supported on the real parts of the zeros. So RH ⇔ that measure sits only at σ = ½. This is a correct reformulation of RH, not a proof.
Dictionary for the repo’s language (my proposal; the repo does not define its terms):
Repo language
Precise object
zero as a wave
term x^ρ/ρ of the explicit formula, ρ = ½ + iz/2
amplitude envelope
x^(σ−½), with offset δ = Re ρ − ½ = −Im z/2
phase velocity
frequency γ in log x; or the drift of the phase under the flow, ȧ = Re ż
repulsive force
Re ż: real zeros repel (Coulomb)
restoring force
Im ż: a conjugate pair attracts, ḃ = −1/b
Littlewood cusp
kink of the averaged log
endomorphic
not mappable
An isolated complex pair reaches the axis at t* = b₀²/2. Real background zeros add attraction of roughly ½ ln(γ/2π), so collisions come sooner.
Gap-lifetime lemma. For consecutive real zeros with gap r, dr/dt = 4/r + τr with τ = −2 Σ 1/((z_{j+1}−z_m)(z_j−z_m)) ≤ 0, so d(r²)/dt ≤ 8. Hence going backward in time a real pair stays apart for at least r²/8 (in γ units, gap²/2). This is the leading term of the Csordas–Smith–Varga theorem behind all Lehmer-pair lower bounds on Λ; their correction involves the tidal sum, which I estimate at about 10⁻⁶ relative for a Odlyzko-2000-type pair. The lemma is exact for finitely many real zeros; I assume it extends to the infinite system.
Consequences.
• RH ⇔ no collision of zeros at any t > 0.
• A uniform lower bound on gaps at t = 0 is heuristically impossible: if RH holds then Λ = 0, so collisions occur arbitrarily close to t = 0 from below, and the available lifetime tends to 0 with height.
• A workable no-collision argument must be height-dependent, and at t = 0 it needs arithmetic input (the analytic approximations that control zeros for t > 0 do not exclude off-line zeros at t = 0).
5. Experiments and results
Experiment
Setup
Result
Candidate energy from zeros
mean of log
ξ(σ+it)/ξ(0)
Fake off-line zero in that energy
quadruple ½ ± δ ± iγ₀ added, δ 0.02–0.40
kink of size 2π/L = 0.0261 at ½ ± δ, matching theory to 4 digits; minimum stays at ½; change in E is finite and tiny; a zero outside the window is invisible
Flow with a fake zero
300 real zeros plus a quadruple
collision times 0.92, 0.86, 0.72, 0.61 × the isolated estimate for δ = 0.05, 0.10, 0.20, 0.30
Backward flow (exact smoothing), mpmath
30 windows, 655 zeros, γ 10² to 10⁶
s_c / (r²/8) between 0.98 and 1.15 (median 1.08, 25 of 26 ≥ 1); mean gap × ln(γ/2π) = 12.5 ≈ 4π
Forward flow (t > 0)
ξ blurred across σ, Gauss–Hermite, winding-number zero count; 8 windows, 36 zeros, γ up to 3000, t = 0.02–0.2
planted non-real pair detected in 32 of 32 runs; none found otherwise; zero positions match the ODE flow to 2×10⁻⁴ (t = 0.05) after correcting for the truncated tail
Riemann–Siegel evaluator
C₀–C₂ Chebyshev tables, long-double phases, fast recurrence
error vs mpmath about 10⁻⁶ at γ = 2×10³, 6×10⁻¹² at 10⁵, 5×10⁻⁹ at 10⁹, 10⁻⁶ at 10¹¹
Survey
2,000,405 zeros (10³–1.13×10⁶) plus about 1.95 million at 10⁷–10¹¹
count vs Riemann–von Mangoldt within ±1.1 in every run
Odlyzko zeros6
2,000,403 table zeros, 1000 < γ < 1,132,490.66
all matched one-to-one; the 25 closest pairs are the same set and their gaps agree to 3.0×10⁻⁹; closest pair in the table: γ = 663318.508, gap 0.002958654; Lehmer pair 7005.0629 has gap 0.037698
Odlyzko zeros3
10,000 zeros near γ = 2.68×10¹¹
all matched; median position difference 3.4×10⁻⁶; the 12 closest pairs are the same set, gaps agree to 5×10⁻⁷
Lemma on extreme pairs
exact backward flow around the closest pairs
21 pairs up to γ = 10⁹: 1.00008–1.0059 (median 1.0002); 9 pairs at 10¹⁰–2.7×10¹¹: 1.0001–1.0010; 150 typical windows: median 1.05
Window test
windows of N₀ = 100 … 30,000 consecutive zeros, 10 datasets, 58,467 windows
smallest gap follows ν_min ∝ N₀^(−1/3) (mean fitted exponent −0.335, theory −0.333); pooled coverage of the 5% / 1% quantiles 0.955 / 0.991 (expected 0.95 / 0.99)
Published zeros near 10²¹ and 10²²
10,000 zeros each, ln(γ/2π) ≈ 45–47
small-gap counts match the random-matrix law within about one standard deviation; window statistics pass (KS p 0.46–0.81, coverage 0.97–1.00; only N₀ = 100 and 300 testable)
Notes on the window test. The small-gap law used is F(x) ≈ (π²/9)x³ − (2π⁴/225)x⁵ (GUE). Windows of different sizes reuse the same zeros, so the pooled counts are not independent: the coverage figures show agreement, not statistical significance. Below about γ = 10⁶ close pairs are rarer than the law predicts (observed/expected about 0.5 at 10³, 0.8–0.85 at 10⁵–10⁶, about 0.9 from 10⁷), so the test rejects it there. The deficit fades with height, which is the known slow convergence of zeta zeros to GUE.
6. Scaling of the first-collision time
For the closest pair in a stretch of N zeros, lifetime ≈ gap²/2, and under the GUE small-gap law its median is
τ_med ≈ (2π²/L²) · (9 ln 2 / (π² N))^(2/3), with L = ln(γ/2π).
For all N ≈ x(L + 2 ln 2 − 1)/2π zeros in [x, 2x] this gives τ ∝ x^(−2/3) L^(−8/3). Typical windows (not the closest pair) scale like 1/L² instead (median s_c × ln²γ is about 1.1–2.1 between 10⁵ and 10⁹).
x
zeros in [x, 2x]
τ, 5% quantile
τ, median
10⁹
3.1×10⁹
3.4×10⁻⁹
1.9×10⁻⁸
10¹¹
3.8×10¹¹
8.9×10⁻¹¹
5.0×10⁻¹⁰
10¹³
4.5×10¹³
2.6×10⁻¹²
1.4×10⁻¹¹
10¹⁵
5.3×10¹⁵
7.9×10⁻¹⁴
4.5×10⁻¹³
10²⁰
7.1×10²⁰
1.6×10⁻¹⁷
9.3×10⁻¹⁷
10²²
7.8×10²²
5.9×10⁻¹⁹
3.3×10⁻¹⁸
Status: the N₀^(−2/3) scaling is tested for N₀ from 100 to 30,000. The step to the 10¹³–10²² zeros in the table is an extrapolation that assumes the GUE small-gap law keeps holding. No published table can test it directly.
7. Consistency with the literature
• Lower bounds on Λ from Lehmer pairs: −4.379×10⁻⁶ (Csordas–Smith–Varga 1994), −5.895×10⁻⁹ (Csordas–Odlyzko–Smith–Varga 1993), −2.7×10⁻⁹ (Odlyzko 2000, pair near zero number 10²⁰), −1.15×10⁻¹¹ (Saouter–Gourdon–Demichel 2011). Since 2018, Λ ≥ 0 (Rodgers–Tao), proved by a different method, supersedes all of these.
• By the lemma, a bound B implies a pair with gap √(2B) in γ. Under GUE such a pair is expected about once per 7.5×10⁹ zeros for the 2000 bound and about once per 10¹³ zeros for the 2011 bound. My extrapolation puts the median closest pair among about 4.5×10¹³ zeros at 1.4×10⁻¹¹. This is an order-of-magnitude consistency check; I could not find the exact height, gap or scanning effort of the 2011 pair.
8. What this does not show
• Nothing here says any zero is, or is not, off the critical line.
• The framework cannot yet reach t = 0: the hard direction needs arithmetic information at every height.
• The x^(−2/3) law is an extrapolation of 7–11 orders of magnitude in N.
• The lemma’s extension from finite to infinite systems is assumed, and all numerics are floating-point.
9. Errors made and corrected during this work
1. I first said the repo never defines an energy. It does (§3); my fake-zero test used a different energy, so it does not test the repo’s.
2. I cited Hardy’s convexity theorem for a truncated mean of log|ξ|, which it does not cover. On a finite window the energy has a small negative curvature (about −0.01).
3. I called a uniform lower bound on gaps at t = 0 “the open core” of the problem; it is heuristically impossible (§4).
4. I conjectured a window time t(x) ∝ 1/log²x. That holds only for typical windows; for the closest pair it falls much faster (§6).
5. I called Λ > −1.15×10⁻¹¹ “the best published lower bound”. It is the best Lehmer-pair bound; Λ ≥ 0 has superseded it.
6. Survey bugs, all found by cross-checking: a Newton root finder that clipped for pairs straddling a grid point (fake pairs with gaps near 10⁻⁵); Gibbs ripples from zero-padded window edges (fake early collisions); and close pairs straddling window boundaries (a 7% gap error against the Odlyzko table). The first survey tables in an earlier README were therefore superseded.
7. An early harness run flagged a collision at 0.57 × the bound; it was a zero drifting across a fixed window edge.
10. Open questions and next steps
• Prove the gap-lifetime lemma for the infinite system (Laguerre–Pólya class), including far-field sums.
• Find the arithmetic input at t = 0 that a barrier argument needs; this is the real difficulty.
• Test the random-matrix small-gap law at the 10¹³-zero scale (the Odlyzko–Schönhage algorithm would be needed to compute zeros that high).
• Extend the forward-flow tool with a fast approximation of ξ so the winding-number test reaches large heights.
• Not explored: the Li-criterion route, and finite-field analogues.
• Read the rest of the repository (RHProposal and the other files).
11. Files
• Reports (Markdown): README_rh_energy, README_fake_zero_test, README_phase_velocity_map, README_no_collision_bound, README_harness, README_forward_flow, README_rs_scaling, README_t_of_x_test, README_literature_extreme_heights. Where numbers differ, README_t_of_x_test and README_literature_extreme_heights supersede README_rs_scaling.
• Code: rh_energy.py, rh_fake_zero_test.py, rh_flow_sketch.py, rh_nocollision_sketch.py, rh_harness.py, rh_forward.py, rs_tables.py, rs_z.py, rh_survey.py, rh_survey2.py, rh_flow2.py, compare_odlyzko.py, window_test.py, final_stats.py, flow_check_high.py, extreme_heights.py, build_report.py, run_all.sh.
• Data and results: results_final.json, forward_results.json, survey_results.json (earlier survey), flow2_results.json, compare_odlyzko.json, window_test.json, final_stats.json, flow_check_high.json, extreme_heights.json, rs_cheb_coef.npy, and the figures.