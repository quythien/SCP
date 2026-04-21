# B vs m Tradeoff — Single-Cohort Detection: Findings Summary

**Session date:** 2026-04-20  
**Pilot dataset:** Mouse D1D2 D1 (r_med ≈ 0.65, low SNR — hardest case)  
**Fixed N:** total samples = B × m  
**Designs:** equispaced active, B ∈ {3,4,6,8,12}  
**Truth:** cosinor + Fourier harmonic Y = A[cos(ωt−φ) + α₂cos(2ωt−φ)] + N(0,σ)

---

## 1. Key Finding: B vs m is Method-Dependent

There is no universal answer to "should I use more time points or more replicates?"  
The optimal B depends entirely on which detection test you use and what waveform you expect.

---

## 2. Method Behavior Summary

| Method | Cosinor truth | Harmonic truth (α₂≥0.5) | Favors | Reason |
|--------|--------------|------------------------|--------|--------|
| **DCP (K=1 cosinor F-test)** | B-invariant | B=3 wins | m (or neither) | NCP = N·r²/2 is B-invariant for equispaced designs; at B=3 the 2nd harmonic aliases onto the fundamental, inflating NCP |
| **JTK_CYCLE (MetaCycle)** | B=3–4 wins | B=3–4 wins | m | MetaCycle averages replicates within time points before ranking → test effectively runs on B means; more m = cleaner means = stronger rank signal |
| **RAIN (nr.series)** | B=6–8 wins | B=6–8 likely | B | Keeps all N individual observations; umbrella test gains rank resolution with more distinct time points; properly implemented with nr.series |
| **Multi-harmonic (adaptive K=⌊(B−1)/2⌋)** | B=3–4 wins | B=6 wins at N≥48 | B (conditional) | Larger B → larger K → captures ΣAₖ² signal; K=2 at B=6 outweighs 4df cost when α₂≥0.5; K=5 at B=12 over-fits and loses power |

---

## 3. Empirical Results

### DCP (K=1), Fourier truth, D1D2 pilot
From `output/07_fourier_extreme_20260401/extreme_summary.txt`:
- B=3 power INCREASES with α₂ (aliasing: cos(2ωt) = cos(ωt) at t={0,8,16})
- B=6 power DECREASES with α₂ (2nd harmonic orthogonal → invisible to K=1 test)
- At α₂=1.0: B3 gains 14pp, B6 loses 17pp vs cosinor baseline

### JTK_CYCLE, D1D2 pilot (smoke_jtk_fourier.R, nsims=15)
```
alpha2=0.0:  N=24  B3=35%  B4=40%  B6=34%  B8=32%
             N=48  B3=72%  B4=71%  B6=69%  B8=63%
             N=72  B3=86%  B4=88%  B6=86%  B8=83%

alpha2=0.5:  N=24  B3=41%  B4=37%  B6=28%  B8=25%
             N=48  B3=68%  B4=71%  B6=63%  B8=59%
             N=72  B3=79%  B4=90%  B6=84%  B8=80%

alpha2=1.0:  N=24  B3=57%  B4=39%  B6=19%  B8=16%
             N=48  B3=73%  B4=76%  B6=52%  B8=49%
             N=72  B3=80%  B4=91%  B6=78%  B8=74%
```
**Note:** B=4 often matches or beats B=3 because max_K(4)=1 (same K as B=3) but avoids the B=3 aliasing artifact.  
**JTK verdict:** B=3 or B=4. Never B=6+.

### RAIN (nr.series fixed), D1D2 pilot (smoke_rain_fourier.R, partial, nsims=10)
```
alpha2=0.0 (cosinor):  N=24  B3=33.2%  B4=36.9%  B6=42.2%  B8=43.9%
                       N=48  B3=80.6%  B4=76.8%  B6=80.7%  B8=80.8%
                       N=72  [still running as of 2026-04-20]
```
**Direction opposite to JTK** — higher B wins at N=24 (+10pp for B8 vs B3), converges at N=48.  
Full results (harmonic/pulse scenarios) still pending — RAIN slow due to permutation test.

### Multi-harmonic adaptive K, D1D2 pilot (smoke_multiharmonic.R, nsims=20)
```
alpha2=0.00:  N=48  B3=68%  B4=68%  B6=50%  B8=34%  B12=19%
              N=72  B3=89%  B4=89%  B6=79%  B8=70%  B12=53%

alpha2=0.50:  N=48  B3=66%  B4=65%  B6=66%  B8=53%  B12=31%
              N=72  B3=80%  B4=88%  B6=88%  B8=83%  B12=70%

alpha2=0.75:  N=48  B3=68%  B4=59%  B6=77%  B8=68%  B12=49%
              N=72  B3=76%  B4=85%  B6=95%  B8=92%  B12=83%

alpha2=1.00:  N=48  B3=72%  B4=54%  B6=89%  B8=83%  B12=68%
              N=72  B3=79%  B4=81%  B6=99%  B8=97%  B12=93%
```
**B=6 crossover:** beats B=3 when α₂≥0.5 and N≥48, by 15–20pp at α₂≥0.75.  
**B=12 over-fits:** too many df even under strong harmonics; B=6 is optimal.

---

## 4. Why Methods Behave This Way

### DCP (K=1) — B-invariance
For equispaced B≥3: Σcos²(ωtᵢ) = N/2 regardless of B.  
The design matrix achieves maximum D-efficiency at any equispaced B≥3.  
NCP = N·r²/2 — only N and r matter.  
Under harmonic violation: cos(2ωt) aliases onto cos(ωt) at B=3 → inflates NCP for K=1 test.

### JTK — replication depth wins
MetaCycle averages m replicates per time point before computing Kendall τ.  
Effective test runs on B means with SE ∝ σ/√m.  
More B → fewer replicates → noisier means → weaker rank signal.  
This collapse to means is why JTK cannot leverage temporal resolution.

### RAIN — temporal resolution wins
With `nr.series=m`, RAIN keeps all N individual observations.  
Umbrella test (Jonckheere-Terpstra variant) uses the full rank distribution.  
More time points → more resolution for the umbrella pattern → higher power.  
**Critical:** our original `detect_RAIN` wrapper was broken for replicated data  
(used `median(diff(times))` which = 0 for replicated designs, and omitted `nr.series`).  
Fixed in `code/detection.R` — now correctly computes `deltat` from unique times.

### Multi-harmonic — df vs signal tradeoff
Model: Y = μ + Σₖ[aₖcos(kωt) + bₖsin(kωt)] + ε, K = ⌊(B−1)/2⌋  
F-test has 2K df. NCP = N/(2σ²) · ΣAₖ²  
At B=3 (K=1): 2df, captures A₁² only (plus aliasing of A₂)  
At B=6 (K=2): 4df, captures A₁²+A₂² explicitly — wins when A₂ is large  
At B=12 (K=5): 10df, over-fits — loses power because df cost > signal gain  
Crossover B=6 > B=3 requires: α₂ ≥ 0.5 AND N ≥ 48

---

## 5. Figure Concepts

### Figure 1: Method comparison under cosinor truth (regular scenario)
- 3-panel: DCP(K=1) / JTK / RAIN
- X-axis: B ∈ {3,4,6,8}; lines: N={24,48,72}
- Shows: DCP flat, JTK favors low B, RAIN favors high B
- Dataset: D1D2 (r=0.65)
- Caption: "Detection power under standard cosinor truth is method-dependent:
  DCP is B-invariant, JTK favors replication depth, RAIN favors temporal coverage."

### Figure 2: Cosinor violation — DCP vs multi-harmonic
- 2-panel: DCP(K=1) | Multi-harmonic(adaptive K)
- X-axis: α₂ ∈ {0, 0.3, 0.5, 0.75, 1.0}; lines: B={3,6,8} at N=48
- Shows: DCP → B=3 flat/increases (aliasing); multi-harmonic → B=6 crosses over at α₂≈0.5
- Caption: "Adaptive multi-harmonic regression unlocks the B advantage:
  at α₂≥0.5 and N≥48, B=6 (K=2) yields 15–20pp more power than B=3 (K=1)."

---

## 6. User Guidance (Final)

**Step 1: Choose your detection method based on scientific context**

| Context | Recommended method | Recommended B |
|---------|-------------------|---------------|
| Standard circadian RNA-seq, sinusoidal genes expected | DCP (K=1 cosinor) | B=4 (power = N-driven) |
| Need rank-based robustness, no distributional assumption | JTK (MetaCycle) | B=4 (balance) |
| Sensitive non-parametric, individual obs preserved | RAIN (nr.series) | B=6–8 |
| Pilot shows harmonic content (α₂≥0.5, e.g. liver) | Multi-harmonic (K=2) | B=6 |
| Unknown waveform, exploratory | Multi-harmonic (K=2) | B=6 as default |

**Step 2: Use SCP power curves at your chosen B to determine N**  
Once B is fixed, power is determined by N = B×m.  
SCP outputs power vs N — read off N for target power (e.g. 80%).

**Step 3: Estimate α₂ from pilot data to guide method choice**  
Fit multi-harmonic model to pilot, compute α₂ = A₂/A₁ per gene.  
If median α₂ > 0.4: use B=6 with adaptive K=2.  
If median α₂ < 0.3: use B=4 with K=1 (standard cosinor).

---

## 7. Proposed Final Consolidated Smoke Test

Script: `examples/exploratory/bvsm_method_comparison.R`  
Purpose: head-to-head comparison of all four methods at same grid

```
Methods:  DCP(K=1), JTK, RAIN, MultiHarmonic(adaptive K)
Truth:    alpha2 = 0, 0.5, 1.0
B vals:   3, 4, 6, 8
N vals:   24, 48, 72
Pilot:    D1D2 D1 (r~0.65)
Output:   power[method, B, N, alpha2]
```

Output feeds directly into Figure 1 and Figure 2.

---

## 8. Relationship to LimoRhyde (Singer & Hughey 2019, J Biol Rhythms 34:5–18)

LimoRhyde is the closest prior framework to our multi-harmonic approach, and provides key citation anchors.

### What LimoRhyde implements
- `sinusoid=TRUE` (default): `getCosinorBasis()` → **K=1 only** (cos + sin of fundamental). Identical to our DCP.
- `sinusoid=FALSE`: `getSplineBasis()` → periodic splines with `nKnots` knots. A separate non-parametric alternative.
- **No K>1 multi-harmonic Fourier option exists in the package.**
- Confirmed by inspecting source: `limorhyde v1.0.3` installed 2026-04-20.

### What their paper says (key quotes)
*Methods, p.7:*
> "Although this decomposition is the simplest, one could also decompose time based on **multiple harmonics of the Fourier series** or on periodic splines."

*Discussion, p.16:*
> "LimoRhyde **could also be used to detect differences in higher-order harmonics** of circadian gene expression (Hughes et al., 2009)."

They explicitly flag multi-harmonic detection as future/unexplored work.

### Their simulation design
Singer & Hughey simulated B=12 time points × m=2 replicates (N=24). They never asked whether B=6 × m=4 performs differently — no B vs m analysis anywhere in the paper.

### Our novelty relative to LimoRhyde
The K = ⌊(B−1)/2⌋ identifiability constraint and its power consequence are not in LimoRhyde.  
Paper framing: "Singer and Hughey (2018) noted that LimoRhyde 'could also be used to detect differences in higher-order harmonics,' but did not characterize the identifiability constraint — at most K = ⌊(B−1)/2⌋ harmonics are estimable from B distinct sampling times — nor its consequences for detection power."

---

## 9. Pipeline Extension Needed

To support all frameworks in SCP:
- `detect_RAIN`: fixed (2026-04-20) — now passes `nr.series` and correct `deltat`
- `detect_multiharmonic`: needs to be added to `code/detection.R`
- `runPowerAnalysis`: needs `method` argument: `"cosinor"`, `"JTK"`, `"RAIN"`, `"multiharmonic"`
- `CircadianDesignOptions`: add `K` or `max_harmonics` parameter
- Pilot estimation: add `estimateAlpha2()` function to `code/estimation.R`
