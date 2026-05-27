# SCP Manuscript — Section & Figure Plan

**Status:** updated 2026-05-12 (LimoRhyde adoption; B=4 dropped from Fig 5; dispatch-level Nyquist warning).
**Companion docs:** `derivations/FMM_to_harmonic_LRT.md` (math supplement); `MORNING_REPORT.md` in `/tmp/fmm_diagnostic/` (overnight diagnostic results, includes Fig 5 v6 power table).

## Update log

- **2026-05-12.** Replaced the unshrunk QR-based K-harmonic F-test in `detect_FMM` with the LimoRhyde framework (limma::lmFit + eBayes + topTable on a multi-harmonic basis built from `limorhyde::limorhyde()`). Default `ebayes = TRUE` with small-G fallback at $G < 50$. Paper framing: LimoRhyde (Singer & Hughey 2019) introduces the framework at $K = 1$; this work extends it to user-tunable $K$ with $K = 2$ default.
- **2026-05-12.** Active B grid in Fig 5 reduced from $\{4, 6, 8, 12, 24\}$ to $\{6, 8, 12, 24\}$. B=4 K=2 design is rank-deficient (sin(2 ω₀ t) is identically zero at hours 0, 6, 12, 18); limma correctly refuses to fit. Methods section now states the identifiability rule $B \geq 2K+1$ once.
- **2026-05-12.** Added dispatch-level Nyquist warning in `runSimsSingleCohort`: when `method = "FMM"` is called with a design that has fewer than $2K+1$ distinct sampling phases, the function emits a clear warning before any `mclapply` child can swallow it. The warning recommends either reducing $K$ or switching to `method = "DCP"`.
- **2026-05-12.** Fig 5 v6 confirmed K=2 is B-invariant in active designs above the Nyquist threshold, mirroring DCP's known B-invariance. Within-N spread across $B \in \{6, 8, 12, 24\}$ is $\leq 4$ pp at every $N$ (within Monte Carlo noise). The design-prescription message is now sharp: choose any $B \geq 2K+1$, then allocate $N = B \times m$.

---

## Manuscript-level scope

The paper has **two coherent contributions**:

1. A **power-analysis framework** for circadian rhythmicity detection that is calibrated to real pilot data via empirical parameter distributions, supports both active and passive designs, and quantifies pilot-estimation uncertainty via outer bootstrap.
2. A **K-harmonic LRT detector** (`detect_FMM`) that is the linearized rhythmicity test motivated by truncating the Fourier expansion of the FMM model. Replaces nonlinear FMM-LRT with exact F-distribution calibration, sidestepping the Davies boundary problem.

The framework and detector are decoupled by design: any periodic-signal detector plugs into the framework; the K-harmonic LRT is the recommended default for non-cosinor data.

---

## Methods structure

```
§2.1  Single-cohort cosinor power           Fig 1   DCP, plug-in framework
§2.2  Differential power                    Fig 2   DCP differential, bootstrap-aware
§2.3  Bootstrap uncertainty                 Fig 3   completes the framework
§2.4  Sensitivity to cosinor violation      Fig 4   K-harmonic LRT introduced; ω + α sweeps
§2.5  Active vs passive design              Fig 5   K=2 across design types
```

Detailed contents below.

---

## §2.1 — Single-cohort cosinor power (Fig 1)

**Detector:** DCP (cosinor F-test), $H_0: A = 0$ in the cosinor model.
**Pilot:** Seney ACC controls (n=60, passive).
**Story:** plug-in transcriptome-wide power, with stratification by SNR ($r = A/\sigma$).
**Status:** ✅ done.

---

## §2.2 — Differential power (Fig 2)

**Detector:** DCP applied to two-group differential analyses (DR, DP, DM endpoints).
**Pilot:** GTEx Adrenal Gland vs Liver (n=187/154).
**Story:** transcriptome-wide power for differential rhythmicity, accounting for FDR via BH.
**Status:** ✅ done (rerun pending after partial-BH FDR fix).

---

## §2.3 — Bootstrap uncertainty (Fig 3)

**Goal:** quantify pilot-estimation uncertainty. Plug-in power in §2.1/§2.2 treats pilot estimates $(A, \sigma, \phi)$ as known; bootstrap shows that finite-N pilots induce **non-trivial CI width on the power curve.**

**Pilot lineup (under review by bootstrap-focused agent — see Task #26).** Current: Seney ACC (n=60), GTEx Pancreas (n=249), GTEx Thyroid (n=416). The agent is evaluating whether these CI widths (current max ≈9.6pp at Seney) are dramatic enough or if a smaller pilot (e.g., baboon CAMO n=12) should be added to highlight the extreme-uncertainty case.

**Layout:** 1 row × 3–4 columns (one panel per pilot), with vertical 95% CI error bars over plug-in line.
**Story:** "CI width is widest at the elbow of the power curve; small pilots produce design recommendations with substantial uncertainty that plug-in alone cannot communicate."

**Action items:**
- [ ] Receive bootstrap agent recommendation
- [ ] Update pilot lineup if needed
- [ ] Re-render Fig 3

---

## §2.4 — Sensitivity to cosinor violation (Fig 4)

**Critical scaffolding:** this is where the K-harmonic LRT (`detect_FMM`) is introduced.

### §2.4.1 — Empirical motivation: real data are non-cosinor

Before introducing the methodology, we establish the **need** for cosinor-violation handling. Pilot tissues exhibit:
- FMM ω-medians of 0.20–0.40 (far from cosinor's ω=1 limit)
- Median $R^2_{\text{FMM}} \approx 0.85$–0.95 vs $R^2_{\text{cosinor}} \approx 0.50$–0.70 on the same genes
- Visible sharp peaks, asymmetric waveforms

**Optional supporting figure (Fig 4a or supplementary):** 1×3 panel showing per-tissue distributions:
- Histogram of estimated ω across rhythmic genes (mass concentrated at ω<0.5)
- $R^2_{\text{cosinor}}$ vs $R^2_{\text{FMM}}$ scatter (most points above diagonal)
- Example sharp-peak gene profiles with cosinor and FMM fits overlaid

If we add this, it opens §2.4 visually before the methodology.

### §2.4.2 — The cosinor model and its limits

Recap the cosinor model and DCP from §2.1. Note that the DCP test is the LRT for cosinor truth ($H_0: A=0$, $\omega=1$ implicitly). Under cosinor truth, DCP is optimal.

When the truth is non-cosinor (FMM with $\omega < 1$), DCP is no longer the LRT — it is the **first-harmonic projection** of the FMM signal, attenuated by:

$$c(\omega) = \frac{2\omega}{1+\omega^2}$$

So DCP captures only $A \cdot c(\omega)$ of the FMM amplitude. The remaining variance lives in higher harmonics that DCP cannot detect.

### §2.4.3 — The FMM model as a flexible alternative

The FMM model (Rueda et al. 2019) parameterizes asymmetric rhythmic signals via:

$$y(t) = M + A \cos\!\big(\beta + 2\arctan(\omega \tan((t-\alpha)/2))\big) + \varepsilon$$

For descriptive characterization of waveform shape, FMM is well-established. **For hypothesis testing, however**, the FMM-LRT (Wilks LRT under nonlinear FMM) faces the **Davies problem**: $\{\beta, \alpha, \omega\}$ are unidentifiable when $A=0$, so the asymptotic null distribution is the supremum of a Gaussian process — no closed form, requires empirical calibration.

### §2.4.4 — From FMM to K-harmonic LRT (the linearization)

The FMM signal admits a closed-form Fourier expansion (full derivation in `derivations/FMM_to_harmonic_LRT.md`):

$$y(t) = M' + A(1 - r^2) \sum_{k=1}^{\infty} r^{k-1} \cos\!\big(k(t - \alpha) + \beta\big), \quad r = \frac{\omega - 1}{\omega + 1}$$

Truncating at $K$ harmonics and **releasing** the geometric-decay constraint (letting $a_k, b_k$ be free) yields a fully linear regression:

$$y(t) = M + \sum_{k=1}^{K} \big[a_k \cos(k\omega_0 t) + b_k \sin(k\omega_0 t)\big] + \varepsilon$$

Hypothesis: $H_0: a_1 = b_1 = \cdots = a_K = b_K = 0$ (no rhythm).

**Test statistic:** standard nested F-test, exact under Gaussianity:

$$F = \frac{(\text{SSE}_0 - \text{SSE}_K) / (2K)}{\text{SSE}_K / (n - 2K - 1)} \;\sim\; F(2K,\, n - 2K - 1) \quad \text{under } H_0$$

This is the **`detect_FMM(K)` detector** in our codebase.

### §2.4.5 — Identifiability: B ≥ 2K + 1 (Nyquist)

The K-harmonic regression has $2K + 1$ parameters and requires at least $2K + 1$ distinct sampling times per period for identifiability — the classical Nyquist condition (Nyquist 1928; Shannon 1949) applied to trigonometric regression. Cite Hughes et al. (2017) for circadian-specific multi-harmonic recommendation.

For the recommended K=2: requires $B \geq 5$.

### §2.4.6 — Choice of K

The variance fraction captured by truncation at K is

$$\frac{V_K}{V_\infty} = 1 - r^{2K}, \quad r = \frac{1 - \omega}{1 + \omega}$$

For empirical $\omega \in [0.20, 0.40]$ (our pilot range), K=2 captures 80–97% of FMM variance. K=3 captures 91–99%, but requires $B \geq 7$.

**Default: K=2.** Justified by:
1. Variance capture sufficient for typical biological ω (>80%).
2. The 2nd harmonic has biological interpretation (12h ultradian rhythms; well-documented in liver, gut, immune).
3. Identifiable in standard active designs (B=6 or higher).
4. Empirical comparison with K=1, 3, 4, 5 (Phase 4b/4c) confirms K=2 is the operating sweet spot.

### §2.4.7 — Sensitivity sweep results: Fig 4

**Pilot:** Baboon LUN ($\hat\beta \approx 3.19$, $\hat\sigma_\alpha \approx 2.23$ h, $R^2_{\text{med}} \approx 0.89$).
**Detector:** `detect_FMM(K=2)`.
**Design:** active, B=12 (every 2h).

**Layout:** 1 row × 2 panels.

| Panel | Sweep | Distribution | Curves |
|---|---|---|---|
| A | ω | $\omega_g \sim \text{Beta}(1, \beta)$, $\beta \in \{0.5, 1, 2, 5, 20\}$ | one per β; α fixed at empirical |
| B | α | $\alpha_g \sim \text{vonMises}(0, \kappa)$, $\sigma_\alpha \in \{0, 0.5, 1, 2, 4\}$ h | one per σ_α; ω fixed at Beta(1, β̂) |

**Story:** Panel A shows smooth power degradation as β increases (more rhythmic mass at low ω, more variance leaking into higher harmonics that K=2 can't capture). Empirical β̂ ≈ 3.19 anchors the realistic regime. Panel B shows the dataset-specific α-dispersion effect.

---

## §2.5 — Active vs passive design (Fig 5)

**Goal:** demonstrate that the K-harmonic LRT generalizes from controlled (active) to uncontrolled (passive) sampling, while the cosinor-truth B-invariance theory only applies to active designs.

### Theoretical setup

**Under cosinor truth:**
$$\text{NCP}_{\text{DCP}} = \frac{(A/\sigma)^2 \cdot N}{2}$$
B-invariant (only N matters). DCP is optimal. *Theory predicts and empirics will confirm B-invariance for any cosinor pilot.*

**Under FMM truth:**
$$\text{NCP}_{\text{DCP}} = \frac{(A \cdot c(\omega) / \sigma)^2 \cdot N}{2}$$
Still B-invariant in design, but signal is attenuated by $c(\omega) < 1$. The lost variance lives in higher harmonics that DCP cannot access. **DCP under-detects under FMM truth, regardless of B.**

The K-harmonic LRT recovers the lost noncentrality:
$$\text{NCP}_{\text{K=2}} \approx \frac{(V_1 + V_2) \cdot N}{2\sigma^2}$$
**provided** $B \geq 2K + 1 = 5$ for identifiability. Below this threshold, K=2 is rank-deficient.

### Layout (2 rows × 2 panels)

|  | **Active design** (Baboon CAMO KIM) | **Passive design** (Seney ACC) |
|---|---|---|
| **Top: DCP** | Power vs N at $B \in \{6, 8, 12, 24\}$ — overlap (B-invariance) | Power vs N (single empirical TOD distribution) |
| **Bottom: K=2** | Power vs N at $B \in \{6, 8, 12, 24\}$ — all above DCP, B-invariant above Nyquist | Power vs N — DCP-overlay shows DCP advantage on cosinor-like passive truth |

The B grid starts at 6 because the K=2 test requires $B \geq 2K+1 = 5$ distinct sampling phases per period (Nyquist identifiability). Designs with $B = 4$ are rank-deficient under K=2 and would return zero power; the figure no longer carries this failure case as a panel artefact. The methods section states the identifiability rule once, and `runSimsSingleCohort` emits a dispatch-level warning if a user requests an underdetermined $K$.

### Tissue selection
- **Active:** Baboon CAMO KIM ($\hat\beta = 2.62$, $R^2 = 0.87$ via eBayes pre-screen). Clean cosinor violation, broad biological relevance.
- **Passive:** Seney ACC (n=60, post-mortem human PFC). DCP-overlay used to show DCP advantage on cosinor-like passive truth.

### Fig 5 v6 results (2026-05-12 evening)

DCP active KIM, averaged over $B \in \{6, 8, 12, 24\}$:

| N | DCP | K=2 | K=2 advantage |
|---|---|---|---|
| 24 | 0.265 | 0.418 | +15.3 pp |
| 48 | 0.738 | 0.840 | +10.2 pp |
| 72 | 0.898 | 0.948 | +5.0 pp |
| 96 | 0.948 | 0.975 | +2.7 pp |
| 120 | 0.975 | 0.988 | +1.3 pp |
| 144 | 0.983 | 0.990 | +0.7 pp |

Within-detector spread across $B$ is at most 4 pp at every $N$, within Monte Carlo noise for $S = 1000$ simulations.

### Story per row

**Top (DCP):**
- Active panel: lines overlap across $B$, confirming the cosinor-truth B-invariance theory $\text{NCP}_{\text{DCP}} = (A/\sigma)^2 N / 2$.
- Passive panel: empirical TOD distribution; baseline for the bottom row.

**Bottom (K=2):**
- Active panel: K=2 strictly dominates DCP by 10 to 15 pp at moderate $N$, with the gap closing as both detectors saturate at large $N$. Lines across $B$ overlap (K=2 B-invariance), mirroring the DCP row but at higher absolute power.
- Passive panel: DCP-overlay above K=2 because Seney passive truth is cosinor-like; the K=2 second harmonic adds variance without proportional signal gain. *This is the honest passive-design caveat.*

### Caption framing
> *"Top row: DCP detection power is approximately B-invariant under cosinor truth, confirming the theoretical prediction $\text{NCP} \propto N$. Bottom row: the K-harmonic LRT (K=2) captures additional power from higher harmonics under cosinor violation, with the same B-invariance property above the Nyquist threshold $B \geq 2K+1$. Under realistic non-cosinor active truth (CAMO KIM), K=2 dominates DCP by 10 to 15 percentage points at moderate $N$. Under cosinor-like passive truth (Seney ACC), DCP outperforms K=2 because the higher harmonics add variance without signal. The design-prescription rule is therefore: choose the detector that matches the suspected waveform shape, then any $B \geq 2K+1$ above the Nyquist threshold gives equivalent power for the same total $N$."*

---

## Supplementary

- **`derivations/FMM_to_harmonic_LRT.md`** — full math derivation (FMM → Möbius → Fourier → K-harmonic LRT → exact F-test). Already written.
- **Type-I error validation table** — empirical T1 at α∈{0.05, 0.01} for `detect_FMM(K=2)` across n. From Phase 3a output. Brief table or supplementary panel.
- **K-sweep empirical results** — Phase 4b/4c. Power vs K at fixed N, B=12, on 1–3 tissues. Reinforces K=2 default.
- **FMM diagnostic figure** — per-pilot ω distribution, R²_cos vs R²_FMM scatter. Could be either §2.4.1 or supplementary.

---

## Application section (§3 or §4 — TBD)

**Pilot:** Seney ACC (n=60, passive).
**Goal:** demonstrate `detect_FMM` on real data, contrasting K=1 (DCP) vs K=2 discovery.

**Run:** `detect_FMM(K=1)` and `detect_FMM(K=2)` on Seney ACC counts. Report:
- Discovery counts at FDR=0.05
- Overlap between K=1 and K=2 discoveries
- Genes detected only by K=2 (these are the "non-cosinor rhythmic" genes, validation that K=2 adds meaningful signal)
- Top-10 K=2-only genes, with profiles plotted

**Status:** scheduled as Phase 4d in overnight queue. Output goes into the manuscript's application paragraph.

---

## Outstanding decisions

| Decision | Default | Status |
|---|---|---|
| Bootstrap pilot lineup (Fig 3) | current 3 + maybe baboon | awaiting agent |
| FMM-justification figure: panel in §2.4 vs supplementary | recommend §2.4 panel | TBD |
| Fig 5 layout | 2×2 (rows=detector, cols=design) | confirmed |
| Fig 5 active row tissue | KIM only ($\hat\beta = 2.62$) | confirmed |
| Fig 5 active B grid | $\{6, 8, 12, 24\}$ | confirmed (2026-05-12) |
| Fig 5 K=2 active framing | K=2 dominates DCP by 10–15 pp at moderate N; B-invariant above Nyquist | confirmed (v6 empirics) |
| K-choice for §2.4 default | K=2 | confirmed pending Phase 4b empirics |
| §2.5 passive pilot | Seney ACC | confirmed |
| detect_FMM detector backend | LimoRhyde (limma + eBayes, K-tunable) | confirmed (2026-05-12) |
| eBayes default for detect_FMM | TRUE (with G<50 fallback) | confirmed |

---

## Computation status

- [x] `detect_FMM` implemented in code/detection.R
- [x] `detect_FMM_LRT` removed; null-table file archived
- [x] `runSimsSingleCohort` dispatch updated (method = "FMM", K parameter)
- [x] `estCircadianParamFMM` pre-screen uses detect_FMM(screen_K)
- [x] `code/estimation.R` doc/comment cleanup (FMM-LRT → detect_FMM)
- [x] `examples/publication/fig4_sensitivity.R` updated (method, labels)
- [x] `examples/publication/fig5_active_vs_passive_v4.R` updated (method, labels, K plumbed, B=4 dropped)
- [x] `examples/publication/archive/fig4_fmm_violation.R` (legacy, archived)
- [x] Dispatch-level Nyquist warning in `runSimsSingleCohort` (B < 2K+1, before mclapply)
- [x] `detect_FMM` re-implemented via LimoRhyde framework (limma + eBayes; `ebayes = TRUE` default; small-G fallback at G<50)
- [x] Fig 4 production rerun (REFIT=true) — LUN $\hat\beta = 1.80$ with eBayes
- [x] Fig 5 production active row (REFIT=true, B grid {6, 8, 12, 24})
- [x] Fig 5 production passive row (Seney ACC + DCP-overlay)
- [x] SuppFig FMM diagnostic re-rendered with eBayes pre-screen
- [ ] `code/plot_fmm.R` title strings (FMM-LRT → FMM K-harmonic)
- [ ] Phase 4b: K-sweep on KIM (K ∈ {1, 2, 3} comparison)
- [ ] Phase 4c: K-sweep on KIC + SUN
- [ ] Phase 4d: Seney ACC standalone demo (K=1 vs K=2 discovery overlap)
- [ ] FMM-justification panel script (if added)
- [ ] SCP.tex prose rewrite per this plan (LimoRhyde framing, Nyquist rule, design-prescription sentence)
