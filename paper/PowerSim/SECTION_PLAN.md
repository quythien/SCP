# SCP Manuscript — Section & Figure Plan

**Status:** updated 2026-05 to match the current manuscript (`SCP.tex`). The
framework is cosinor-only: a single-harmonic cosinor F-test and a K-harmonic
generalization (default K = 2, the 12-hour second harmonic), unified through
`detect_cosinor(K)`. FMM is not used; JTK_CYCLE, RAIN, and BIO_CYCLE are cited
only as prior single-condition detectors. All figures use human GTEx and the
Ketchesin Putamen cohort (GSE160521); no baboon or mouse pilots appear in the
main paper. `SCP.tex` is the source of truth.

## Methods (Section 2)

- **2.1 Single-cohort power** (`sec:single_cohort`). Cosinor model, rhythmicity
  F-test, noncentral-F power, the effect size r = A/sigma, the design factor
  d(phi), and the four-step pilot-calibrated simulation pipeline with BH-FDR.
- **2.2 Differential power** (`sec:differential`). DR, DP, DM endpoints via a
  data-derived joint rhythmicity classification; targeted power.
- **2.3 Bootstrap uncertainty** (`sec:bootstrap`). Outer subject bootstrap that
  propagates pilot-estimation uncertainty into the projected power curves.
- **2.4 K-harmonic cosinor extension** (`sec:cosinor_violation`). Extends 2.1 by
  adding harmonics: model, nested F-test, noncentrality
  lambda = (N/2) r_eff^2 with r_eff^2 = sum_k r_k^2, identifiability
  B >= 2K + 1 (B >= 5 for K = 2), and K = 2 as the default. No FMM / Mobius
  content.

## Results (Section 3)

| Fig | Subsection | Pilot(s) | Content |
|---|---|---|---|
| 1 | 3.2 `sec:sc_results` (`fig:sc_adrliv`) | GTEx Adrenal, Liver | Single-cohort power driven by the pilot effect-size distribution |
| 2 | 3.3 `sec:diff_results` (`fig:diff_adrliv`) | GTEx Adrenal vs Liver | Differential DR / DP / DM power |
| 3 | 3.4 `sec:bootstrap_results` (`fig:bootstrap_sc`) | Putamen SCZ controls (n=28), GTEx Thyroid (n=416) | Bootstrap uncertainty vs pilot size |
| 4 | 3.5 `sec:twoharm_results` (`fig:twoharm_demo`) | GTEx Liver | K-harmonic recovers non-sinusoidal rhythms (exemplars, Venn, KEGG) |
| 5 | 3.5 `sec:twoharm_power` (`fig:twoharm_power`) | GTEx Liver | K-harmonic operating characteristics (FDR sweep, r_eff strata) |
| 6 | 3.7 `sec:active_passive` (`fig:active_vs_passive`) | Putamen control (n=59) | Active-design power governed by total N above identifiability |

Section 3.1 (`sec:pilot_datasets`) introduces the pilots and the simulation
setup.

## Discussion (Section 4)

Synthesis of Figs 1-6 against the Section 2 derivations, then a four-step
decision tree for prospective design (`sec:decision_tree`): detector (K),
sampling grid (B >= 2K + 1), sample size (N_80 at the pilot r-tilde), and
bootstrap when n_pilot is small.

## Reproduction

Figure scripts live in `examples/publication/` and
`examples/publication/two_harmonic/`; see the figure table in the top-level
`README.md` for the script and output PDF of each figure.
