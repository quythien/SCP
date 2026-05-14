# SCP: Simulation-Based Circadian Power Analysis

SCP is an R framework for sample-size planning and study design in circadian
transcriptomic experiments. Given a small pilot RNA-seq dataset, it answers:
*how many samples are needed, and how should the time points be allocated, to
detect rhythmic genes or differentially rhythmic genes at a target power and
controlled false discovery rate?*

The framework supports both active (controlled animal) and passive
(post-mortem human, clinical observational) designs, single-cohort
rhythmicity, two-group differential analysis (DR, DP, DM endpoints), and
non-parametric bootstrap for pilot-driven uncertainty.

For non-cosinor waveforms (sharp peaks, asymmetric rhythms common in real
data), SCP provides the K-harmonic likelihood-ratio test `detect_FMM`, a
linearised rhythmicity test motivated by truncating the Fourier expansion of
the Frequency Modulated Mobius model. The test is calibrated by an exact
F-distribution (no Davies-supremum approximation), so it sidesteps the
boundary-identifiability problem of the nonlinear FMM-LRT.

## Status

Source-and-run scripts plus a v0.4.0 R-package skeleton (DESCRIPTION,
NAMESPACE, R/). Companion manuscript in preparation (Pham et al., 2026).

## Install

Two install paths are supported.

**(a) As an R package (recommended for users):**
```r
remotes::install_local("path/to/PowerSim")
library(SCP)
```

**(b) Source-style for development (no install):**
```r
old_wd <- setwd("code"); source("setup.R"); setwd(old_wd)
```

Required: R >= 4.4. Imports: ggplot2, reshape2, limma, minpack.lm, nloptr,
parallel. Optional: FMM, MetaCycle, rain, circacompare, limorhyde,
patchwork, cowplot, gtable, Rcpp, RcppArmadillo, readxl.

## Quick start

```r
library(SCP)   # or: old_wd <- setwd("code"); source("setup.R"); setwd(old_wd)

# Load a bundled pilot
bio <- readRDS("data/gtex_Liver_single_pilot.rds")

# Specify design and analysis
des  <- CircadianDesignOptions(sample_sizes = c(40, 80, 120, 160),
                                nsims = 100L,
                                design = "passive", cts = bio$cts)
aopt <- CircadianAnalysisOptions(alpha = 0.05)

# Run single-cohort power (DCP for cosinor truth, FMM for cosinor violation)
res  <- runSingleCohortPower(bio, des, aopt, methods = "DCP", plot = FALSE)

# Sample size for 80 percent genome-wide power
npower(res, target_power = 0.80, fdr = 0.05)
```

## When to use DCP vs detect_FMM

| Pilot characteristic | Recommended detector |
|---|---|
| Cosinor-like waveforms (estimated omega-median > 0.6) | `detect_DCP` |
| Sharp or asymmetric peaks (omega-median < 0.5) AND design has B >= 5 distinct time points | `detect_FMM(K = 2)` |
| Sharp peaks but design has B < 5 | `detect_DCP` (K=2 is not identifiable; Nyquist condition violated) |
| Passive design with typical human cohorts (Seney-ACC-like) | `detect_DCP` (cosinor-truth is usually a good approximation; K=2 pays a df cost for no harmonic gain) |
| Active animal cohorts with strong non-cosinor signal | `detect_FMM(K = 2)` |

To diagnose your pilot's omega distribution, run
`estCircadianParamFMM(expr, times)` and inspect
`bio$diagnostics$omega_median`.

## Reproducing the paper figures

| Figure | Section | Script | Output PDF |
|---|---|---|---|
| Fig 1 | 2.1 single cohort | `examples/publication/fig1_single_cohort_power.R` | `output/main_figures/Fig1*.pdf` |
| Fig 2 | 2.2 differential | `examples/publication/fig2_differential_power_adr_liv.R` | `output/main_figures/Fig2_*.pdf` |
| Fig 3 | 2.3 bootstrap | `examples/publication/fig5_bootstrap_sc.R` | `output/main_figures/Fig3_bootstrap_singlecohort.pdf` |
| Fig 4 | 2.4 sensitivity | `examples/publication/fig4_sensitivity.R` | `output/main_figures/Fig4_sensitivity.pdf` |
| Fig 5 | 2.5 active vs passive | `examples/publication/fig5_active_vs_passive_v4.R` | `output/main_figures/Fig5_active_vs_passive.pdf` |
| Supp | FMM diagnostic | `examples/publication/supp_fmm_diagnostic.R` | `output/main_figures/SuppFig_FMM_diagnostic.pdf` |

For details on the methodology and section structure, see
[`paper/PowerSim/SECTION_PLAN.md`](paper/PowerSim/SECTION_PLAN.md).
The Fourier-expansion derivation behind the K-harmonic LRT lives in
[`paper/PowerSim/derivations/FMM_to_harmonic_LRT.md`](paper/PowerSim/derivations/FMM_to_harmonic_LRT.md).

## Key functions

| Function | Purpose |
|---|---|
| `estCircadianParam()` | Cosinor-based pilot fit; returns `CircadianBioOptions` |
| `estCircadianParamFMM()` | FMM-aware pilot fit (recommended for non-cosinor data) |
| `estCircadianParamTwoGroup()` | Two-group pilot for differential power |
| `prepCircadianData()` | Pre-normalisation helper (counts / CPM / log2) |
| `runSingleCohortPower()` | Single-cohort rhythmicity power; methods: DCP, FMM, MH, JTK, RAIN |
| `runDifferentialPower()` | Differential (DR / DP / DM) power |
| `runSingleCohortGrid()` | B vs m grid sweep (active designs only) |
| `runBootstrapDesignGrid()` | Outer bootstrap for pilot uncertainty |
| `detect_FMM(K = 2)` | K-harmonic LRT (the recommended non-cosinor detector) |
| `detect_DCP()` | Cosinor F-test (1-harmonic, the recommended cosinor detector) |
| `npower()` | Linear-interpolated sample size for a target power |
| `simCircadianFMM()` | Simulate single-cohort data under FMM truth |

JTK_CYCLE, RAIN, and the multi-harmonic detector remain in the package
via `detect_JTK`, `detect_RAIN`, and `detect_MH` for benchmarking, but are
not the primary recommended methods.

## Tests

```bash
Rscript tests/test_detect_FMM_typeI.R
```

verifies type-I error calibration for `detect_DCP` and `detect_FMM` at
K in {1, 2, 3} across n in {12, 24, 48, 72, 144}, plus the rank-deficient
edge case (B < 2K+1 must be conservative).

## Documentation

- [`doc/TUTORIAL.md`](doc/TUTORIAL.md): scenario-based usage guide
- [`vignettes/SCP_vignette.md`](vignettes/SCP_vignette.md): full user guide
- [`paper/PowerSim/SECTION_PLAN.md`](paper/PowerSim/SECTION_PLAN.md): paper structure and figure plan
- [`paper/PowerSim/derivations/FMM_to_harmonic_LRT.md`](paper/PowerSim/derivations/FMM_to_harmonic_LRT.md): derivation of K-harmonic LRT from FMM Fourier expansion
- [`NEWS.md`](NEWS.md): release notes

## Citation

A `CITATION.cff` is provided. Until the manuscript has a DOI, please cite
the software as:

> Pham, T. (2026). SCP: Simulation-Based Circadian Power Analysis with
> K-harmonic LRT (v0.4.0).
> https://github.com/qtp/SCP

BibTeX:

```bibtex
@software{Pham_SCP_2026,
  author  = {Pham, Thien},
  title   = {{SCP: Simulation-Based Circadian Power Analysis with
              K-harmonic LRT}},
  year    = {2026},
  version = {0.4.0},
  url     = {https://github.com/qtp/SCP}
}
```

## License

MIT. See [`LICENSE.md`](LICENSE.md).

## Contact

Thien Pham, University of Pittsburgh, <quythien14@gmail.com>.
