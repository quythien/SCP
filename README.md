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

For non-sinusoidal waveforms (sharp peaks, asymmetric or 12-hour ultradian
rhythms common in real data), SCP provides a two-harmonic cosinor test through
the unified detector `detect_cosinor(K)`. With `K = 2` it adds the 12-hour
second harmonic and tests the joint harmonic block by an exact nested
F-test; with `K = 1` it reduces to the standard single-harmonic cosinor
F-test. Identifiability requires `B >= 2K + 1` distinct sampling times per
period.

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

Required: R >= 4.4. Imports: ggplot2, reshape2, limma, limorhyde,
minpack.lm, nloptr, parallel, grid, grDevices, Rcpp (with LinkingTo
Rcpp, RcppArmadillo). Optional: MetaCycle, rain, circacompare, patchwork,
cowplot, gtable, readxl, tidyr.

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

# Run single-cohort power (K=1 single-harmonic; K=2 adds the 12-hour harmonic)
res  <- runSingleCohortPower(bio, des, aopt, methods = "DCP", plot = FALSE)

# Sample size for 80 percent genome-wide power
npower(res, target_power = 0.80, fdr = 0.05)
```

## Choosing the harmonic order K

Both detectors are exposed through `detect_cosinor(K)`: `K = 1` is the
single-harmonic cosinor F-test, `K = 2` adds the 12-hour second harmonic.

| Pilot characteristic | Recommended order |
|---|---|
| Approximately sinusoidal waveforms | `detect_cosinor(K = 1)` |
| Non-sinusoidal or 12-hour ultradian structure AND design has B >= 5 distinct time points | `detect_cosinor(K = 2)` |
| Non-sinusoidal but design has B < 5 | `detect_cosinor(K = 1)` (K = 2 is not identifiable; B >= 2K + 1 violated) |
| Typical passive human cohorts (near-cosinor) | `detect_cosinor(K = 1)` (K = 2 pays a degrees-of-freedom cost for little gain) |

## Reproducing the paper figures

| Figure | Content | Script | Output PDF |
|---|---|---|---|
| Fig 1 | Single-cohort power (GTEx Adrenal, Liver) | `examples/publication/fig1_single_cohort_power.R` | `submission/figures/Fig1*.pdf` |
| Fig 2 | Differential power (GTEx Adrenal vs Liver) | `examples/publication/fig2_paired_sims.R` | `submission/figures/Fig2_*.pdf` |
| Fig 3 | Bootstrap uncertainty (Putamen, GTEx Thyroid) | `examples/publication/fig3_bootstrap_2panel.R` | `submission/figures/Fig3_bootstrap_singlecohort.pdf` |
| Fig 4 | Two-harmonic discoveries (GTEx Liver) | `examples/publication/two_harmonic/fig4_twoharm_demo.R` | `submission/figures/Fig4_twoharm_demo.pdf` |
| Fig 5 | Two-harmonic operating characteristics (GTEx Liver) | `examples/publication/two_harmonic/fig5_twoharm_framework.R` | `submission/figures/Fig5_twoharm_framework.pdf` |
| Fig 6 | Active vs passive design (Putamen control) | `examples/publication/two_harmonic/fig6_v9_putamen_paired.R` | `submission/figures/Fig6_active_BvsM.pdf` |

For details on the methodology and section structure, see
[`paper/PowerSim/SECTION_PLAN.md`](paper/PowerSim/SECTION_PLAN.md).

## Data and reproducibility

The raw expression matrices (GTEx; Ketchesin striatum, GEO accession
GSE160521) are access-controlled or large and are not redistributed here.
Instead, the small pilot summaries derived from them are bundled so the
figures regenerate without the raw data:

- **Figures 1, 2, 6** run directly from bundled pilots in `data/`
  (`gse160521_*_ctrl_pilot.rds`, `gtex_adr_vs_liv_pilot_paired.rds`,
  `gtex_*_single_pilot.rds`). Figure 2 rebuilds its pilot from the raw GTEx
  matrix only if the bundled pilot is absent.
- **Figure 4** regenerates from cached per-gene fits, the discovery overlap,
  and the KEGG enrichment in `output/two_harmonic/results/`
  (`tissue_Liver_fits.rds`, `tissue_Liver_venn.rds`, `tissue_Liver_kegg.rds`).
- **Figure 5** uses the bundled two-harmonic pilot
  (`output/two_harmonic/results/pilot_2h_GTExLiver_topK300.rds`); the full
  power-sweep cache is large and not committed, but is regenerated by setting
  `RUN_SIMS <- TRUE` in the script.
- **Figure 3** uses a large bootstrap cache that is not committed; it is
  regenerated by rerunning the bootstrap script.

Regenerating a pilot from scratch (rather than from the bundled summaries)
requires the corresponding raw dataset under its original data-use terms.

## Key functions

| Function | Purpose |
|---|---|
| `estCircadianParam()` | Single-harmonic cosinor pilot fit; returns `CircadianBioOptions` |
| `estCircadianParam2H()` | Two-harmonic pilot fit (amplitudes and phases of both harmonics) |
| `estCircadianParamTwoGroup()` | Two-group pilot for differential power |
| `prepCircadianData()` | Pre-normalisation helper (counts / CPM / log2) |
| `runSingleCohortPower()` | Single-cohort rhythmicity power |
| `runDifferentialPower()` | Differential (DR / DP / DM) power |
| `runSingleCohortGrid()` | B vs m grid sweep (active designs only) |
| `runBootstrapDesignGrid()` | Outer bootstrap for pilot uncertainty |
| `detect_cosinor(K)` | Unified cosinor detector: K = 1 single-harmonic, K = 2 two-harmonic |
| `simCircadianSingleCohort()` | Simulate single-cohort data (single harmonic) |
| `simCircadianSingleCohort2H()` | Simulate single-cohort data with a 12-hour second harmonic |
| `npower()` | Linear-interpolated sample size for a target power |

`detect_DCP` and `detect_FMM` remain as the single- and multi-harmonic
backends called by `detect_cosinor`. JTK_CYCLE, RAIN, and the
multi-harmonic detector (`detect_JTK`, `detect_RAIN`, `detect_MH`) remain
in the package for benchmarking but are not the primary recommended methods.

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
- [`NEWS.md`](NEWS.md): release notes

## Citation

A `CITATION.cff` is provided. Until the manuscript has a DOI, please cite
the software as:

> Pham, T. Q. (2026). SCP: Simulation-Based Circadian Power Analysis
> (v0.4.0). https://github.com/quythien/SCP

BibTeX:

```bibtex
@software{Pham_SCP_2026,
  author  = {Pham, Thien Quy},
  title   = {{SCP: Simulation-Based Circadian Power Analysis}},
  year    = {2026},
  version = {0.4.0},
  url     = {https://github.com/quythien/SCP}
}
```

## License

MIT. See [`LICENSE.md`](LICENSE.md).

## Contact

Thien Quy Pham, University of Pittsburgh, <quythien14@gmail.com>.
