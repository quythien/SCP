# Application figures

Run from a repository checkout:

```sh
Rscript paper_reproduction/updated_application/redraw_all.R
```

This draws Figures 1–5 and the supplementary harmonic-framework figure in
`generated/`. Base R is sufficient. The plotting scripts and compact data preserve
the manuscript's fonts, colors and panel geometry. Each `redraw.R` can also be run
from this directory. Figure 4 retains the observed-data vector panels and redraws
the overlap and enrichment panels from saved results.

| Figure | Example | Plotting inputs | Full analysis |
|---|---|---|---|
| 1 | NAc and Caudate, passive | `Figure1/data/` | `Figure1/run_analysis.R` |
| 2 | Caudate–Putamen DR, DP and DM | `Figure2/data/` | `Figure2/run_analysis.R` |
| 3 | Caudate and BA47 donor bootstrap | `Figure3/data/` | `Figure3/run_analysis.R` |
| 4 | Observed Liver K = 1 and K = 2 fits | `Figure4/data/` and `Figure4/final/` | `Figure4/code/observed_liver.R` |
| 5 | Balanced Liver and Caudate–Putamen designs | `Figure5/data/`, `Figure2/data/` | `Figure5/run_analysis.R`, active Figure 2 cases |
| Supplement | Harmonic framework and matched detector comparison | `Supplement/data/` | `Supplement/run_analysis.R` |

## Recompute the simulations

The fast redraw does not rerun Monte Carlo experiments. Full runs use the saved
calibrations, recorded grids and seeds. They can take hours. Use a new output
directory for each run. `digest` and `jsonlite` are required for input verification.

```sh
Rscript paper_reproduction/updated_application/Figure1/run_analysis.R /path/to/new_figure1_results
Rscript paper_reproduction/updated_application/Figure1/summarize.R /path/to/new_figure1_results /path/to/new_figure1_panels

mkdir -p /path/to/frozen_SCP_library
R CMD INSTALL --library=/path/to/frozen_SCP_library paper_reproduction/updated_application/Figure2/analysis/SCP_v2_source.tar.gz
Rscript paper_reproduction/updated_application/Figure2/run_analysis.R /path/to/frozen_SCP_library /path/to/new_differential_results
Rscript paper_reproduction/updated_application/Figure2/summarize.R /path/to/new_differential_results /path/to/new_differential_report

Rscript paper_reproduction/updated_application/Figure3/run_analysis.R Caudate /path/to/Caudate_Control.rds /path/to/new_caudate_bootstrap
Rscript paper_reproduction/updated_application/Figure3/run_analysis.R BA47 /path/to/BA47_full146.rds /path/to/new_ba47_bootstrap
Rscript paper_reproduction/updated_application/Figure3/summarize.R /path/to/new_caudate_bootstrap /path/to/new_ba47_bootstrap /path/to/new_bootstrap_report

Rscript paper_reproduction/updated_application/Figure5/run_analysis.R /path/to/new_liver_balanced_results
Rscript paper_reproduction/updated_application/Figure5/analysis/03_summarize.R /path/to/new_liver_balanced_results /path/to/new_liver_report

Rscript paper_reproduction/updated_application/Supplement/run_framework.R /path/to/frozen_SCP_library /path/to/new_framework_results
Rscript paper_reproduction/updated_application/Supplement/run_analysis.R /path/to/new_matched_results
Rscript paper_reproduction/updated_application/Supplement/matched_analysis/03_summarize.R /path/to/new_matched_results /path/to/new_matched_report
```

Figure 1 uses 600 experiments per sample size and 5,000 simulated genes. Set
`SCP_REPLICATES` and `SCP_N` only for a smaller runtime check. The full run uses the
recorded settings. Figure 2 includes 600 experiments in each of 54 cells, with N
measured per group. An optional list of case IDs limits the run. The frozen engine
keeps this numerical calculation independent of later package changes.
Figure 5 and Supplement accept optional job IDs after the output path.

Figure 3 needs the prepared full-cohort expression matrices. Its default runs the
original pilot and 100 donor bootstraps with 100 inner experiments per grid point.
Optional outer indices after the output path allow individual checkpoints to be
checked. The expression, time and gene hashes in `Figure3/data/configuration.rds`
identify the inputs without distributing donor records.

`Figure3/prepare_BA47.R` prepares GSE71620 from its phenotype table and all 420 CEL
files. It applies joint RMA, maps probes with `hugene11sttranscriptcluster.db`, keeps
the maximum-IQR probe per gene and selects the 146 BA47 donors with recorded TOD.
Adapted source preparation blocks are retained in `Figure3/code/`. Annotation
versions must match the recorded 19,998-gene input for exact numerical replay.
Caudate uses the full GSE160521 control matrix with 49,135 fitted genes and 59 donors.
The source preparation is `examples/pilots/native/02_human_local.R`: match control
donors to corrected TOD, remove zero-variance genes, calculate library-size CPM,
and use log2(CPM + 1). The smaller bundled expression fixture is not this full matrix.

`Figure4/run_analysis.R` takes an authorized GTEx v10 CPM RData file and a new
output directory. It uses the saved fit and enrichment results. No individual GTEx input
matrix or donor time file is distributed here. The current package's calibration
settings need not reproduce the fixed Liver tuple distribution used in these
figures. Those tuples are included in the simulation inputs.

The supplementary harmonic-framework figure uses its own calibration and N80
of 248. Both panels preserve paired harmonic parameters and noise. Input times
use up to 24 equally spaced points, with remaining samples allocated across
those points. Panel A retains an unavailable N = 300 cell. Panel B displays the
two occupied signal-to-noise ratio strata and reports the valid replicate counts
in its companion CSV. The separate matched-detector table uses the same two-harmonic truth for
both tests. Its rounded N80 values are 283/253 for passive sampling and 278/249
for balanced sampling. These are different analyses.
For a balanced design, round the planned sample size up to a multiple of B.
At B = 6, the matched comparison therefore gives feasible totals of 282 and 252.

Power and FDP summaries retain every completed cell. BH levels are nominal
analysis thresholds. In particular, the Caudate–Putamen DM result at N = 20 has
mean FDP about 0.393 at the 0.20 threshold. The figures do not establish empirical
FDR control at every grid point.

`source_manifest.json` records source identities. `SHA256SUMS` verifies the bundled
files. Earlier paper figures remain in the parent directory.

The framework replay retains the recorded seed formula. Its N = 300 seeds exceed
the integer range accepted by `set.seed`, so that cell is recorded as invalid.
The N80 crossing uses the valid N = 200 and N = 250 results.
