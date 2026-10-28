# GTEx v10 builders

These scripts prepare the active GTEx pilots for SCP. CHIRAL estimates donor
input time during data preparation. SCP estimates gene acrophase relative to
that time.

## Inputs

Use the official v10 gene-read GCT, public sample attributes, and public subject
phenotypes listed in `public_inputs.json`. The manifest records their URLs and
hashes. The GTEx study accession is phs000424.

Set `SCP_GTEX_V10_READS`, `SCP_GTEX_V10_SAMPLE_QC`, and
`SCP_GTEX_V10_SUBJECTS` to the source files. Set `SCP_GTEX_V10_FULL` to the
within-tissue log-CPM cache. The first normalization run creates this cache from
the gene reads. Set `SCP_CHIRAL_SOURCE` to the pinned upstream source directory.
Install the package in its `Pkg/CHIRAL` directory and set
`SCP_CHIRAL_LIBRARY` to the installation library. Set `SCP_GTEX_V10_STAGE` to
choose the output directory. Its default is `output/gtex_v10_chiral_protocol`.

Preprocessing follows the published CHIRAL workflow. Within each tissue, genes
with mean raw count above 10 undergo TMM normalization and
`edgeR::cpm(log = TRUE, prior.count = 0.25)`. Mean log-CPM thresholds are applied
before sample QC. The strict threshold is 3 and the broader threshold is 0.
The author functions adjust expression for ischemic-time group, public age
group, sex, and death type. The strict matrices supply phase inference. The
broader matrices supply the pilot fits.

Set `SCP_GTEX_V10_STRICT` and `SCP_GTEX_V10_BROAD` to the adjusted outputs.
The local season-adjusted TOD table is supplied through `SCP_GTEX_PHENO` for
paired TOD fitting. Recorded TOD does not enter CHIRAL inference or orientation.

The original CHIRAL study used GTEx v8. These scripts estimate DIP from v10
expression using CHIRAL 0.1.0 and the published preprocessing and consensus
functions. Source revision: `1f78a9a60b6d18facfedf5923d43bae9fa79f0b0`.

The source manifest maps recorded recipe hashes to the matching files in
`source_snapshots/`. Use the numbered scripts for new builds.

## Run the analysis

Work from the SCP source root and run the scripts in order:

```r
source("examples/pilots/gtex_v10/00_normalize_inputs.R")
source("examples/pilots/gtex_v10/01_infer_DIP.R")
source("examples/pilots/gtex_v10/02_build_pilots.R")
source("examples/pilots/gtex_v10/03_stage_full_data.R")
source("examples/pilots/gtex_v10/04_diagnostics.R")
source("examples/pilots/gtex_v10/05_compare_sample_sizes.R")
source("examples/pilots/gtex_v10/06_render_report.R")
source("examples/pilots/gtex_v10/07_build_demo.R")
```

The scripts write to a staging directory and record source hashes and analysis
settings. Existing paper inputs and manuscript figures remain unchanged.

## Input time

**CHIRAL-estimated DIP** is the default input time. CHIRAL uses clock-gene
expression across tissues to estimate donor-internal phase. The published
human muscle reference sets direction, with the PER3 peak anchored at nine
hours. Each tissue uses one seeded fit with up to 500 iterations. Native
convergence is recorded. Donor consensus uses the author's
`Phi.from.phi_mat(ct = 1.95)` defaults.

Cultured and transformed cells are excluded. Whole Blood receives donor DIP
from other tissues and does not contribute to inference. Samples without an
inferred donor DIP retain a missing value and a reason.

CHIRAL is an offline dependency with its own license. SCP does not redistribute
CHIRAL code under its MIT license.

**Season-adjusted TOD** is the second input-time example. Matching uses exact
v10 donor keys. Sample barcode fields are not interpreted as collection times.

## Pilots and summaries

Each paired comparison uses the same eligible samples, expression values, and
genes for DIP and TOD. K = 1 and K = 2 are fit separately. All eligible genes
enter fitting and multiple-testing calculations. Pilot selection is recorded
separately from that tested gene set.

The active manifest has one row per tissue, condition, and input-time basis.
The signal summary is `summary/gtex_v10_signal_summary.csv`. Pilot objects
contain fitted gene parameters and aggregate time-sampling distributions.
Individual recorded TOD values and source identifiers are omitted. Matching
tables remain local.

Diagnostics exclude inference clock genes and repeat donor consensus without
each target tissue. Comparisons report matched donor coverage. Tissues outside
the inference set are identified separately.

## Full expression archive

The archive retains all genes and samples in the adjusted analysis matrices
before pilot selection. It contains expression, gene annotation, tissue,
project keys, DIP in hours, status, provenance, and checksums. Samples without
DIP remain in the expression matrices. Recorded TOD, clinical fields, original
identifiers, and the key bridge are omitted.

## References

- Talamanca L, Gobet C, Naef F. 2023. Sex-dimorphic and age-dependent organization of 24-hour gene expression rhythms in humans. Science 379:478-483. https://doi.org/10.1126/science.add0846
- GTEx Consortium. 2020. The GTEx Consortium atlas of genetic regulatory effects across human tissues. Science 369:1318-1330. https://doi.org/10.1126/science.aaz1776
- CHIRAL source: https://github.com/naef-lab/CHIRAL/tree/1f78a9a60b6d18facfedf5923d43bae9fa79f0b0
- GTEx data portal: https://www.gtexportal.org/home/downloads/adult-gtex
