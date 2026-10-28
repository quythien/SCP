# GTEx v10 builders

The current GTEx pilots use normalized v10 expression and season-adjusted time
of death (TOD). The same pilots supply the package and app. SCP fits gene
acrophase relative to TOD and uses the fitted pilot to simulate the planned
sampling design.

## 1. Set the source files

Set `SCP_GTEX_V10_FULL` to the within-tissue normalized log-CPM cache and
`SCP_GTEX_PHENO` to the local donor phenotype file. The phenotype table must
contain unique `SUBJID` values and seasonal `CT` in hours. Sample donor keys
are matched to that table. The GTEx study accession is phs000424.

```r
Sys.setenv(
  SCP_GTEX_V10_FULL = "path/to/v10_CPM_full.RData",
  SCP_GTEX_PHENO = "path/to/gtex_TOD_pheno.rds",
  SCP_GTEX_TOD_STAGE = "output/gtex_v10_normalized_TOD"
)
```

When the cache needs to be created, also set `SCP_GTEX_V10_READS` and
`SCP_GTEX_V10_SAMPLE_QC` to the official v10 gene-read GCT and sample attributes.
Their public URLs and hashes are recorded in `public_inputs.json`.

```r
Sys.setenv(
  SCP_GTEX_V10_READS = "path/to/v10_gene_reads.gct.gz",
  SCP_GTEX_V10_SAMPLE_QC = "path/to/v10_sample_attributes.txt"
)
```

Within each tissue, genes with mean raw count above 10 undergo TMM
normalization and `edgeR::cpm(log = TRUE, prior.count = 0.25)`. This is the
expression normalization step of Talamanca et al. (2023). An existing verified
cache is used directly.

## 2. Build the normalized TOD pilots

Work from the SCP source root and run:

```r
source("examples/pilots/gtex_v10/10_build_normalized_TOD.R")
```

The builder retains all samples with finite seasonal TOD. It keeps the first
gene for each symbol and applies mean log-CPM greater than 3 within that
sample cohort. Cosinor models are fitted to normalized expression.

Both harmonic variants are built for 46 tissues. All retained genes enter
fitting and multiple-testing calculations. The pilot effect-size distribution
uses up to 300 genes with the smallest p-values among genes passing nominal
p < 0.01. The summary records the number selected separately from the tested
gene count.

`SCP_GTEX_TOD_STAGE` defaults to `output/gtex_v10_normalized_TOD`. The output
contains `inst/extdata/pilots/`, the signal summary in `summary/`, and private
matching records. Pilot objects contain fitted gene parameters and aggregate
sampling distributions. The hourly TOD histogram records counts in 24 bins.
Individual recorded TOD values and source sample identifiers stay in the
local private inputs.

## 3. Compare sample sizes

Run the saved pilots through the sample-size comparison:

```r
source("examples/pilots/gtex_v10/11_compare_normalized_TOD_sample_sizes.R")
```

The simulation records its pilot identities, selection settings, random seeds,
sample-size grid, and Monte Carlo counts. Results describe the fitted pilot,
input time, model, and specified sampling design.

## 4. Render the report

Build the report from the saved summaries and simulation results:

```r
source("examples/pilots/gtex_v10/12_render_normalized_TOD_report.R")
```

The report presents the tissue signal summaries and sample-size comparisons.
The builder records source hashes and output checksums for future reproduction.
Existing paper inputs and manuscript figures remain unchanged.

After integrating these pilots and installing SCP, rebuild the app demo:

```r
source("examples/pilots/gtex_v10/07_build_demo.R")
```

## Earlier paired research analysis

Scripts 00 through 06 retain the earlier paired DIP and TOD workflow. Scripts
08 and 09 retain its TOD publication and report. The research workflow uses
CHIRAL covariate adjustment and v10 DIP inference, with the strict expression
matrices supplying inference and the broader matrices supplying paired pilots.
Its settings use `SCP_GTEX_V10_STAGE`, `SCP_GTEX_V10_STRICT`,
`SCP_GTEX_V10_BROAD`, `SCP_CHIRAL_SOURCE`, and `SCP_CHIRAL_LIBRARY`.
The default research output directory is `output/gtex_v10_chiral_protocol`.

The original CHIRAL study used GTEx v8. The retained v10 research workflow uses
CHIRAL 0.1.0 at source revision
`1f78a9a60b6d18facfedf5923d43bae9fa79f0b0`. Recorded TOD does not enter its
CHIRAL inference or orientation. The source manifest maps recipe hashes to
files in `source_snapshots/`. CHIRAL has its own license.

The earlier full-expression and DIP archive remains available with SCP 0.5.0.
It contains adjusted expression, gene annotation, tissue, project keys, DIP,
status, provenance, and checksums. Recorded TOD, clinical fields, original
identifiers, and the key bridge are omitted. Earlier package versions are
available through GitHub releases.

## References

- Talamanca L, Gobet C, Naef F. 2023. Sex-dimorphic and age-dependent organization of 24-hour gene expression rhythms in humans. Science 379:478-483. https://doi.org/10.1126/science.add0846
- GTEx Consortium. 2020. The GTEx Consortium atlas of genetic regulatory effects across human tissues. Science 369:1318-1330. https://doi.org/10.1126/science.aaz1776
- CHIRAL source: https://github.com/naef-lab/CHIRAL/tree/1f78a9a60b6d18facfedf5923d43bae9fa79f0b0
- GTEx data portal: https://www.gtexportal.org/home/downloads/adult-gtex
