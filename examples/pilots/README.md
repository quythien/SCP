# GTEx pilot analyses

Run these analyses from the SCP source directory. Each script is divided into
numbered sections so the inputs, fitting, and checks can be followed separately.

## 1. Choose the local source files

These analyses use the R package `digest` to record source hashes.
Set `SCP_GTEX_CPM_NORM` to the v8 `CPM.all.norm.RData` file and
`SCP_GTEX_PHENO` to the matching `gtex_TOD_pheno.rds` file. Set
`SCP_GTEX_DONOR_ALIASES` to the local donor-alias CSV with `donor_original` and
`donor_resolved` columns. Its default path is
`data/local/gtex_v8_donor_aliases.csv`. This input stays local and is required so
the donor mapping cannot silently be omitted. The normalized
expression is used as supplied. Sample identifiers provide the donor key for
joining the phenotype table. The `CT` column supplies season-adjusted time in
hours. Tissue-site fields in sample identifiers are not time measurements.

The source adapter applies the documented donor aliases, keeps expression columns
in their original order, and excludes samples without verified CT from fitting.
It rejects duplicate phenotype donor keys and CT outside the 24-hour range.

## 2. Rebuild the database pilots

```r
source("examples/pilots/rebuild_gtex_database.R")
```

This stages K=1 and K=2 pilots for 31 v8 tissues. The Adrenal, Liver, and Muscle
files used by the paper remain unchanged. Existing v10 pilots retain their own
source release and time information. This analysis does not overwrite paper
figures, their caches, or their input files.

Outputs are staged under `output/pilot_rebuild/gtex_v8/` by default.
`SCP_GTEX_STAGE` can select a different staging directory. The staged manifest
records GTEx release and time source. The CSV update changes only rows for rebuilt
pilots. Other rows retain their existing values and each CSV keeps its own schema.
The public package summary is required. The two local paper CSVs are also updated
when they are present.

The fit uses seed 2025, period 24 hours, raw rhythm-screen p < 0.01, paired amplitude
and noise, and the existing harmonic-specific top-K setting. It uses the full
eligible tissue matrix. The build stages outputs for review before replacement. The older GTEx builder
entrypoints call this same staged analysis, including symbol annotations,
tested-gene universes, and source metadata. A GTEx-only legacy call does not save
its old manifest over the staged result. Mixed-species calls retain their other
species workflows. Existing v10 outputs and the three paper-linked tissues are
preserved by the legacy GTEx routes too.

The single-harmonic signal summaries use up to 300 genes ranked by raw cosinor
p-value among those passing p < 0.01. Signal/noise is amplitude divided by residual
standard deviation. Discovery counts use the full tested gene universe for BH
adjustment. The supplementary CSV phase columns contain arithmetic quantiles in
hours for the same selected genes.

## 3. Create local CPM alternatives

Set `SCP_GTEX_CHIRAL_PHASES` to the public CHIRAL `DIPs.RData` file, then run:

```r
source("examples/pilots/create_local_gtex_cpm.R")
```

This creates separate Adrenal, Liver, and Muscle files under
`data/local/gtex_v8_seasonal_CT/`. The original expression objects and identifiers
are preserved. `CT.by.tissue` contains named seasonal CT vectors aligned to the
expression columns. Sample metadata also contains raw and aligned CHIRAL phases.
Missing times remain missing. Consumers must read the explicit CT vector rather
than extracting time from sample identifiers.

CHIRAL phase origin and direction are aligned using all matched phenotype donors.
These inferred phases are a separate sensitivity analysis. They do not replace the
recorded seasonal CT used by the rebuilt database pilots.

Raw expression, donor-linked time metadata, and these local alternatives are not
included in the public repository. Only the intended derived pilot outputs and
reproduction code are prepared for the package update.
