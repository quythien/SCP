# Native-time pilot recipes

These scripts prepare candidates for the unchanged non-GTEx database. The source index covers 114 base pilots and their 212 existing K = 1 and K = 2 files. It records source recipes, original script hashes, input hashes when known, the current pilot hashes, and the stored fitting settings. The new recipes have a synthetic smoke check. They have not been used to rebuild the full native database in this update.

Run from the repository root with SCP installed. Every output goes to `SCP_NATIVE_STAGE`, which defaults to `output/pilot_rebuild/native`. The scripts refuse to use the installed or source pilot directory as their staging directory. They do not install candidates into the package.

## Inputs

| Dataset | Configuration | Required object |
|---|---|---|
| Baboon GSE98965 | `SCP_BABOON_INPUT` | RData containing `baboon_withTOD`, with tissue matrices in `$baboon` and matched ZT vectors in `$tod` |
| Mouse GSE54651 | `SCP_MOUSE_RNASEQ_INPUT` | RDS containing tissue matrices in `$count_clean` and ZT vectors in `$tod` |
| Human GSE160521 | `SCP_STRIATUM_COUNTS`, `SCP_STRIATUM_METADATA` | Gzipped count CSV and the authorized local metadata table described in `02_human_local.R` |
| Human GSE71620 | `SCP_PFC_INPUT` | RDS containing `$expr` and matched `$pheno` |
| Other GEO inputs | `SCP_NATIVE_SAMPLE_PLAN`, `SCP_NATIVE_CACHE` | Reviewed sample/time plan and downloaded or cached GEO expression |

Local metadata and expression files are supplied by the user. These recipes do not grant access or redistribution permission for those inputs.

## Analysis blocks

1. Set the input variables and optionally `SCP_NATIVE_DATASETS` to a comma-separated dataset list.
2. Run `01_baboon_mouse.R`, `02_human_local.R`, or `03_geo_inputs.R` to prepare the selected inputs.
3. Run `04_fit_pilots.R` to fit both recorded harmonic variants and write a comparison against the unchanged package files.
4. Inspect the comparison before considering any later replacement. A traced recipe is not a verified reproduction.

The fitting script uses `estCircadianParam()` or `estCircadianParam2H()`, with paired amplitude and residual noise and seed 2025. The source index preserves the actual file settings. Most K = 1 files use raw p < 0.01 and top 300 genes. Eight human K = 1 files use p < 0.1. Existing K = 2 files use p < 0.01 and top 500 genes. Prepared normalized atlases are not transformed again.

Run `Rscript examples/pilots/native/05_smoke_fixture.R` for a small synthetic check. This creates its own temporary directory and leaves it available for inspection.

## GEO sample selection

The earlier GEO helper selected platforms, time columns, and tissue or condition arms using heuristics. The source index identifies that helper and its hash. `legacy_metadata_selection.R` retains its time and arm selection functions for tracing candidate plans. The numbered scripts do not apply those heuristics automatically. Those selections have not been independently recovered for all 35 generic GEO arms in this update. A reviewed sample plan is therefore required for the new preparation script.

The plan is a CSV with `base_file`, `sample`, `time`, and `eset_index`. Each row identifies one expression column and its input time in hours. `base_file` must match the source index. Add `expression_file` when the selected series matrix is empty and supplementary expression is needed. The expression matrix must use the same sample keys. The preparation script fails on unmatched columns instead of assuming an order.

GSE10644 originally parsed am/pm labels from sample titles. GSE112660 originally stripped nonnumeric characters from collection time and divided by 100. Its exact source format needs confirmation before that convention is adopted in a new sample plan. Record any changed interpretation and compare the resulting fits. Do not claim numerical identity without that check.

Gene symbols come from supplied or GEO platform annotation when available. Unmapped identifiers remain unchanged. Historical annotation mappings have not all been reverified.

Current GTEx builders are in `../gtex_v10/`. The native recipes do not route through an older GTEx builder.
