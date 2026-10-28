# Pilot builders

Run the scripts from the SCP repository root. Each analysis has numbered blocks
for inputs, preparation, fitting, and output checks. Candidates are written to a
separate staging directory before any package files are replaced.

## GTEx v10

The active GTEx pilots use v10. The two input-time examples are CHIRAL-estimated
time (DIP) and season-adjusted time of death (TOD). The app uses DIP.

The [v10 recipe](gtex_v10/README.md) covers public-source normalization,
CHIRAL inference, paired pilot fits, diagnostics, sample-size comparisons,
and the demo bundle. Input versions, hashes, settings, and output checks are
recorded with the build. Source expression and donor-linked recorded times are
not placed in the code repository.

## Other datasets

The [native-time recipes](native/README.md) and
[source index](builder_source_index.csv) trace the unchanged non-GTEx pilots.
They record the original helper hashes and the settings stored in each pilot.
The new staging workflow has a synthetic check. It has not been used to rebuild
the native database during this GTEx update. The index records remaining sample
selection and annotation checks explicitly.

## Earlier versions

Earlier package commits preserve the previous inputs and code. The
[earlier GTEx builders](https://github.com/quythien/SCP/tree/418728f05bf32ba5c1150549f93016486b8af68f/examples/pilots)
remain at their original revision. Local paper inputs and figures remain
unchanged. They are separate from the current v10 staging workflow.

## References

- Talamanca L, Gobet C, Naef F. 2023. Sex-dimorphic and age-dependent organization of 24-hour gene expression rhythms in humans. Science 379:478-483. https://doi.org/10.1126/science.add0846
- GTEx Consortium. 2020. The GTEx Consortium atlas of genetic regulatory effects across human tissues. Science 369:1318-1330. https://doi.org/10.1126/science.aaz1776

The manuscript identifies the GTEx study as dbGaP accession phs000424. Public
v10 expression and QC file names and hashes are recorded separately so study
accession and expression release remain clear.
