# Paper figure code and data

The current application figures are in [updated_application](updated_application/README.md).
Run `Rscript paper_reproduction/updated_application/redraw_all.R` from the repository root.

The folders below retain the earlier paper figure inputs and plotting code. Their
`redraw_all.R` entrypoint continues to write the earlier figures to `generated/`.

Run from this folder in a terminal:

```sh
Rscript redraw_all.R
```

The scripts write to `generated/`. They use the supplied panel results and saved figure sources.

| Figure | Plotting data | Redraw script | Figure files |
| --- | --- | --- | --- |
| Figure 1 | `Figure1/data/panel_data.rds` | `Figure1/redraw.R` | `Figure1/final/` |
| Figure 2 | `Figure2/data/panel_data.rds` | `Figure2/redraw.R` | `Figure2/final/` |
| Figure 3 | `Figure3/data/fig3_values.csv` | `Figure3/redraw.R` | `Figure3/final/` |
| Figure 4 | `Figure4/data/` and the vector figure in `Figure4/final/` | `Figure4/redraw.R` | `Figure4/final/` |
| Figure 5 | `Figure5/data/` | `Figure5/redraw.R` | `Figure5/final/` |
| Figure S1 | `Supplement/FigureS1/data/panel_data.rds` | `Supplement/FigureS1/redraw.R` | `Supplement/FigureS1/final/` |
| Tables S1 and S2 | `Supplement/Tables/data/` | `Supplement/Tables/redraw.R` | `Supplement/Tables/code/` |

Each script can also be run separately, for example `Rscript Figure1/redraw.R`. Figure 4 writes the supplied vector figure as `generated/Figure4A.pdf` and redraws B/C as `generated/Figure4BC.pdf`.

Figure 5 retains the `Fig6` filenames used by its source scripts. Figure S1 retains the `Fig5` filenames. The `code/original/` folders hold the corresponding source scripts. The shared plotting settings are in `code/onecol_common.R`.

Saved GTEx v8 pilots are in `data/GTEx_v8/`. The GTEx v10 Liver two-harmonic pilot parameters and sampling distribution are in `data/GTEx_v10/`. The caudate control pilot is in `data/native/`.

Panel data are provided as RDS files, with CSV copies for inspection. For earlier package source, see [SCP v0.4.60](https://github.com/quythien/SCP/tree/v0.4.60).
