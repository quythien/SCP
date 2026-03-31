# Dataset Survey — Circadian DR Power Analysis
## Date: 2026-03-19

---

## Overview

This document summarizes all datasets explored for DR (Differential Rhythmicity) power
analysis, their circadian signal strength, and which tissue pairs are good candidates
for power study examples beyond the human aging reference case.

**Key metrics:**
- `r = A/sigma` — effect size (amplitude / noise); drives power most
- `prop_rhythmic` — fraction of genes rhythmic at p < 0.05
- `prop_DR` — fraction of genes rhythmic in one group but not the other
- `design` — active (controlled ZT) vs passive (post-mortem irregular TOD)

---

## 1. Human Aging Brain — BA11/BA47 (REFERENCE CASE)

**File:** `combined_data.rds`
(`Projects/Collaborative/Paper/Congruence/PNAS_aging/data/`)

**Structure:** 19,998 genes × 136 samples (62 young, 74 old)
- Regions: BA11 + BA47 (prefrontal cortex)
- Age groups: younger (16–45 yr) and older (60–96 yr)
- Design: **passive** — irregular post-mortem TOD, range [-5.8, 17.4] h

**Signal (all samples pooled, younger as pilot):**
- prop_rhythmic: ~20% at p < 0.10 (from `03_power_core.R`)
- r_median: ~0.56 (moderate)
- Amplitude: 0.126 ± 0.065

**DR — Young vs Old (data-driven, p < 0.10):**
- Young: prop_rhy = 5.4%, r_median = 0.51
- Old:   prop_rhy = 9.8%
- prop_DR = 13.7% (discordant rhythmicity)

**DR Power (young vs old, passive design, n per group):**
| N per group | Power |
|---|---|
| 30  | 0.5%  |
| 60  | 7.7%  |
| 100 | 20.0% |
| 150 | 41.1% |
| 200 | 53.1% |
| 300 | 70.8% |

**Conclusion:** Human aging DR is **low power** even at large N (n80 >> 300). The
weak signal (r ~0.5) and moderate prop_DR (~14%) combine to make this a challenging
detection problem. The aging comparison is biologically important but statistically
underpowered at realistic sample sizes.

**Usage in `03_power_core.R`:** Parameters estimated from younger subjects only;
prop_DR = 0.15 assumed hypothetically (not estimated from young vs old comparison).

---

## 2. Baboon Multi-Tissue (STRONG SIGNAL — RECOMMENDED)

**File:** `CAMO_PRC_hmb.RData` (object: `baboon_withTOD`)
(`Projects/Collaborative/GTEXdata/data/`)

**Structure:** 4,938 genes (human orthologs) × 12 samples per tissue
- 61 tissues total
- Design: **active** — controlled ZT, time points: 0, 2, 4, ..., 22 h (n=1 per ZT)
- Data: raw CPM → log2(CPM+1) transform needed
- Gene IDs: human ENSEMBL (orthologs shared with GTEx)

**Note:** TOD is a named list (per tissue), not a single vector. Each tissue has
its own `baboon_withTOD$tod[[tissue]]` vector (always 0:2:22).

**Top tissues by effect size r (scan of all 61 tissues, p < 0.05):**
| Tissue | prop_rhythmic | r_median | Notes |
|---|---|---|---|
| THR (thyroid)       | 58% | 2.08 | Highest effect size |
| STF (stomach fundus)| 54% | 1.92 | |
| ILE (ileum)         | 54% | 1.95 | |
| HEA (heart)         | 61% | 1.84 | Highest prop_rhythmic |
| ADM (adrenal medulla)| 46% | 1.88 | |
| PVN                 | 41% | 1.87 | |
| LUN (lung)          | 38% | 1.71 | |
| LIV (liver)         | 22% | 1.63 | |
| CER (cerebellum)    | 7%  | 1.57 | |

**Best DR pairs — all tissues (top by delta prop_rhythmic):**
| Pair | prop_rhy (t1/t2) | prop_DR | Expected n80 |
|---|---|---|---|
| HEA vs SMM | 61% / 10% | 53% | ~36 |
| HEA vs HIP | 61% / 12% | 51% | ~36 |
| THR vs SMM | 58% / 10% | ~48% | ~36 |

**Candidate pairs for publication (biologically motivated):**
| Pair | prop_rhy (t1/t2) | prop_DR | Notes |
|---|---|---|---|
| LUN vs CER | 38% / 7%  | **41%** | Lung vs brain — large diff |
| LUN vs LIV | 38% / 22% | **32%** | Both peripheral, meaningful |
| LIV vs CER | 22% / 7%  | **27%** | Liver vs brain |

**Recommendation:** LUN vs CER or LUN vs LIV for publication — these are
biologically interpretable pairs that span peripheral (lung/liver) vs CNS (cerebellum),
mirroring the mouse comparison below. Prop_DR 27–41% with r ~1.6–1.7 → n80 likely
in the 24–48 range.

---

## 3. Mouse Multi-Tissue GSE54651 — Zhang et al. 2014 (STRONG SIGNAL)

**File:** `mice_GSE54651_CPM.RData`
(`Projects/Collaborative/GTEXdata/data/`)

**Load:** `readRDS()` (not `load()` — file has non-standard magic number)

**Structure:** 12 tissues, 8 samples per tissue
- Design: **active** — Circadian Time (CT) points: 22, 28, 34, 40, 46, 52, 58, 64 h
  (every 6h over 2 days; folds to ZT 22, 4, 10, 16 × 2 = 4 unique ZT, 2 replicates each)
- Data: already log2-scaled (range −5.96 to 15.93), use `count_clean`
- Gene IDs: **ENSMUSG** (mouse Ensembl) — NOT human orthologs
- n=1 per CT time point (n=2 per ZT when folded mod 24)

**Tissues:** ADR, AOR, BFAT, BSTM, CER, HEA, HYP, KIC, LIV, LUN, MUS, WFAT

**Signal by tissue (p < 0.05, subset 2000 genes):**
| Tissue | prop_rhythmic | r_median | Notes |
|---|---|---|---|
| LIV (liver)    | 25% | 2.89 | **Strongest — well-known circadian tissue** |
| KIC (kidney)   | 23% | 2.92 | |
| HEA (heart)    | 21% | 2.93 | |
| MUS (muscle)   | 15% | 2.77 | |
| ADR (adrenal)  | 14% | 2.82 | |
| AOR (aorta)    | 14% | 2.76 | |
| LUN (lung)     | 12% | 2.80 | |
| CER (cerebellum)| 10% | 2.73 | |
| HYP (hypothalamus)| 8% | 2.58 | |
| BFAT           | 0%  | — | Too few samples (n=8) to detect |
| BSTM           | 0%  | — | Too few samples |
| WFAT           | 0%  | — | Too few samples |

**Best DR pairs:**
| Pair | prop_rhy (t1/t2) | Approx prop_DR | Notes |
|---|---|---|---|
| LIV vs CER | 25% / 10% | ~20–25% | **Primary recommendation** |
| LIV vs HYP | 25% / 8%  | ~22–25% | Liver vs hypothalamus |
| LUN vs CER | 12% / 10% | ~10–15% | Lung vs cerebellum |
| KIC vs CER | 23% / 10% | ~18–22% | |

**Key advantage:** r_median = 2.7–3.4 (much stronger than baboon r ~1.7 and human
r ~0.56). Even with modest N (likely 12–24 per group), expect high power.

**Limitation:** Only n=8 samples per tissue (n=1 per CT point). Parameter estimates
from pilot data will have high uncertainty. For power simulation, this means CI bands
from bootstrap will be wide.

**Recommendation:** **LIV vs CER** as primary mouse pair — mirrors the well-known
finding that liver is the most rhythmic peripheral tissue while cerebellum has weak
circadian expression. LUN vs CER is also valid for comparing with baboon.

---

## 4. GTEx Human (WEAK — NOT RECOMMENDED FOR DR)

**File:** `CAMO_PRC_hmb.RData` (object: `gtex`)
**Structure:** PRC1 (n=209) and PRC2 (n=255)
- Design: **passive** — irregular post-mortem TOD
- Gene IDs: human ENSEMBL (same 4,938 as baboon)

**Signal:**
- GTEx PRC1: prop_rhythmic ~0%, r ~NA (no detectable signal)
- GTEx PRC2: prop_rhythmic ~0%, r ~NA

**Conclusion:** Post-mortem timing noise overwhelms circadian signal. GTEx data
is not suitable for DR analysis in this framework.

---

## 5. Mouse D1 vs D2 Cell Types (MODERATE SIGNAL — LOCAL DATA)

**Files:** `data/mouse_clinicalinfo_03082021_rmOutliers.csv` + `data/mouse_D1D2_logCPMfiltered_counts.csv`

**Structure:** 14,355 genes × 91 samples (45 D1, 46 D2)
- D1 vs D2 = two cell populations (likely striatal dopamine receptor-expressing neurons)
- Design: **active** — controlled ZT, 6 time points: 2, 6, 10, 14, 18, 22h (every 4h)
- ~7–8 biological replicates per ZT per group
- Expression: raw counts → apply `log2(CPM + 1)` normalization (see note below)
- Gene IDs: ENSMUSG

**Note on expression file:** Despite "logCPMfiltered" in the filename, values are raw counts
(range 0–662,435). The filter was applied when selecting genes, not to the count values.

**Signal (p < 0.05, log2 CPM+1):**
- D1: prop_rhythmic = ~10–13%, r_median = ~0.66
- D2: prop_rhythmic = ~13%, r_median = ~0.64
- prop_DR (D1 vs D2) = ~18–22%

**B vs m tradeoff:** B=6 (full 4h coverage) consistently outperforms B=3 (every 8h)
for the same N. More time-point coverage > more replicates at r~0.65. Expected n80: ~80–150.

**Comparison with other datasets:** Signal is between baboon (r~1.7) and human aging (r~0.5).
Uniquely has true biological replicates per ZT (most datasets have n=1 per ZT).

---

## 6. Seney et al. — MDD vs Control ACC (PASSIVE — RECOMMENDED FOR PASSIVE EXAMPLE)

**Files:** `data/ACC_RNA_filtered_normalized.csv` + `data/MD5_MetaData_1-15-25.xlsx` + `data/TOD.xlsx`

**Structure:** 14,455 genes × 120 samples (60 Control, 60 MDD)
- Region: ACC (anterior cingulate cortex)
- Groups: Control (Disease=1) vs MDD (Disease=2)
- Sex: 30 Male + 30 Female per disease group (GENDER: 1=Male, 2=Female)
- Design: **Passive** — irregular post-mortem TOD from OFFC_TIME (official pronounced time)
- Expression: already log-normalized (range −5.22 to 16.08); do NOT re-normalize
- TOD extraction: `as.numeric(format(OFFC_TIME, "%H")) + as.numeric(format(OFFC_TIME, "%M"))/60`
- Sample ID matching: strip trailing letter from colnames: `gsub("[A-Za-z]+$", "", colnames(expr))`

**Signal (p < 0.05, combined n=60/group):**
- Control: prop_rhythmic = ~9%, r_median = ~0.55
- MDD: prop_rhythmic = ~12%, r_median = ~0.54
- prop_DR (xor, via estCircadianParamTwoGroup) = ~19%
- n80 estimate: **~280–300**

**Sex-stratified signal (p < 0.05):**
| Subgroup | n/group | r_median (CTL) | prop_DR | n80 estimate |
|---|---|---|---|---|
| Combined | 60 | 0.55 | ~19% | ~280–300 |
| Male-only | 30 | 0.82 | ~19% | ~200–220 |
| Female-only | 30 | 0.76 | ~15% | ~200–220 |

**Why sex stratification improves power:**
Sex is a well-documented source of heterogeneity in circadian gene expression in the ACC and in MDD biology (Seney et al.). When male and female subjects are pooled, two rhythm profiles with different phase alignments and amplitudes are averaged into a single cosinor estimate. This sex-driven heterogeneity:
- Reduces apparent amplitude A — phase incoherence between sexes causes partial cancellation of the rhythmic signal
- Inflates residual variance σ — cosinor residuals absorb between-sex variation
Together these reduce r = A/σ from ~0.82 (male-only) to ~0.55 (pooled), a ~50% drop. Since power scales as r² × N, a 50% drop in r requires ~4× the N to recover the same power. Stratifying by sex restores within-sex phase coherence, yielding n80 ~200 (male-only) vs ~300 (combined). This is a practical recommendation: controlling a known biological confounder in a passive study can substantially reduce n80, even though the timing heterogeneity cannot be fixed.

**Scripts:**
- `examples/exploratory/05_multi_dataset_DR.R` — Section 3 (combined Seney as passive example)
- `examples/exploratory/07_seney_sex_DR.R` — sex-stratified comparison (combined vs male vs female)

**Data loading:**
```r
library(readxl)
meta_s  <- read_excel("data/MD5_MetaData_1-15-25.xlsx")
tod_s   <- read_excel("data/TOD.xlsx")
expr_s  <- as.matrix(read.csv("data/ACC_RNA_filtered_normalized.csv", row.names=1, check.names=FALSE))
col_ids <- gsub("[A-Za-z]+$", "", colnames(expr_s))
meta_idx <- match(col_ids, as.character(meta_s$HU_NUM))
tod_idx  <- match(col_ids, as.character(tod_s$HU_NUM))
tod_hour <- as.numeric(format(tod_s$OFFC_TIME[tod_idx], "%H")) +
            as.numeric(format(tod_s$OFFC_TIME[tod_idx], "%M")) / 60
disease <- meta_s$Disease[meta_idx]  # 1=Control, 2=MDD
gender  <- meta_s$GENDER[meta_idx]   # 1=Male, 2=Female
```

**Comparison with human aging:**
Both are passive designs with r~0.5 and n80 >300 when combined. Seney is preferred as the paper's passive example because:
1. Sex stratification adds a concrete design lesson (controlling confounders reduces n80 by ~30%)
2. MDD + circadian disruption is a highly active research area
Human aging is retained as a commented-out alternative in `05_multi_dataset_DR.R`.

---

## 7. Summary Comparison

| Dataset | Design | r_median | prop_DR range | n80 estimate | Recommended? |
|---|---|---|---|---|---|
| Mouse GSE54651 LIV vs CER | Active CT | 2.9 | ~20–25% | **< 24** | YES — primary |
| Baboon LUN vs CER | Active ZT | ~1.6 | ~41% | ~24–36 | YES |
| Baboon LUN vs LIV | Active ZT | ~1.7 | ~32% | ~36–48 | YES |
| Baboon HEA vs SMM | Active ZT | ~1.8 | ~53% | ~36 | YES (exploratory) |
| Mouse D1 vs D2 (local) | Active ZT | ~0.66 | ~18–22% | ~80–150 | YES — unique multi-replicate |
| Seney MDD vs Control (ACC) | Passive | ~0.55 | ~19% | ~280–300 | YES — primary passive |
| Seney male-only MDD vs CTL | Passive | ~0.82 | ~19% | ~200–220 | YES — stratified passive |
| Human aging (young vs old) | Passive | ~0.50 | ~14% | >300 | Alt (see 05_multi_dataset_DR.R) |
| GTEx PRC1 vs PRC2 | Passive | ~0 | ~0% | N/A | NO |

---

## 7. Recommended Publication Examples

### Primary examples (progression from strong to moderate signal):

1. **Mouse cross-tissue: LIV vs CER** — highest effect size (r ~2.9), clean active
   design, well-known biological contrast (Zhang 2014). Expected n80 < 24.

2. **Baboon cross-tissue: LUN vs CER** — strong effect (r ~1.6), active ZT design,
   parallels mouse comparison (both lung/peripheral vs cerebellum/CNS). Expected n80 ~24–36.

3. **Mouse D1 vs D2 cell types** — moderate effect (r ~0.66), multi-replicate active
   design. Demonstrates B vs m tradeoff with real data. Expected n80 ~80–150.

4. **Human aging: young vs old brain** — weakest signal (r ~0.5), passive design,
   biologically motivated but underpowered. Illustrates the challenge of passive
   human post-mortem studies. n80 >> 300.

### Script locations:
- **Unified 3-dataset DR comparison**: `examples/exploratory/05_multi_dataset_DR.R` (mouse GSE54651 / baboon / human)
- **Mouse D1 vs D2 DR power**: `examples/exploratory/06_mouse_D1D2_DR.R` ← new
- Baboon tissue scan (all 61 tissues): `examples/exploratory/baboon_tissue_scan.R`
- DR cross-species/tissue (baboon vs GTEx): `examples/exploratory/dr_camo_power.R`
- Core publication pipeline: `examples/publication/03_power_core.R` (younger pilot only)

---

## 7. Data Notes

### Baboon data loading
```r
load("CAMO_PRC_hmb.RData")   # loads: baboon_withTOD, gtex, mice
bab_expr  <- baboon_withTOD$baboon   # named list of 61 tissues (data.frame, raw CPM)
bab_tod   <- baboon_withTOD$tod      # named list of TOD vectors per tissue
# Convert to matrix:  data.matrix(log2(mat + 1))
```

### Mouse data loading
```r
dat <- readRDS("mice_GSE54651_CPM.RData")  # must use readRDS, not load()
mat <- data.matrix(dat$count_clean[[tissue]])  # already log2-scaled; MUST wrap data.matrix()
tod <- dat$tod[[tissue]]           # CT values 22,28,...,64 (use raw CT with period=24)
# Note: gene IDs are ENSMUSG (mouse), not human ENSG — within-mouse DR only
# CT design: n=8 per tissue (1 animal per CT), folds to 4 unique ZT × 2 replicates
# DO NOT analyze day1 and day2 separately — only 1 residual df each; use all 8 CT
```

### Mouse D1 vs D2 loading (local data/)
```r
pheno    <- read.csv("data/mouse_clinicalinfo_03082021_rmOutliers.csv", row.names=1)
expr_raw <- as.matrix(read.csv("data/mouse_D1D2_logCPMfiltered_counts.csv", row.names=1))
# IMPORTANT: values are RAW COUNTS despite "logCPMfiltered" in filename
# Always normalize: log2(CPM + 1)
lib_size <- colSums(expr_raw)
log_mat  <- log2(sweep(expr_raw, 2, lib_size, "/") * 1e6 + 1)
# Key pheno columns: sample, time (ZT: 2,6,10,14,18,22), sex, cell (D1/D2)
d1_idx <- pheno$cell == "D1"   # n=45
d2_idx <- pheno$cell == "D2"   # n=46
# Replicates per ZT: D1 ~7-8/ZT, D2 ~7-8/ZT (true biological replicates)
```

### Human aging loading
```r
dat   <- readRDS("combined_data.rds")
expr  <- dat$expr    # 19998 x 136 matrix (genes x samples)
pheno <- dat$pheno   # 136 x 13 data.frame
# Key columns: AgeGroup ("younger"/"older"), TOD.y, region ("BA11"/"BA47")
```
