#' Curate the pilot subset for the shinyapps.io free-tier demo
#'
#' Free tier caps each app at 50 MB. With ~5 MB for the SCP package
#' itself we have ~45 MB for pilot data. Strategy:
#'   - all human (broad tissue + design coverage, ~11 MB)
#'   - all baboon (61 tissues, ~13 MB, the Mure 2018 atlas)
#'   - selected mouse (Zhang 2014 12-organ atlas + Hughes liver + a few
#'     KO/treatment arms, ~10 MB)
#'   - all rat (only 6 pilots, ~6 MB)
#' = ~95 pilots, ~38-40 MB, ~7 MB headroom.
#'
#' Writes scripts/shinyapps_pilots.txt (one pilot file per line, relative
#' to inst/extdata/pilots/) for consumption by build_shinyapps_demo.R.
#'
#' Usage:
#'   Rscript scripts/curate_shinyapps_pilots.R

setwd("/home/qtp1/Projects/Circadian/Kyle/Circadian-analysis-main/R/v1/PowerSim")

m <- read.csv("inst/extdata/pilots/manifest.csv",
              stringsAsFactors = FALSE, check.names = FALSE)
m <- m[m$status == "ingested", , drop = FALSE]

# All human and all rat fit comfortably
keep_human <- m[m$species == "human", , drop = FALSE]
keep_rat   <- m[m$species == "rat",   , drop = FALSE]

# Baboon: keep ALL brain / CNS / eye regions (especially SCN), plus a
# representative selection of peripheral tissues. Mure 2018 codes:
baboon_brain <- c(
  "SCN",  # suprachiasmatic nucleus (the master clock)
  "AMY", "ARC", "BLA", "CER", "COR", "HAB", "HIP", "LH",
  "LGP", "MGP", "MMB", "OLB", "PIN", "PIT", "PON", "PRA",
  "PRC", "PUT", "PVN", "SON", "STF", "SUN", "THA", "VMH",
  "ONH", "RET", "RPE"  # eye / CNS-adjacent
)
baboon_peripheral <- c(
  "LIV", "HEA", "ADC", "ADM", "KIC", "KIM", "LUN", "PAN",
  "THR", "SPL", "WAM", "WAP", "MEL", "TES", "OMF", "SKI",
  "AOR", "MUA"
)
keep_baboon <- m[m$species == "baboon" &
                 m$tissue %in% c(baboon_brain, baboon_peripheral),
                 , drop = FALSE]

# Mouse: select pilots representing distinct tissues/designs without
# over-sampling the GSE54650 atlas
mouse_keep_files <- c(
  # Zhang 2014 atlas: 8 tissues out of 12 (uses 3-letter codes)
  "mouse/GSE54650_LIV_All.rds",
  "mouse/GSE54650_KID_All.rds",
  "mouse/GSE54650_CER_All.rds",
  "mouse/GSE54650_HEA_All.rds",
  "mouse/GSE54650_ADR_All.rds",
  "mouse/GSE54650_HYP_All.rds",
  "mouse/GSE54650_LUN_All.rds",
  "mouse/GSE54650_MUS_All.rds",
  # Hughes 2009 liver (gold-standard 1-hour grid)
  "mouse/GSE11923_Liver_All.rds",
  # Constant darkness liver
  "mouse/GSE11516_Liver_DD.rds",
  # KO arms paired with WT
  "mouse/GSE27366_Kidney_WT.rds",
  "mouse/GSE27366_Kidney_ClockKO.rds",
  "mouse/GSE35026_Adipose_WT.rds",
  "mouse/GSE35026_Adipose_Bmal1KO.rds",
  # Cardiomyopathy heart (Bray 2008)
  "mouse/GSE10045_Heart_VentricleWT.rds",
  "mouse/GSE10045_Heart_VentricleCCMTg.rds"
)
keep_mouse <- m[m$species == "mouse" & m$file %in% mouse_keep_files,
                , drop = FALSE]

curated <- rbind(keep_human, keep_baboon, keep_mouse, keep_rat)

# Compute size budget
sz <- function(file) {
  p <- file.path("inst/extdata/pilots", file)
  if (file.exists(p)) file.info(p)$size else 0
}
curated$size_mb <- sapply(curated$file, sz) / 1048576

cat("Per-species pilot counts and sizes for the shinyapps demo:\n")
for (sp in unique(curated$species)) {
  ms <- curated[curated$species == sp, ]
  cat(sprintf("  %-15s %3d pilots, %5.1f MB\n",
              sp, nrow(ms), sum(ms$size_mb)))
}
cat(sprintf("  TOTAL          %3d pilots, %5.1f MB\n",
            nrow(curated), sum(curated$size_mb)))

cap_mb <- 50
overhead_mb <- 6
budget_mb <- cap_mb - overhead_mb
if (sum(curated$size_mb) > budget_mb) {
  warning(sprintf("Bundle exceeds %.0f MB budget (%.1f MB). Trim mouse subset.",
                  budget_mb, sum(curated$size_mb)))
}

writeLines(curated$file, "scripts/shinyapps_pilots.txt")
cat(sprintf("\nWrote scripts/shinyapps_pilots.txt with %d entries.\n",
            nrow(curated)))
cat("Next: Rscript scripts/build_shinyapps_demo.R\n")
