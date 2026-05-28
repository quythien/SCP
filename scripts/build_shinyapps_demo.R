#' Build the slim shinyapps.io deploy bundle for SCP
#'
#' shinyapps.io free tier caps each app at 50 MB. The full pilot database
#' (~80 MB across 127 rds) exceeds that, so the public live demo ships a
#' curated subset. Power users install the full SCP package from GitHub
#' to get all 127 pilots.
#'
#' Strategy:
#'   1. Read the curated pilot list at scripts/shinyapps_pilots.txt
#'   2. Copy app.R and just-those-pilots into inst/shiny_demo/
#'   3. Sub-set the manifest to match.
#'   4. Report the resulting bundle size.
#'
#' Deploy with:
#'   rsconnect::deployApp("inst/shiny_demo/", appName = "SCP")
#'
#' Usage:
#'   Rscript scripts/build_shinyapps_demo.R

setwd("/home/qtp1/Projects/Circadian/Kyle/Circadian-analysis-main/R/v1/PowerSim")

# ---------- inputs ----------------------------------------------------------
curated_list  <- "scripts/shinyapps_pilots.txt"
pilots_src    <- "inst/extdata/pilots"
app_src       <- "inst/shiny/app.R"
demo_dir      <- "inst/shiny_demo"
pilots_dst    <- file.path(demo_dir, "pilots")

# ---------- read curated list -----------------------------------------------
if (!file.exists(curated_list)) {
  stop(sprintf("Curated pilot list not found at %s. Create it first (one ",
               "pilot file path per line, relative to inst/extdata/pilots/).",
               curated_list))
}
chosen <- readLines(curated_list)
chosen <- trimws(chosen)
chosen <- chosen[nzchar(chosen) & !startsWith(chosen, "#")]
cat(sprintf("Curated list: %d pilots\n", length(chosen)))

# ---------- prepare destination ---------------------------------------------
if (dir.exists(demo_dir)) unlink(demo_dir, recursive = TRUE)
dir.create(pilots_dst, recursive = TRUE)
for (sub in c("human", "baboon", "other_primate", "mouse", "rat")) {
  dir.create(file.path(pilots_dst, sub), showWarnings = FALSE)
}

# ---------- copy app.R + slim manifest + chosen rds ------------------------
file.copy(app_src, file.path(demo_dir, "app.R"), overwrite = TRUE)

m <- read.csv(file.path(pilots_src, "manifest.csv"),
              stringsAsFactors = FALSE, check.names = FALSE)
m_slim <- m[m$file %in% chosen, , drop = FALSE]
m_slim <- m_slim[m_slim$status == "ingested", , drop = FALSE]
write.csv(m_slim, file.path(pilots_dst, "manifest.csv"), row.names = FALSE)
cat(sprintf("Slim manifest: %d rows (of %d total ingested)\n",
            nrow(m_slim), sum(m$status == "ingested")))

copied <- 0L; missing <- character(0)
for (rel in chosen) {
  src <- file.path(pilots_src, rel)
  dst <- file.path(pilots_dst, rel)
  if (!file.exists(src)) { missing <- c(missing, rel); next }
  ok <- file.copy(src, dst, overwrite = TRUE)
  if (ok) copied <- copied + 1L
}
cat(sprintf("Copied %d / %d pilot rds.\n", copied, length(chosen)))
if (length(missing) > 0) {
  cat("WARNING: missing source rds (will be excluded from demo):\n")
  for (x in missing) cat("  ", x, "\n")
}

# ---------- patch app.R for demo-local pilot path ---------------------------
app_text <- readLines(file.path(demo_dir, "app.R"))
# tell the app it is in shinyapps demo mode (loads ./pilots/manifest.csv)
app_text <- c('Sys.setenv(SCP_SHINY_DEMO = "1")', app_text)
writeLines(app_text, file.path(demo_dir, "app.R"))

# ---------- size report -----------------------------------------------------
sizes <- file.info(list.files(demo_dir, recursive = TRUE,
                              full.names = TRUE))$size
total_mb <- sum(sizes, na.rm = TRUE) / 1048576
cat(sprintf("\nBundle total: %.1f MB across %d files.\n",
            total_mb, length(sizes)))
if (total_mb > 50) {
  warning(sprintf("Bundle exceeds shinyapps.io free-tier 50 MB cap (%.1f MB).",
                  total_mb))
} else {
  cat(sprintf("OK: under shinyapps.io free-tier 50 MB cap (%.1f MB headroom).\n",
              50 - total_mb))
}

cat("\nNext step:\n")
cat("  rsconnect::deployApp(\"inst/shiny_demo/\", appName = \"SCP\")\n")
