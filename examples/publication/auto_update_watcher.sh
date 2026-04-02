#!/usr/bin/env bash
# =======================================================================
# auto_update_watcher.sh — Monitor completion of 03d / 08c / 08e and
#   run all downstream figure + manuscript updates automatically.
#
# Triggers:
#   03d done  → verify bm_boot RDS, copy bm_tradeoff.pdf to paper/figures/
#   08c done  → run 08d_bootstrap_summary.R (08b already done),
#               run post_completion_report.R to extract Seney numbers,
#               copy bootstrap figure to paper/figures/,
#               patch §3.2 table in manuscript
#   08e done  → run post_completion_report.R for Mouse GSE numbers,
#               re-run updated 08d (adds 4th panel), copy figure
#
# Usage (launch in screen before sleeping):
#   screen -S watcher -dm bash examples/publication/auto_update_watcher.sh
# =======================================================================

set -euo pipefail

ROOT="${POWERSIM_ROOT:-/home/qtp1/Projects/Circadian/Kyle/Circadian-analysis-main/R/v1/PowerSim}"
cd "$ROOT"

RUN_TAG_08="20260330_1126"   # tag used by 08b/08c output dir
RUN_TAG_08E="20260401"       # tag used by 08e output dir
OUT_08="${ROOT}/output/08_two_stage_vs_bootstrap_${RUN_TAG_08}"
OUT_08E="${ROOT}/output/08_two_stage_vs_bootstrap_${RUN_TAG_08E}"
PAPER_FIGS="${ROOT}/paper/PowerSim/figures"
LOG="${ROOT}/output/auto_update_watcher.log"
REPORT="${ROOT}/output/auto_update_report.txt"

DONE_03D="${ROOT}/output/03d_mouse_gse_bm_boot.rds"
DONE_08C="${OUT_08}/s3_seney_comparison.rds"
DONE_08E="${OUT_08E}/s4_mouse_gse_comparison.rds"

POLL_INTERVAL=60   # seconds between checks

_log() { echo "[$(date '+%H:%M:%S')] $*" | tee -a "$LOG"; }
_log "=== auto_update_watcher started ==="
_log "Watching: 03d | 08c | 08e"

DONE_FLAG_03D=false
DONE_FLAG_08C=false
DONE_FLAG_08E=false

# -----------------------------------------------------------------------
# Helper: run R script with POWERSIM_ROOT set
# -----------------------------------------------------------------------
_rscript() {
  POWERSIM_ROOT="$ROOT" RUN_TAG="$RUN_TAG_08" Rscript "$@" >> "$LOG" 2>&1
}

# -----------------------------------------------------------------------
# Action: 03d finished — verify + copy bm_tradeoff.pdf
# -----------------------------------------------------------------------
_handle_03d() {
  _log "03d DONE — verifying bootstrap RDS..."
  POWERSIM_ROOT="$ROOT" Rscript --vanilla - >> "$LOG" 2>&1 <<'REOF'
  root <- Sys.getenv("POWERSIM_ROOT")
  setwd(root)
  rds <- "output/03d_mouse_gse_bm_boot.rds"
  bm  <- readRDS(rds)
  stopifnot(is.list(bm), !is.null(bm$power_mean), !is.null(bm$N_values))
  n_B <- length(bm$B_values); n_N <- length(bm$N_values)
  pm  <- matrix(bm$power_mean[,,1], nrow = n_N, ncol = n_B)
  cat(sprintf("03d verified: B=%s  N=%s\n",
    paste(bm$B_values, collapse=","), paste(bm$N_values, collapse=",")))
  for (b in seq_len(n_B)) {
    n80 <- bm$N_values[which(pm[,b] >= 0.80)[1]]
    cat(sprintf("  B=%d: n80=%s  power@N24=%.1f%% N48=%.1f%% N72=%.1f%%\n",
      bm$B_values[b], ifelse(is.na(n80),">max",n80),
      pm[match(24,bm$N_values),b]*100,
      pm[match(48,bm$N_values),b]*100,
      pm[match(72,bm$N_values),b]*100))
  }
REOF

  # bm_tradeoff.pdf is regenerated automatically by 03d_bvsm_mouse_gse_only.R
  # Copy to paper figures
  BM_PDF="${ROOT}/paper/PowerSim/figures/bm_tradeoff.pdf"
  if [ -f "$BM_PDF" ]; then
    _log "bm_tradeoff.pdf already updated by 03d script"
  else
    _log "WARNING: bm_tradeoff.pdf not found at $BM_PDF"
  fi
  _log "03d update complete."
}

# -----------------------------------------------------------------------
# Action: 08c finished — run 08d, extract Seney numbers, patch manuscript
# -----------------------------------------------------------------------
_handle_08c() {
  _log "08c DONE — running 08d_bootstrap_summary.R..."
  POWERSIM_ROOT="$ROOT" RUN_TAG="$RUN_TAG_08" \
    Rscript examples/publication/08d_bootstrap_summary.R >> "$LOG" 2>&1 \
    && _log "08d complete" \
    || _log "WARNING: 08d failed — check $LOG"

  # Copy combined bootstrap figure to paper figures
  S4_PDF="${OUT_08}/s4_ci_width_summary.pdf"
  if [ -f "$S4_PDF" ]; then
    cp "$S4_PDF" "${PAPER_FIGS}/bootstrap_summary.pdf"
    _log "Copied bootstrap_summary.pdf → paper/figures/"
  fi

  # Extract Seney numbers and write manuscript patch
  _log "Extracting Seney bootstrap numbers for manuscript..."
  POWERSIM_ROOT="$ROOT" Rscript --vanilla - >> "$LOG" 2>&1 <<REOF
root   <- Sys.getenv("POWERSIM_ROOT")
setwd(root)
source_dir <- file.path(root, "code")
old_wd <- setwd(source_dir); source("setup.R"); setwd(old_wd)

rds_08c <- "${DONE_08C}"
rds_08b <- "${OUT_08}/s2_d1d2_comparison.rds"
rds_08a <- "${OUT_08}/s1_baboon_comparison.rds"

read_comp <- function(f) {
  if (!file.exists(f)) return(NULL)
  dat <- readRDS(f)
  comp <- dat\$comparison\$comparison
  list(
    n      = comp\$n,
    ts     = comp\$two_stage_power,
    bm     = comp\$boot_power_mean,
    lo     = comp\$boot_ci_lo,
    hi     = comp\$boot_ci_hi,
    n80_ts = dat\$comparison\$n80_two_stage,
    n80_bm = dat\$comparison\$n80_boot_median,
    n80_lo = dat\$comparison\$n80_boot_lo,
    n80_hi = dat\$comparison\$n80_boot_hi
  )
}

seney  <- read_comp(rds_08c)
d1d2   <- read_comp(rds_08b)
baboon <- read_comp(rds_08a)

lines <- c(
  "=== AUTO-GENERATED MANUSCRIPT NUMBERS ===",
  sprintf("Generated: %s", Sys.time()), "",
  "--- SENEY (MDD vs CTL ACC, passive, n=60) ---"
)
if (!is.null(seney)) {
  lines <- c(lines,
    sprintf("  Two-stage n80: %s", ifelse(is.na(seney\$n80_ts), ">max(N)", seney\$n80_ts)),
    sprintf("  Bootstrap n80: %s  [95%% CI: %s, %s]",
      ifelse(is.na(seney\$n80_bm), ">max(N)", round(seney\$n80_bm)),
      ifelse(is.na(seney\$n80_lo), "NA", seney\$n80_lo),
      ifelse(is.na(seney\$n80_hi), "NA", seney\$n80_hi)),
    "  Power table (FDR 5%):  N | Two-stage | Bootstrap [95% CI]"
  )
  for (i in seq_along(seney\$n)) {
    lines <- c(lines, sprintf("    N=%3d  | %5.1f%%  | %5.1f%% [%5.1f%%, %5.1f%%]",
      seney\$n[i], 100*seney\$ts[i], 100*seney\$bm[i],
      100*seney\$lo[i], 100*seney\$hi[i]))
  }
}
lines <- c(lines, "",
  "--- D1D2 (D1 vs D2 striatum, active, n=45) ---")
if (!is.null(d1d2)) {
  lines <- c(lines,
    sprintf("  Two-stage n80: %s", ifelse(is.na(d1d2\$n80_ts), ">max(N)", d1d2\$n80_ts)),
    sprintf("  Bootstrap n80: %s  [95%% CI: %s, %s]",
      ifelse(is.na(d1d2\$n80_bm), ">max(N)", round(d1d2\$n80_bm)),
      ifelse(is.na(d1d2\$n80_lo), "NA", d1d2\$n80_lo),
      ifelse(is.na(d1d2\$n80_hi), "NA", d1d2\$n80_hi))
  )
  for (i in seq_along(d1d2\$n)) {
    lines <- c(lines, sprintf("    N=%3d  | %5.1f%%  | %5.1f%% [%5.1f%%, %5.1f%%]",
      d1d2\$n[i], 100*d1d2\$ts[i], 100*d1d2\$bm[i],
      100*d1d2\$lo[i], 100*d1d2\$hi[i]))
  }
}
lines <- c(lines, "",
  "=== MANUSCRIPT PATCH — §3.2 TABLE ROWS TO ADD/REPLACE ===",
  "Copy these rows into tab:bootstrap_results in PowerSim_Paper2.tex:",
  ""
)
if (!is.null(d1d2)) {
  for (i in seq_along(d1d2\$n)) {
    lines <- c(lines, sprintf(
      "%d & %.1f\\%% & %.1f\\%% [%.1f\\%%, %.1f\\%%] & D1D2 \\\\",
      d1d2\$n[i], 100*d1d2\$ts[i], 100*d1d2\$bm[i],
      100*d1d2\$lo[i], 100*d1d2\$hi[i]))
  }
}
if (!is.null(seney)) {
  lines <- c(lines, "")
  for (i in seq_along(seney\$n)) {
    lines <- c(lines, sprintf(
      "%d & %.1f\\%% & %.1f\\%% [%.1f\\%%, %.1f\\%%] & Seney \\\\",
      seney\$n[i], 100*seney\$ts[i], 100*seney\$bm[i],
      100*seney\$lo[i], 100*seney\$hi[i]))
  }
}
writeLines(lines, "${REPORT}")
cat(paste(lines, collapse="\n"), "\n")
REOF
  _log "Seney/D1D2 numbers written to: $REPORT"
  _log "08c update complete."
}

# -----------------------------------------------------------------------
# Action: 08e finished — extract Mouse GSE numbers, re-run 08d (4-panel)
# -----------------------------------------------------------------------
_handle_08e() {
  _log "08e DONE — extracting Mouse GSE bootstrap numbers..."
  POWERSIM_ROOT="$ROOT" Rscript --vanilla - >> "$LOG" 2>&1 <<REOF
root   <- Sys.getenv("POWERSIM_ROOT")
setwd(root)
source_dir <- file.path(root, "code")
old_wd <- setwd(source_dir); source("setup.R"); setwd(old_wd)

rds <- "${DONE_08E}"
dat  <- readRDS(rds)
comp <- dat\$comparison\$comparison
report <- "${REPORT}"
existing <- if (file.exists(report)) readLines(report) else character(0)

lines <- c("", "--- MOUSE LIV vs CER (GSE54651, active, n=8) ---",
  sprintf("  Two-stage n80: %s", ifelse(is.na(dat\$comparison\$n80_two_stage), ">max(N)", dat\$comparison\$n80_two_stage)),
  sprintf("  Bootstrap n80: %s  [95%% CI: %s, %s]",
    ifelse(is.na(dat\$comparison\$n80_boot_median), ">max(N)", round(dat\$comparison\$n80_boot_median)),
    ifelse(is.na(dat\$comparison\$n80_boot_lo), "NA", dat\$comparison\$n80_boot_lo),
    ifelse(is.na(dat\$comparison\$n80_boot_hi), "NA", dat\$comparison\$n80_boot_hi)),
  "  Power table:",
  "  N  | Two-stage | Bootstrap [95% CI]"
)
for (i in seq_along(comp\$n)) {
  lines <- c(lines, sprintf("  N=%3d  | %5.1f%%  | %5.1f%% [%5.1f%%, %5.1f%%]",
    comp\$n[i], 100*comp\$two_stage_power[i], 100*comp\$boot_power_mean[i],
    100*comp\$boot_ci_lo[i], 100*comp\$boot_ci_hi[i]))
}
writeLines(c(existing, lines), report)
cat(paste(lines, collapse="\n"), "\n")
REOF

  # Copy per-dataset PDF
  S4E_PDF="${OUT_08E}/s4_mouse_gse_comparison.pdf"
  [ -f "$S4E_PDF" ] && cp "$S4E_PDF" "${PAPER_FIGS}/bootstrap_mouse_gse.pdf" \
    && _log "Copied bootstrap_mouse_gse.pdf → paper/figures/"

  # Re-run 08d if 08c is also done (to include 4th panel)
  if [ "$DONE_FLAG_08C" = true ]; then
    _log "Re-running 08d with 4-panel (Baboon + D1D2 + Seney + Mouse GSE)..."
    POWERSIM_ROOT="$ROOT" RUN_TAG="$RUN_TAG_08" \
      Rscript examples/publication/08d_bootstrap_summary.R >> "$LOG" 2>&1 \
      && _log "08d 4-panel complete" \
      || _log "WARNING: 08d re-run failed"
    S4_PDF="${OUT_08}/s4_ci_width_summary.pdf"
    [ -f "$S4_PDF" ] && cp "$S4_PDF" "${PAPER_FIGS}/bootstrap_summary.pdf" \
      && _log "Updated bootstrap_summary.pdf in paper/figures/"
  else
    _log "08c not yet done — will re-run 08d once 08c also completes"
  fi

  _log "08e update complete."
}

# -----------------------------------------------------------------------
# Main poll loop
# -----------------------------------------------------------------------
while true; do
  # Check 03d
  if [ "$DONE_FLAG_03D" = false ] && [ -f "$DONE_03D" ]; then
    DONE_FLAG_03D=true
    _handle_03d
  fi

  # Check 08c
  if [ "$DONE_FLAG_08C" = false ] && [ -f "$DONE_08C" ]; then
    DONE_FLAG_08C=true
    _handle_08c
    # If 08e was already done, re-run 08d now too
    if [ "$DONE_FLAG_08E" = true ]; then
      _log "08e was already done — re-running 08d for 4-panel..."
      POWERSIM_ROOT="$ROOT" RUN_TAG="$RUN_TAG_08" \
        Rscript examples/publication/08d_bootstrap_summary.R >> "$LOG" 2>&1
    fi
  fi

  # Check 08e
  if [ "$DONE_FLAG_08E" = false ] && [ -f "$DONE_08E" ]; then
    DONE_FLAG_08E=true
    _handle_08e
  fi

  # Exit once all three are done
  if [ "$DONE_FLAG_03D" = true ] && [ "$DONE_FLAG_08C" = true ] && [ "$DONE_FLAG_08E" = true ]; then
    _log "=== All jobs complete. Watcher exiting. ==="
    _log "Key outputs:"
    _log "  paper/figures/bm_tradeoff.pdf       — updated B vs m (Mouse+Baboon+D1D2)"
    _log "  paper/figures/bootstrap_summary.pdf — 3-panel (or 4-panel) bootstrap vs two-stage"
    _log "  output/auto_update_report.txt        — all n80 numbers for manuscript"
    _log "  Check report for §3.2 LaTeX table rows to copy into PowerSim_Paper2.tex"
    break
  fi

  sleep "$POLL_INTERVAL"
done
