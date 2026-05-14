#' =======================================================================
#' suppFig_clock_gene_profiles.R
#'
#' Supplementary gene-profile panel: raw expression vs time for four
#' core clock genes in baboon PAN (GSE98965, Mure et al. 2018), with
#' cosinor and FMM fits overlaid. Visual demonstration that cosinor
#' under-fits the sharp peaks/troughs that FMM captures.
#'
#' Genes chosen by largest amplitude / clearest rhythm in PAN:
#'   NR1D1 (REV-ERBalpha), PER1, ARNTL (BMAL1), DBP.
#'
#' Output: output/main_figures/SuppFig_clock_gene_profiles.pdf
#' =======================================================================

setwd("/home/qtp1/Projects/Circadian/Kyle/Circadian-analysis-main/R/v1/PowerSim")
old_wd <- setwd("code"); source("setup.R"); setwd(old_wd)
source("examples/publication/_pub_style.R")

if (!requireNamespace("FMM", quietly = TRUE))
  stop("FMM package required")

# ---- Load PAN data ----
load("data/CAMO_PRC_hmb.RData")
m_df <- baboon_withTOD$baboon$PAN
m    <- as.matrix(m_df); storage.mode(m) <- "numeric"
tods <- baboon_withTOD$tod$PAN %% 24
ord  <- order(tods)
m    <- m[, ord]
tods <- tods[ord]

# ---- Genes to plot ----
clock_genes <- list(
  ARNTL = "ENSG00000133794",   # BMAL1
  NR1D1 = "ENSG00000126368",   # REV-ERBalpha
  PER1  = "ENSG00000179094",
  DBP   = "ENSG00000105516"
)

# ---- Cosinor OLS fit ----
fit_cosinor <- function(y, t_h, period = 24) {
  w <- 2 * pi / period
  X <- cbind(1, cos(w * t_h), sin(w * t_h))
  b <- solve(crossprod(X), crossprod(X, y))
  M <- b[1]
  A <- sqrt(b[2]^2 + b[3]^2)
  phi <- atan2(b[3], b[2])           # acrophase in radians
  list(M = M, A = A, phi = phi,
       predict = function(tt) M + A * cos(w * tt - phi))
}

# ---- FMM fit (Rueda 2019 nonlinear) ----
fit_fmm <- function(y, t_h, period = 24) {
  t_rad <- t_h * (2 * pi / period)
  yc    <- y - mean(y)
  fit   <- tryCatch(
    FMM::fitFMM(yc, timePoints = t_rad,
                lengthAlphaGrid = 24, lengthOmegaGrid = 12,
                showProgress = FALSE),
    error = function(e) NULL
  )
  if (is.null(fit)) return(NULL)
  M  <- fit@M
  Av <- fit@A
  al <- fit@alpha
  be <- fit@beta
  om <- fit@omega
  list(M = M, A = Av, alpha = al, beta = be, omega = om,
       predict = function(tt) {
         tt_r <- tt * (2 * pi / period)
         M + mean(y) +
           Av * cos(be + 2 * atan(om * tan((tt_r - al) / 2)))
       })
}

# ---- Plot ----
out_pdf <- "output/main_figures/SuppFig_clock_gene_profiles.pdf"
fig_local <- "output/supplementary/figures/supp_clock_gene_profiles.pdf"
dir.create("output/supplementary/figures", recursive = TRUE, showWarnings = FALSE)
dir.create("output/main_figures", recursive = TRUE, showWarnings = FALSE)

t_fine <- seq(0, 24, length.out = 401)

for (path in c(out_pdf, fig_local)) {
  cairo_pdf(path, width = 7.0, height = 5.6)
  pub_par(mfrow = c(2, 2), mar = c(4.0, 4.2, 2.4, 1.0))

  letters_seq <- letters
  pal_DT      <- pub_palette_detector()
  col_cos     <- unname(pal_DT["DCP"])    # blue
  col_fmm     <- unname(pal_DT["FMM"])    # orange

  for (i in seq_along(clock_genes)) {
    sym <- names(clock_genes)[i]
    id  <- clock_genes[[sym]]
    if (!(id %in% rownames(m))) {
      plot.new(); title(main = sprintf("%s: not measured", sym)); next
    }
    y <- as.numeric(m[id, ])

    cos_fit <- fit_cosinor(y, tods)
    fmm_fit <- fit_fmm(y, tods)

    yr <- range(y, cos_fit$predict(t_fine),
                if (!is.null(fmm_fit)) fmm_fit$predict(t_fine) else y)
    yspan <- diff(yr)
    yr    <- yr + c(-0.05, 0.05) * yspan

    plot(tods, y, pch = 19, col = "grey30", cex = 0.9,
         xlim = c(0, 24), ylim = yr,
         xlab = "Time (h)", ylab = "Expression",
         main = sym)
    panel_label(letters_seq[i])
    lines(t_fine, cos_fit$predict(t_fine),
          col = col_cos, lwd = 1.6, lty = 2)
    if (!is.null(fmm_fit)) {
      lines(t_fine, fmm_fit$predict(t_fine),
            col = col_fmm, lwd = 2.0, lty = 1)
      legend_lab <- c("data", "cosinor", "FMM")
      legend_col <- c("grey30", col_cos, col_fmm)
      legend_pch <- c(19, NA, NA)
      legend_lty <- c(NA, 2, 1)
      legend_lwd <- c(NA, 1.6, 2.0)
    } else {
      legend_lab <- c("data", "cosinor")
      legend_col <- c("grey30", col_cos)
      legend_pch <- c(19, NA)
      legend_lty <- c(NA, 2)
      legend_lwd <- c(NA, 1.6)
    }
    pub_legend("topright",
               legend = legend_lab,
               col    = legend_col,
               pch    = legend_pch,
               lty    = legend_lty,
               lwd    = legend_lwd,
               cex    = 0.7)
  }

  dev.off()
  cat(sprintf("Saved: %s\n", path))
}
