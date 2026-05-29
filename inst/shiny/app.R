# SCP Shiny app: simulation-based circadian power analysis (single-cohort MVP)
# Run locally: SCP::launchShiny()
# Deploy:      rsconnect::deployApp(system.file("shiny", package = "SCP"))

suppressPackageStartupMessages({
  library(shiny)
  library(SCP)
})

`%||%` <- function(a, b) if (!is.null(a) && length(a) > 0L) a else b


# When SCP_SHINY_DEMO=1, the app reads pilots from ./pilots/ next to app.R
# instead of system.file("extdata","pilots", package = "SCP"). This is the
# slim 86-pilot subset bundled for the shinyapps.io free-tier demo.
.demo_mode    <- nzchar(Sys.getenv("SCP_SHINY_DEMO"))
.demo_pilots  <- file.path(getwd(), "pilots")

.read_manifest <- function() {
  if (.demo_mode && file.exists(file.path(.demo_pilots, "manifest.csv"))) {
    m <- utils::read.csv(file.path(.demo_pilots, "manifest.csv"),
                         stringsAsFactors = FALSE)
    m[m$status == "ingested", , drop = FALSE]
  } else {
    SCP::scp_pilots(status = "ingested")
  }
}

.read_pilot <- function(species, dataset, tissue, condition, K = 1) {
  K <- as.integer(K)
  if (.demo_mode) {
    m <- .read_manifest()
    hit <- m$species == species & m$dataset == dataset &
           m$tissue == tissue & m$condition == condition
    if (!any(hit)) stop("Pilot not in demo bundle.")
    f <- m$file[which(hit)[1]]
    if (K == 2L) {
      f2 <- sub("\\.rds$", "_2H.rds", f)
      p2 <- file.path(.demo_pilots, f2)
      if (!file.exists(p2)) stop("No two-harmonic (K=2) variant for this pilot.")
      f <- f2
    }
    readRDS(file.path(.demo_pilots, f))
  } else {
    SCP::scp_load_pilot(species, dataset, tissue, condition, K = K)
  }
}

# Does a K=2 (two-harmonic) variant exist for this pilot?
.has_k2_variant <- function(species, dataset, tissue, condition) {
  out <- tryCatch(.read_pilot(species, dataset, tissue, condition, K = 2),
                  error = function(e) NULL)
  !is.null(out)
}

# Safe worker count for parallel::mclapply. Caps at `cap`, leaves one core
# free, and is robust to environments where detectCores() returns NA (common
# in containers / shinyapps.io). mclapply uses fork(), which is unsupported on
# Windows, so force serial (1 core) there.
.safe_cores <- function(cap = 4L) {
  if (.Platform$OS.type == "windows") return(1L)
  nc <- tryCatch(parallel::detectCores(logical = FALSE), error = function(e) NA_integer_)
  if (is.na(nc) || nc < 1L)
    nc <- tryCatch(parallel::detectCores(), error = function(e) NA_integer_)
  if (is.na(nc) || nc < 1L) nc <- 1L
  max(1L, min(as.integer(cap), nc - 1L))
}

# Re-select the rhythmicity threshold from a pilot's stored per-gene table
# (rhythm_fit). Self-contained so it works in both demo mode (raw readRDS) and
# packaged mode. Mirrors SCP:::.reslice_pilot: rhythm_fit is p-sorted and capped,
# so thresholding is a prefix slice; the top-300 candidates drive the effect-size
# distribution. Pilots built before rhythm_fit existed are returned unchanged.
.app_apply_threshold <- function(p, stat = "p", thresh = 0.01, paired_sigma = TRUE) {
  if (is.null(p)) return(p)
  rf <- p$rhythm_fit
  if (is.null(rf) || nrow(rf) == 0L) return(p)   # legacy pilot: keep baked values
  cap <- p$pilot_cap %||% 0.2
  if (identical(stat, "q")) {
    n_total <- p$ngenes %||% nrow(rf)
    i <- seq_len(nrow(rf))
    sval <- rev(cummin(rev(rf$pvalue * n_total / i)))   # BH q-values
  } else {
    sval <- rf$pvalue
    if (thresh > cap) thresh <- cap                      # clamp raw-p to stored cap
  }
  cand  <- rf[sval < thresh, , drop = FALSE]
  K     <- min(as.integer(p$pilot_top_k %||% 300L), nrow(cand))
  estim <- if (K > 0L) cand[seq_len(K), , drop = FALSE] else cand[0, , drop = FALSE]
  p$prop_rhythmic  <- nrow(cand) / (p$ngenes %||% nrow(rf))
  if ("A2" %in% names(rf)) {
    # Two-harmonic pilot: keep all five fields gene-index-aligned (length K).
    p$amplitude      <- estim$A
    p$phase          <- estim$phi
    p$sigma_rhythmic <- estim$sigma
    p$amplitude2     <- estim$A2
    p$phase2         <- estim$phi2
    p$paired_2h      <- TRUE
  } else {
    p$amplitude      <- estim$A
    p$sigma_rhythmic <- estim$sigma
    p$phase          <- estim$phi[is.finite(estim$phi)]
  }
  # paired_sigma=TRUE keeps (A, sigma) gene-paired (realistic narrow r-tilde,
  # the marginal / Panel-A behavior). FALSE decouples sigma -> wide r-tilde
  # (stratified Panel-B/C). The app defaults to paired for recommended-N.
  if (!isTRUE(paired_sigma) && length(p$sigma_rhythmic) > 1L)
    p$sigma_rhythmic <- sample(p$sigma_rhythmic)
  p$paired_sigma <- isTRUE(paired_sigma)
  p
}

# ---------------------------------------------------------------------------
# UI
# ---------------------------------------------------------------------------

ui <- fluidPage(
  tags$head(tags$style(HTML("
    .center-titles h1, .center-titles h2, .center-titles h3, .center-titles h4, .center-titles p {
      text-align: center;
    }
  "))),
  div(class = "center-titles",
      titlePanel("Power Evaluation and Study Design for Circadian Biomarker Detection"),
      p(em("Pilot-calibrated power analysis for circadian transcriptomics studies."),
        br(),
        "Pick a bundled pilot, set the target power and FDR, and read off the recommended sample size.")
  ),
  sidebarLayout(
    sidebarPanel(
      width = 3,
      h4("Pilot rhythmicity threshold"),
      radioButtons("rhy_stat", NULL,
                   choices = c("Adjusted FDR (BH)" = "q",
                               "Raw p-value"       = "p"),
                   selected = "q", inline = TRUE),
      selectInput("rhy_thresh", "Threshold (alpha_pilot)",
                  choices  = c("0.20" = 0.20, "0.15" = 0.15,
                               "0.10" = 0.10, "0.05" = 0.05, "0.01" = 0.01,
                               "0.001" = 0.001, "0.0001" = 0.0001),
                  selected = 0.05),
      helpText(em("Defines which pilot genes count as 'rhythmic' for the simulation. The top-300 by p-value among genes passing this threshold drive the empirical effect-size distribution.")),
      hr(),
      h4("Pilot"),
      radioButtons(
        "pilot_source", NULL,
        choices  = c("Use a bundled pilot" = "bundled",
                     "Upload my own pilot" = "upload"),
        selected = "bundled", inline = TRUE
      ),
      conditionalPanel(
        condition = "input.pilot_source == 'bundled'",
        selectInput("species",    "Species",       choices = NULL),
        selectInput("design_filter", "Design",     choices = c("active", "passive")),
        selectInput("tissue",     "Tissue",        choices = NULL),
        selectInput("dataset",    "Study/dataset", choices = NULL),
        selectInput("condition",  "Condition",     choices = NULL)
      ),
      conditionalPanel(
        condition = "input.pilot_source == 'upload'",
        fileInput("upload_expr", "Expression CSV (genes x samples)",
                  accept = c(".csv", ".tsv", ".txt")),
        fileInput("upload_tod", "Times-of-day CSV (one value per sample, in hours)",
                  accept = c(".csv", ".tsv", ".txt")),
        helpText(em("The expression file's first column should be the gene name. ",
                    "The TOD file should be a single column of numeric hours, length = number of samples. ",
                    "Once both files load, the pilot is fit on the fly (5-15 sec).")),
        verbatimTextOutput("upload_status")
      ),
      hr(),
      h4("Detector"),
      radioButtons(
        "K", NULL,
        choices  = c("Single-harmonic (K = 1)" = "1",
                     "Two-harmonic (K = 2)"    = "2"),
        selected = "1"
      ),
      checkboxInput("eff_sens", "Effect-size sensitivity", value = FALSE),
      helpText(em("Also show power stratified by effect size r̃ = A/σ.")),
      hr(),
      h4("Sampling design"),
      radioButtons(
        "design", NULL,
        choices  = c("Active (controlled timecourse)" = "active",
                     "Passive (TOD-of-death)"         = "passive"),
        selected = "active"
      ),
      conditionalPanel(
        condition = "input.design == 'active'",
        sliderInput("active_step", "Sampling interval (h)",
                    min = 1, max = 12, value = 4, step = 1),
        helpText(em("Active grid: 0, step, 2*step, ... up to 24 h."))
      ),
      hr(),
      h4("Sample size grid"),
      sliderInput("n_min", "Minimum N", min = 4, max = 100,
                  value = 20, step = 2),
      sliderInput("n_max", "Maximum N", min = 20, max = 500,
                  value = 120, step = 10),
      sliderInput("n_step", "Step", min = 2, max = 60,
                  value = 20, step = 2),
      helpText(em("Grid: N_min, N_min + step, ..., up to N_max.")),
      hr(),
      h4("Targets"),
      sliderInput("target_power", "Target power",
                  min = 0.50, max = 0.95, value = 0.80, step = 0.05),
      selectInput("target_fdr",  "Target FDR",
                  choices  = c("0.01" = 0.01, "0.05" = 0.05, "0.10" = 0.10,
                               "0.15" = 0.15, "0.20" = 0.20),
                  selected = 0.05),
      hr(),
      actionButton("run", "Run simulation",
                   class = "btn-primary btn-block", width = "100%"),
      uiOutput("stale_banner")
    ),
    mainPanel(
      width = 9,
      fluidRow(
        column(7,
          h4("Pilot summary"),
          verbatimTextOutput("pilot_summary"),
          uiOutput("pilot_links")
        ),
        column(5,
          h4("Sampling design (TOD distribution)"),
          plotOutput("tod_plot", height = "180px")
        )
      ),
      conditionalPanel(
        condition = "input.run > 0",
        hr(),
        h4("Power curve"),
        plotOutput("power_curve", height = "560px"),
        hr(),
        div(style = "padding: 12px; background:#f4f8fc; border-left: 4px solid #2c7fb8;",
            h3(textOutput("recommended_n_text"), style = "margin: 0;"))
      ),
      conditionalPanel(
        condition = "input.run == 0",
        br(),
        helpText(em("Click 'Run simulation' to compute the power curve."))
      )
    )
  )
)

# ---------------------------------------------------------------------------
# Server
# ---------------------------------------------------------------------------

server <- function(input, output, session) {

  manifest <- reactive({
    m <- tryCatch(.read_manifest(), error = function(e) NULL)
    if (is.null(m) || nrow(m) == 0L) return(NULL)
    if ("sample_type" %in% names(m))
      m <- m[m$sample_type %in% "in_vivo" | is.na(m$sample_type), , drop = FALSE]
    m
  })

  observe({
    m <- manifest()
    req(m)
    sp <- sort(unique(m$species))
    updateSelectInput(session, "species", choices = sp,
                      selected = if ("human" %in% sp) "human" else sp[1])
  })

  # Use tissue_canonical for the user-facing tissue dropdown when available
  .tissue_display <- function(m) {
    if ("tissue_canonical" %in% names(m)) m$tissue_canonical else m$tissue
  }

  # Cascade: species -> design -> tissue (canonical) -> dataset -> condition
  observeEvent(input$species, {
    m <- manifest()
    req(m, input$species)
    designs <- sort(unique(m$design[m$species == input$species]))
    updateSelectInput(session, "design_filter", choices = designs,
                      selected = if ("active" %in% designs) "active" else designs[1])
  })

  observeEvent(list(input$species, input$design_filter), {
    m <- manifest()
    req(m, input$species, input$design_filter)
    sub <- m[m$species == input$species & m$design == input$design_filter, , drop = FALSE]
    ts <- sort(unique(.tissue_display(sub)))
    updateSelectInput(session, "tissue", choices = ts, selected = ts[1])
  })

  observeEvent(list(input$species, input$design_filter, input$tissue), {
    m <- manifest()
    req(m, input$species, input$design_filter, input$tissue)
    sub <- m[m$species == input$species & m$design == input$design_filter, , drop = FALSE]
    sub <- sub[.tissue_display(sub) == input$tissue, , drop = FALSE]
    ds <- sort(unique(sub$dataset))
    updateSelectInput(session, "dataset", choices = ds, selected = ds[1])
  })

  observeEvent(list(input$species, input$design_filter, input$tissue, input$dataset), {
    m <- manifest()
    req(m, input$species, input$design_filter, input$tissue, input$dataset)
    sub <- m[m$species == input$species & m$design == input$design_filter &
             m$dataset == input$dataset, , drop = FALSE]
    sub <- sub[.tissue_display(sub) == input$tissue, , drop = FALSE]
    cd <- sort(unique(sub$condition))
    updateSelectInput(session, "condition", choices = cd, selected = cd[1])
  })

  # ---- custom pilot upload ------------------------------------------------
  uploaded_pilot_state <- reactiveValues(pilot = NULL, msg = "Awaiting upload.")

  observe({
    req(input$pilot_source == "upload")
    req(input$upload_expr, input$upload_tod)
    expr_file <- input$upload_expr$datapath
    tod_file  <- input$upload_tod$datapath

    res <- tryCatch({
      sep <- if (grepl("\\.tsv$|\\.txt$", input$upload_expr$name, TRUE)) "\t" else ","
      ex  <- utils::read.csv(expr_file, row.names = 1, check.names = FALSE,
                             stringsAsFactors = FALSE, sep = sep)
      tod <- as.numeric(utils::read.csv(tod_file, header = FALSE,
                                         stringsAsFactors = FALSE)[, 1])
      if (ncol(ex) != length(tod))
        stop(sprintf("Sample-count mismatch: expression has %d columns, TOD has %d values.",
                     ncol(ex), length(tod)))
      mat <- as.matrix(ex)
      mode(mat) <- "numeric"
      withProgress(message = "Fitting pilot...", value = 0.3, {
        bio <- SCP::estimate_circadian_params(mat, tod, verbose = FALSE)
        incProgress(1)
        bio$times <- tod
        bio
      })
    }, error = function(e) e)

    if (inherits(res, "error")) {
      uploaded_pilot_state$pilot <- NULL
      uploaded_pilot_state$msg <- sprintf("Upload error: %s", conditionMessage(res))
    } else {
      uploaded_pilot_state$pilot <- res
      uploaded_pilot_state$msg <- sprintf(
        "OK: %d genes x %d samples fit. Median r-tilde = %.2f.",
        ncol(res$raw$pvalue %||% rep(NA, res$ngenes)),
        length(input$upload_tod$datapath),
        stats::median((res$amplitude %||% NA) / (res$sigma_rhythmic %||% NA),
                      na.rm = TRUE)
      )
    }
  })

  uploaded_pilot <- reactive({ uploaded_pilot_state$pilot })

  output$upload_status <- renderText({ uploaded_pilot_state$msg })

  # Resolve canonical tissue back to the actual manifest tissue value
  .resolve_actual_tissue <- function() {
    m <- manifest()
    sub <- m[m$species  == input$species  &
             m$design   == input$design_filter &
             m$dataset  == input$dataset  &
             m$condition == input$condition, , drop = FALSE]
    sub <- sub[.tissue_display(sub) == input$tissue, , drop = FALSE]
    if (nrow(sub) > 0L) sub$tissue[1] else input$tissue
  }

  # ---- bundled pilot ------------------------------------------------------
  pilot <- reactive({
    if (isTRUE(input$pilot_source == "upload")) return(uploaded_pilot())
    req(input$species, input$dataset, input$tissue, input$condition)
    actual_tissue <- .resolve_actual_tissue()
    K_req <- as.integer(input$K %||% 1L)
    p <- tryCatch(
      .read_pilot(input$species, input$dataset,
                  actual_tissue, input$condition, K = K_req),
      error = function(e) {
        if (K_req == 2L)
          showNotification(
            paste0("No two-harmonic (K=2) pilot has been built for this dataset. ",
                   "Switch to K = 1, or build a K=2 variant."),
            type = "warning", duration = 8)
        else
          showNotification(
            sprintf("Could not load pilot [%s / %s / %s / %s]: %s",
                    input$species, input$dataset, actual_tissue, input$condition,
                    conditionMessage(e)),
            type = "error", duration = 12)
        NULL
      }
    )
    if (!is.null(p) && !inherits(p, "CircadianBioOptions"))
      class(p) <- c("CircadianBioOptions", class(p))
    # Apply the user-selected rhythmicity threshold so prop_rhythmic and the
    # effect-size distribution reflect alpha_pilot. Never let a threshold
    # hiccup blank out a successfully-loaded pilot: fall back to the loaded
    # object and surface the error.
    p <- tryCatch(
      .app_apply_threshold(p, input$rhy_stat, as.numeric(input$rhy_thresh)),
      error = function(e) {
        showNotification(sprintf("Threshold step failed (%s); showing pilot defaults.",
                                 conditionMessage(e)), type = "warning", duration = 8)
        p
      }
    )
    p
  })

  # Track what was last simulated so we can warn when settings drift
  last_run_state <- reactiveVal(NULL)

  current_state <- reactive({
    list(species = input$species, dataset = input$dataset,
         tissue = input$tissue, condition = input$condition,
         rhy_stat = input$rhy_stat, rhy_thresh = input$rhy_thresh,
         K = input$K, eff_sens = input$eff_sens, design = input$design,
         active_step = input$active_step,
         n_min = input$n_min, n_max = input$n_max, n_step = input$n_step,
         target_fdr = input$target_fdr)
  })

  output$stale_banner <- renderUI({
    if (is.null(last_run_state())) return(NULL)
    if (identical(current_state(), last_run_state())) return(NULL)
    div(style = "margin-top: 8px; padding: 8px; background: #fff3cd; border-left: 4px solid #d68f00; font-size: 0.85em;",
        strong("Settings changed."), br(),
        "Click Run again to refresh the curve and recommended N.")
  })

  output$pilot_links <- renderUI({
    req(input$species, input$dataset, input$tissue, input$condition)
    m <- manifest()
    sub <- m[m$species == input$species & m$design == input$design_filter &
             m$dataset == input$dataset & m$condition == input$condition, , drop = FALSE]
    sub <- sub[.tissue_display(sub) == input$tissue, , drop = FALSE]
    if (nrow(sub) == 0L) return(NULL)
    acc <- sub$accession[1] %||% ""
    cit <- sub$citation[1] %||% ""
    url <- ""
    if (grepl("^GSE[0-9]+", acc)) {
      gse <- regmatches(acc, regexpr("GSE[0-9]+", acc))
      url <- sprintf("https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=%s", gse)
    } else if (grepl("GTEx", acc, ignore.case = TRUE)) {
      url <- "https://gtexportal.org/"
    } else if (grepl("CMC|Seney", acc, ignore.case = TRUE)) {
      url <- "https://www.synapse.org/#!Synapse:syn2759792"
    } else if (grepl("E-MTAB", acc)) {
      url <- sprintf("https://www.ebi.ac.uk/biostudies/arrayexpress/studies/%s", acc)
    }
    HTML(sprintf(
      "<div style='margin-top:6px; font-size:0.9em;'><b>Accession:</b> %s%s<br/><b>Citation:</b> <i>%s</i></div>",
      acc,
      if (nzchar(url)) sprintf(" &nbsp;<a href='%s' target='_blank'>open page &rarr;</a>", url) else "",
      cit
    ))
  })

  output$pilot_summary <- renderText({
    p <- pilot()
    if (is.null(p)) return("Pick a (species, dataset, tissue, condition) above.")
    # Prefer cached r-tilde summaries baked into the pilot (correct gene pairing)
    r_med <- p$r_median %||% NA_real_
    r_q25 <- p$r_q25    %||% NA_real_
    r_q75 <- p$r_q75    %||% NA_real_
    if (!is.finite(r_med)) {
      amp <- p$amplitude %||% numeric(0); sig <- p$sigma_rhythmic %||% numeric(0)
      n_min <- min(length(amp), length(sig))
      if (n_min > 0L) {
        r <- amp[seq_len(n_min)] / sig[seq_len(n_min)]
        r <- r[is.finite(r) & r > 0]
        if (length(r) > 0) {
          r_med <- median(r); r_q25 <- quantile(r, 0.25); r_q75 <- quantile(r, 0.75)
        }
      }
    }
    # prop_rhythmic is already re-derived at the selected threshold in pilot()
    # via the stored rhythm_fit table; fall back to the per-gene p-vector or the
    # baked value for pilots built before rhythm_fit existed.
    prop_at_fdr <- p$prop_rhythmic %||% NA_real_
    if (is.null(p$rhythm_fit) && !is.null(p$raw$pvalue)) {
      pv <- as.numeric(p$raw$pvalue)
      pv <- pv[is.finite(pv)]
      stat <- if (input$rhy_stat == "q") p.adjust(pv, "BH") else pv
      prop_at_fdr <- mean(stat < as.numeric(input$rhy_thresh), na.rm = TRUE)
    }
    rhy_label <- if (input$rhy_stat == "q") sprintf("FDR %.0f%%", 100 * as.numeric(input$rhy_thresh)) else sprintf("p < %.3g", as.numeric(input$rhy_thresh))
    sprintf(
      paste0(
        "Pilot     : %s / %s / %s / %s\n",
        "Design    : %s\n",
        "Samples (n) : %s\n",
        "Genes     : %s\n",
        "Prop. rhythmic (%s) : %s\n",
        "Median r-tilde (A / sigma) : %s   IQR [%s, %s]"
      ),
      input$species, input$dataset, input$tissue, input$condition,
      {
        mrow <- manifest()[manifest()$species == input$species &
                           manifest()$dataset == input$dataset &
                           .tissue_display(manifest()) == input$tissue   &
                           manifest()$condition == input$condition, , drop = FALSE]
        if (nrow(mrow) > 0L) mrow$design[1] else "unknown"
      },
      {
        mrow <- manifest()[manifest()$species == input$species &
                           manifest()$dataset == input$dataset &
                           .tissue_display(manifest()) == input$tissue   &
                           manifest()$condition == input$condition, , drop = FALSE]
        nn <- if (nrow(mrow) > 0L && "n" %in% names(mrow)) mrow$n[1] else NA_integer_
        if (is.null(nn) || is.na(nn)) "NA" else format(nn, big.mark = ",")
      },
      {
        ng <- p$ngenes
        if (is.null(ng) || is.na(ng)) {
          mrow <- manifest()[manifest()$species == input$species &
                             manifest()$dataset == input$dataset &
                             .tissue_display(manifest()) == input$tissue   &
                             manifest()$condition == input$condition, , drop = FALSE]
          ng <- if (nrow(mrow) > 0L && "ngenes" %in% names(mrow)) mrow$ngenes[1] else NA_integer_
        }
        if (is.null(ng) || is.na(ng)) "NA" else format(ng, big.mark = ",")
      },
      rhy_label,
      if (!is.finite(prop_at_fdr)) "NA" else sprintf("%.1f%%", 100 * prop_at_fdr),
      if (!is.finite(r_med)) "NA" else sprintf("%.2f", r_med),
      if (!is.finite(r_q25)) "NA" else sprintf("%.2f", r_q25),
      if (!is.finite(r_q75)) "NA" else sprintf("%.2f", r_q75)
    )
  })

  output$tod_plot <- renderPlot({
    p <- pilot()
    # balanced left/right margins so the title visually centers on the plot
    par(mar = c(3.4, 3.4, 2.2, 3.0), mgp = c(2.2, 0.6, 0))
    if (is.null(p)) {
      plot(0, 0, type = "n", axes = FALSE, xlab = "", ylab = "")
      text(0, 0, "No pilot loaded", cex = 1.2, col = "grey40")
      return()
    }
    # Prefer raw times when preserved; otherwise reconstruct from manifest
    times  <- p$times %||% p$raw$times %||% numeric(0)
    phases <- numeric(0)
    n_samp <- NA_integer_
    design_label <- NA_character_

    if (length(times) > 0L) {
      phases <- times %% 24
      n_samp <- length(times)
      design_label <- "raw TOD"
    } else if (input$pilot_source == "bundled") {
      mrow <- manifest()[manifest()$species == input$species &
                         manifest()$dataset == input$dataset &
                         .tissue_display(manifest()) == input$tissue   &
                         manifest()$condition == input$condition, , drop = FALSE]
      if (nrow(mrow) > 0L && "tod_phases" %in% names(mrow) &&
          !is.na(mrow$tod_phases[1]) && nzchar(mrow$tod_phases[1]) &&
          !grepl("phases\\)$", mrow$tod_phases[1])) {
        ph <- suppressWarnings(as.numeric(strsplit(mrow$tod_phases[1], ",")[[1]]))
        ph <- ph[is.finite(ph)]
        if (length(ph) > 0) {
          phases <- ph
          n_samp <- mrow$n[1] %||% length(ph)
          design_label <- "design summary"
        }
      }
    }

    if (length(phases) == 0L) {
      plot(0, 0, type = "n", axes = FALSE, xlab = "", ylab = "")
      text(0, 0, "TOD metadata not available\n(no raw times, no manifest summary)",
           cex = 1.05, col = "grey40")
      return()
    }

    n_u <- length(unique(round(phases, 1)))
    h <- hist(phases, breaks = seq(0, 24, by = 1), plot = FALSE)
    plot(h, col = "#a6cee3", border = "#1f78b4", xlim = c(0, 24),
         xlab = "Time of day (h, mod 24)", ylab = "Samples", main = "",
         xaxt = "n")
    axis(side = 1, at = seq(0, 24, by = 4))
    title(main = sprintf("n = %s samples, %d distinct time point%s",
                         ifelse(is.na(n_samp), "?", as.character(n_samp)),
                         n_u, if (n_u == 1L) "" else "s"),
          cex.main = 1.0, line = 0.6, adj = 0.5, outer = FALSE)
  })

  sim_result <- eventReactive(input$run, {
    p <- pilot()
    req(p)
    withProgress(message = "Running simulation...", value = 0.1, {
      n_cores <- .safe_cores(4L)
      incProgress(0.2, detail = "Setting up design")

      # Defensive numeric coercion: some pilots may store fields as character / list
      .as_num <- function(x) suppressWarnings(as.numeric(x))
      for (fld in c("amplitude", "sigma_rhythmic", "phase",
                    "lBaselineExpr", "lOD")) {
        if (!is.null(p[[fld]])) p[[fld]] <- .as_num(p[[fld]])
      }
      # Fixed simulation gene count: matches the manuscript figure scripts
      # (NGENES_SIM = 2000) and keeps the interactive app responsive. The
      # runner resamples the per-gene baseline vectors to this length.
      p$ngenes <- 2000L
      if (!is.null(p$period)) p$period <- .as_num(p$period)
      # Apply the user's pilot rhythmicity threshold to redefine the rhythmic
      # gene set (top-K = 300 by p-value among genes passing the threshold).
      # This means the sidebar threshold actually drives which genes' empirical
      # (amplitude, sigma, phase) populate the simulation distribution.
      if (!is.null(p$raw) && !is.null(p$raw$A) && !is.null(p$raw$pvalue)) {
        thresh <- as.numeric(input$rhy_thresh)
        stat <- if (isTRUE(input$rhy_stat == "q"))
          p.adjust(p$raw$pvalue, "BH") else p$raw$pvalue
        good <- !is.na(p$raw$A) & p$raw$A > 0 &
                is.finite(p$raw$sigma) & p$raw$sigma > 0
        cand <- which(good & !is.na(stat) & stat < thresh)
        if (length(cand) == 0L) cand <- which(good)   # fallback
        K <- min(300L, length(cand))
        topK <- cand[order(p$raw$pvalue[cand])][seq_len(K)]
        p$amplitude      <- p$raw$A[topK]
        p$sigma_rhythmic <- p$raw$sigma[topK]
        p$phase          <- p$raw$phi[topK]
        p$lOD            <- log(p$raw$sigma[good])
        p$lBaselineExpr  <- p$raw$M
      }
      # Harmonize amplitude / phase to match sigma_rhythmic length when the
      # pilot was built with mismatched gene-sets (caught in some agent-
      # ingested GTEx pilots: amplitude length 1122 vs sigma 300).
      if (!is.null(p$amplitude) && !is.null(p$sigma_rhythmic)) {
        n_sigma <- length(p$sigma_rhythmic)
        if (length(p$amplitude) != n_sigma && n_sigma > 0L) {
          # Truncate or recycle to match sigma
          if (length(p$amplitude) > n_sigma)
            p$amplitude <- p$amplitude[seq_len(n_sigma)]
          else
            p$amplitude <- rep_len(p$amplitude, n_sigma)
        }
        if (!is.null(p$phase) && length(p$phase) != n_sigma) {
          if (length(p$phase) > n_sigma) p$phase <- p$phase[seq_len(n_sigma)]
          else p$phase <- rep_len(p$phase, n_sigma)
        }
      }
      # Use the pilot's full gene count (matches the manuscript figure scripts)
      # Build sample-size grid from user-chosen min/max/step
      n_grid <- seq(as.integer(input$n_min),
                    as.integer(input$n_max),
                    by = as.integer(input$n_step))
      if (length(n_grid) < 2L) {
        n_grid <- as.integer(input$n_min) +
                  as.integer(input$n_step) * 0:3
      }
      # cts: q4h grid for active; for passive use the pilot's empirical times
      # if preserved, else fall back to the q4h grid as a placeholder
      design_type <- input$design
      cts_vec <- if (design_type == "passive") {
        tt <- p$times %||% p$raw$times %||% NULL
        if (is.null(tt) || length(tt) == 0L) {
          # Pilot has no raw times; try the manifest's tod_phases as fallback
          m <- manifest()
          mrow <- m[m$species == input$species & m$design == input$design_filter &
                    m$dataset == input$dataset & m$condition == input$condition, , drop = FALSE]
          mrow <- mrow[.tissue_display(mrow) == input$tissue, , drop = FALSE]
          ph <- if (nrow(mrow) > 0L && "tod_phases" %in% names(mrow) &&
                    !is.na(mrow$tod_phases[1]) && nzchar(mrow$tod_phases[1]) &&
                    !grepl("phases\\)$", mrow$tod_phases[1])) {
            suppressWarnings(as.numeric(strsplit(mrow$tod_phases[1], ",")[[1]]))
          } else NULL
          if (!is.null(ph) && length(ph) > 0L) ph else seq(0, 22, by = 4)
        } else as.numeric(tt)
      } else {
        step <- suppressWarnings(as.numeric(input$active_step))
        if (!is.finite(step) || step <= 0) step <- 4
        seq(0, 24 - step, by = step)
      }
      cts_vec <- as.numeric(cts_vec)
      cts_vec <- cts_vec[is.finite(cts_vec)]
      if (length(cts_vec) == 0L) cts_vec <- seq(0, 22, by = 4)
      design <- SCP::CircadianDesignOptions(
        sample_sizes = n_grid,
        nsims        = 20L,   # fixed for a responsive, reproducible interactive run
        design       = design_type,
        cts          = cts_vec
      )
      analysis <- SCP::CircadianAnalysisOptions(
        alpha    = as.numeric(input$target_fdr),
        DCmethod = "DCP"
      )

      incProgress(0.3, detail = "Simulating")
      # Use the same seed as the manuscript figure scripts (GLOBAL_SEED = 2025)
      K_val <- as.integer(input$K)
      run_sim <- function(bio) tryCatch(
        { set.seed(2025)
          SCP::runSimsSingleCohort(bio.opts = bio, design.opts = design,
                                   analysis.opts = analysis, K = K_val,
                                   mc.cores = n_cores) },
        error = function(e) {
          showNotification(sprintf("Simulation failed: %s", conditionMessage(e)),
                           type = "error", duration = 8)
          structure(list(error = conditionMessage(e)), class = "sim_error")
        }
      )
      # Panel A: paired (A, sigma) draw -> realistic marginal power (p is paired).
      res <- run_sim(p)
      # Effect-size sensitivity (Panels B/C): unpaired draw -> wide r-tilde, so
      # the stratified panels span the full effect-size range. Re-derived from
      # the same rhythm_fit table at the same threshold (no extra fitting).
      if (isTRUE(input$eff_sens) && !inherits(res, "sim_error")) {
        incProgress(0.3, detail = "Effect-size sensitivity")
        p_unpaired <- .app_apply_threshold(p, input$rhy_stat,
                                           as.numeric(input$rhy_thresh),
                                           paired_sigma = FALSE)
        res$.bc <- run_sim(p_unpaired)
      }
      incProgress(1, detail = "Done")
      last_run_state(isolate(current_state()))
      res
    })
  })

  output$power_curve <- renderPlot({
    res <- sim_result()
    req(res)
    tryCatch({
      tfd <- as.numeric(input$target_fdr)
      tpw <- as.numeric(input$target_power)
      # If effect-size sensitivity is on, draw the stratified panels from the
      # unpaired run (.bc) while keeping Panel A from the paired run (panel_a_res).
      bc <- res$.bc
      if (!is.null(bc) && !inherits(bc, "sim_error")) {
        SCP::plotSingleCohortPower(bc, panel_a_res = res,
                                    fdr = tfd, fdr_thresholds = tfd,
                                    panel_fdr = tfd, vline_fdr = tfd,
                                    vline_power = tpw, r_max = NULL)
      } else {
        SCP::plotSingleCohortPower(res,
                                    fdr            = tfd,
                                    fdr_thresholds = tfd,
                                    panel_fdr      = tfd,
                                    vline_fdr      = tfd,
                                    vline_power    = tpw,
                                    r_max          = NULL)
      }
      },
      error = function(e) {
        plot.new()
        title(main = sprintf("Plot error: %s", conditionMessage(e)))
      }
    )
  })


  output$recommended_n_text <- renderText({
    res <- sim_result()
    req(res)
    tpw <- as.numeric(input$target_power)
    tfd <- as.numeric(input$target_fdr)
    # Use the runner's stored marginal_power matrix [N x nsims] — this is
    # exactly what the runner prints to the console at each N and what
    # plotSingleCohortPower's Panel A renders, so the text cannot drift
    # from the plot.
    mp <- res$marginal_power
    if (is.null(mp)) {
      return("Simulation did not return a marginal_power matrix.")
    }
    mean_pow <- rowMeans(mp, na.rm = TRUE)
    sizes    <- res$sample_sizes
    above    <- which(mean_pow >= tpw)
    if (length(above) > 0L) {
      j2 <- min(above)
      if (j2 == 1L) {
        rec_n <- sizes[1L]
      } else {
        j1 <- j2 - 1L
        rec_n <- ceiling(sizes[j1] +
                         (tpw - mean_pow[j1]) /
                         (mean_pow[j2] - mean_pow[j1]) *
                         (sizes[j2] - sizes[j1]))
      }
      return(sprintf("Recommended N = %d for %.0f%% power at FDR %.0f%% (from simulation grid)",
                     rec_n, 100 * tpw, 100 * tfd))
    }
    n_largest <- max(sizes)
    last_pow  <- mean_pow[which(sizes == n_largest)[1]]
    sprintf("The largest sample size simulated (N = %d) reached %.0f%% power at FDR %.0f%%. A larger N is needed to hit the %.0f%% target.",
            n_largest, 100 * last_pow, 100 * tfd, 100 * tpw)
  })
}

shinyApp(ui, server)
