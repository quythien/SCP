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

.read_pilot <- function(species, dataset, tissue, condition) {
  if (.demo_mode) {
    m <- .read_manifest()
    hit <- m$species == species & m$dataset == dataset &
           m$tissue == tissue & m$condition == condition
    if (!any(hit)) stop("Pilot not in demo bundle.")
    readRDS(file.path(.demo_pilots, m$file[which(hit)[1]]))
  } else {
    SCP::scp_load_pilot(species, dataset, tissue, condition)
  }
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
      selectInput("rhy_thresh", "Threshold",
                  choices  = c("0.25" = 0.25, "0.20" = 0.20, "0.15" = 0.15,
                               "0.10" = 0.10, "0.05" = 0.05, "0.01" = 0.01,
                               "0.001" = 0.001, "0.0001" = 0.0001),
                  selected = 0.05),
      helpText(em("Sanity-check lens on the bundled pilot's 'Prop. rhythmic' line; not used in the simulation (the sim uses the Target FDR slider below).")),
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
      sliderInput("n_min", "Minimum N", min = 10, max = 100,
                  value = 20, step = 10),
      sliderInput("n_max", "Maximum N", min = 40, max = 400,
                  value = 120, step = 20),
      sliderInput("n_step", "Step", min = 10, max = 60,
                  value = 20, step = 10),
      helpText(em("Grid: N_min, N_min + step, ..., up to N_max.")),
      hr(),
      h4("Targets"),
      sliderInput("target_power", "Target power",
                  min = 0.50, max = 0.95, value = 0.80, step = 0.05),
      sliderInput("target_fdr",  "Target FDR",
                  min = 0.01, max = 0.20, value = 0.05, step = 0.01),
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
    p <- tryCatch(
      .read_pilot(input$species, input$dataset,
                  actual_tissue, input$condition),
      error = function(e) NULL
    )
    if (!is.null(p) && !inherits(p, "CircadianBioOptions"))
      class(p) <- c("CircadianBioOptions", class(p))
    p
  })

  # Track what was last simulated so we can warn when settings drift
  last_run_state <- reactiveVal(NULL)

  current_state <- reactive({
    list(species = input$species, dataset = input$dataset,
         tissue = input$tissue, condition = input$condition,
         K = input$K, design = input$design,
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
    # Threshold-aware prop_rhythmic from raw p-values if present
    prop_at_fdr <- p$prop_rhythmic %||% NA_real_
    if (!is.null(p$raw$pvalue)) {
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
      n_cores <- max(1L, min(4L, parallel::detectCores(logical = FALSE) - 1L))
      incProgress(0.2, detail = "Setting up design")

      # Defensive numeric coercion: some pilots may store fields as character / list
      .as_num <- function(x) suppressWarnings(as.numeric(x))
      for (fld in c("amplitude", "sigma_rhythmic", "phase",
                    "lBaselineExpr", "lOD")) {
        if (!is.null(p[[fld]])) p[[fld]] <- .as_num(p[[fld]])
      }
      if (!is.null(p$ngenes)) p$ngenes <- as.integer(.as_num(p$ngenes))
      if (!is.null(p$period)) p$period <- .as_num(p$period)
      # Sanity warning: too few rhythmic candidates -> simulation power is unreliable
      n_rhy_pilot <- length(p$sigma_rhythmic %||% numeric(0))
      if (n_rhy_pilot < 50L) {
        showNotification(
          sprintf(paste0("This pilot has only %d genes passing its internal ",
                         "rhythmic threshold. Power numbers may be unreliable. ",
                         "Try a stronger pilot (larger n or higher r-tilde) or ",
                         "interpret with caution."), n_rhy_pilot),
          type = "warning", duration = 12
        )
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
      # Cap genes at 2000 for snappier Shiny response
      if (!is.null(p$ngenes) && is.finite(p$ngenes) && p$ngenes > 2000L) {
        p$ngenes <- 2000L
      }
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
        nsims        = 10L,
        design       = design_type,
        cts          = cts_vec
      )
      analysis <- SCP::CircadianAnalysisOptions(
        alpha    = input$target_fdr,
        DCmethod = "DCP"
      )

      incProgress(0.3, detail = "Simulating")
      # Fixed seed so repeated clicks on the same pilot/grid give the same N80
      set.seed(2026)
      K_val <- as.integer(input$K)
      res <- tryCatch(
        SCP::runSimsSingleCohort(
          bio.opts      = p,
          design.opts   = design,
          analysis.opts = analysis,
          K             = K_val,
          mc.cores      = n_cores
        ),
        error = function(e) {
          msg <- sprintf("Simulation failed: %s", conditionMessage(e))
          showNotification(msg, type = "error", duration = 8)
          structure(list(error = msg), class = "sim_error")
        }
      )
      incProgress(1, detail = "Done")
      last_run_state(isolate(current_state()))
      res
    })
  })

  output$power_curve <- renderPlot({
    res <- sim_result()
    req(res)
    tryCatch(
      SCP::plotSingleCohortPower(res, fdr = input$target_fdr),
      error = function(e) {
        plot.new()
        title(main = sprintf("Plot error: %s", conditionMessage(e)))
      }
    )
  })


  output$recommended_n_text <- renderText({
    res <- sim_result()
    req(res)
    np <- tryCatch(
      SCP::npower(res, target_power = input$target_power,
                  fdr = input$target_fdr),
      error = function(e) list(n = NA_real_)
    )
    if (!is.na(np$n)) {
      return(sprintf("Recommended N = %.0f for %.0f%% power at FDR %.0f%% (from simulation grid)",
                     ceiling(np$n), 100 * input$target_power, 100 * input$target_fdr))
    }
    # Grid did not reach target. Report the max observed power AND a
    # closed-form per-gene-at-median-r-tilde estimate as separate numbers
    # so the user sees both honestly.
    max_pow <- if (length(np$power) > 0L) max(np$power, na.rm = TRUE) else NA_real_
    n_at_max <- if (is.finite(max_pow))
      res$sample_sizes[which.max(np$power)] else max(res$sample_sizes)
    p <- pilot()
    n_cf <- tryCatch(
      SCP::circaPowerApproxN80(
        bio.opts     = p,
        alpha        = input$target_fdr,
        target_power = input$target_power,
        n_search     = seq(20, 2000, by = 5)
      ),
      error = function(e) NA_real_
    )
    sim_msg <- if (is.finite(max_pow))
      sprintf("The largest sample size simulated (N = %d) reached only %.0f%% power at FDR %.0f%%. A larger N is needed to hit the %.0f%% target.",
              n_at_max, 100 * max_pow,
              100 * input$target_fdr, 100 * input$target_power)
    else
      sprintf("The simulation did not reach %.0f%% power at any N in your grid.",
              100 * input$target_power)
    cf_msg <- if (is.finite(n_cf))
      sprintf(" A rough closed-form estimate (based on the pilot's median r-tilde and ignoring multiple testing) suggests N around %.0f, but this number is usually optimistic for genome-wide power.",
              ceiling(n_cf))
    else
      " A closed-form estimate is not available for this pilot."
    paste0(sim_msg, cf_msg)
  })
}

shinyApp(ui, server)
