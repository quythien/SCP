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
      width = 4,
      h4("Pilot"),
      radioButtons(
        "pilot_source", NULL,
        choices  = c("Use a bundled pilot" = "bundled",
                     "Upload my own pilot" = "upload"),
        selected = "bundled", inline = TRUE
      ),
      conditionalPanel(
        condition = "input.pilot_source == 'bundled'",
        selectInput("species",   "Species",   choices = NULL),
        selectInput("dataset",   "Dataset",   choices = NULL),
        selectInput("tissue",    "Tissue",    choices = NULL),
        selectInput("condition", "Condition", choices = NULL)
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
      h4("Targets"),
      sliderInput("target_power", "Target power",
                  min = 0.50, max = 0.95, value = 0.80, step = 0.05),
      sliderInput("target_fdr",  "Target FDR",
                  min = 0.01, max = 0.20, value = 0.05, step = 0.01),
      hr(),
      actionButton("run", "Run simulation",
                   class = "btn-primary btn-block", width = "100%"),
      helpText("A single run takes ~2-5 seconds.")
    ),
    mainPanel(
      width = 8,
      fluidRow(
        column(7,
          h4("Pilot summary"),
          verbatimTextOutput("pilot_summary")
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
        plotOutput("power_curve", height = "420px"),
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

  observeEvent(input$species, {
    m <- manifest()
    req(m, input$species)
    ds <- sort(unique(m$dataset[m$species == input$species]))
    updateSelectInput(session, "dataset", choices = ds, selected = ds[1])
  })

  observeEvent(list(input$species, input$dataset), {
    m <- manifest()
    req(m, input$species, input$dataset)
    ts <- sort(unique(m$tissue[m$species == input$species &
                               m$dataset == input$dataset]))
    updateSelectInput(session, "tissue", choices = ts, selected = ts[1])
  })

  observeEvent(list(input$species, input$dataset, input$tissue), {
    m <- manifest()
    req(m, input$species, input$dataset, input$tissue)
    cd <- sort(unique(m$condition[m$species == input$species &
                                  m$dataset == input$dataset &
                                  m$tissue  == input$tissue]))
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

  # ---- bundled pilot ------------------------------------------------------
  pilot <- reactive({
    if (isTRUE(input$pilot_source == "upload")) return(uploaded_pilot())
    req(input$species, input$dataset, input$tissue, input$condition)
    p <- tryCatch(
      .read_pilot(input$species, input$dataset,
                  input$tissue, input$condition),
      error = function(e) NULL
    )
    # ensure the class is present so runSimsSingleCohort accepts it
    if (!is.null(p) && !inherits(p, "CircadianBioOptions"))
      class(p) <- c("CircadianBioOptions", class(p))
    p
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
    # FDR-aware prop_rhythmic (recompute from raw p-values if present)
    prop_at_fdr <- p$prop_rhythmic %||% NA_real_
    if (!is.null(p$raw$pvalue)) {
      q <- p.adjust(p$raw$pvalue, "BH")
      prop_at_fdr <- mean(q < input$target_fdr, na.rm = TRUE)
    }
    sprintf(
      paste0(
        "Pilot     : %s / %s / %s / %s\n",
        "Design    : %s\n",
        "Genes     : %s\n",
        "Prop. rhythmic (FDR %.0f%%) : %s\n",
        "Median r-tilde (A / sigma) : %s   IQR [%s, %s]"
      ),
      input$species, input$dataset, input$tissue, input$condition,
      {
        mrow <- manifest()[manifest()$species == input$species &
                           manifest()$dataset == input$dataset &
                           manifest()$tissue  == input$tissue   &
                           manifest()$condition == input$condition, , drop = FALSE]
        if (nrow(mrow) > 0L) mrow$design[1] else "unknown"
      },
      {
        ng <- p$ngenes
        if (is.null(ng) || is.na(ng)) {
          mrow <- manifest()[manifest()$species == input$species &
                             manifest()$dataset == input$dataset &
                             manifest()$tissue  == input$tissue   &
                             manifest()$condition == input$condition, , drop = FALSE]
          ng <- if (nrow(mrow) > 0L && "ngenes" %in% names(mrow)) mrow$ngenes[1] else NA_integer_
        }
        if (is.null(ng) || is.na(ng)) "NA" else format(ng, big.mark = ",")
      },
      100 * input$target_fdr,
      if (!is.finite(prop_at_fdr)) "NA" else sprintf("%.1f%%", 100 * prop_at_fdr),
      if (!is.finite(r_med)) "NA" else sprintf("%.2f", r_med),
      if (!is.finite(r_q25)) "NA" else sprintf("%.2f", r_q25),
      if (!is.finite(r_q75)) "NA" else sprintf("%.2f", r_q75)
    )
  })

  output$tod_plot <- renderPlot({
    p <- pilot()
    par(mar = c(3.4, 3.4, 1.8, 0.6), mgp = c(2.2, 0.6, 0))
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
                         manifest()$tissue  == input$tissue   &
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
         xlab = "Time of day (h, mod 24)", ylab = "Samples", main = "")
    title(main = sprintf("n = %s samples, %d distinct phase%s (%s)",
                         ifelse(is.na(n_samp), "?", as.character(n_samp)),
                         n_u, if (n_u == 1L) "" else "s", design_label),
          cex.main = 1.0, line = 0.6)
    rug(unique(round(phases, 2)), side = 3, col = "#e31a1c", lwd = 2)
  })

  sim_result <- eventReactive(input$run, {
    p <- pilot()
    req(p)
    withProgress(message = "Running simulation...", value = 0.1, {
      n_cores <- max(1L, min(4L, parallel::detectCores(logical = FALSE) - 1L))
      incProgress(0.2, detail = "Setting up design")

      design <- SCP::CircadianDesignOptions(
        sample_sizes = c(20, 40, 60, 80),
        nsims        = 10L,
        cts          = seq(0, 22, by = 4),
        ngenes       = 2000L
      )
      analysis <- SCP::CircadianAnalysisOptions(
        alpha    = input$target_fdr,
        DCmethod = if (as.integer(input$K) == 1L) "DCP" else "FMM"
      )

      incProgress(0.3, detail = "Simulating")
      res <- SCP::runSimsSingleCohort(
        bio.opts      = p,
        design.opts   = design,
        analysis.opts = analysis,
        mc.cores      = n_cores
      )
      incProgress(1, detail = "Done")
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
    if (is.na(np$n)) {
      sprintf("Target %.0f%% power at FDR %.0f%% not reached in the simulated grid (N = 20-80). Try a larger N or stronger pilot.",
              100 * input$target_power, 100 * input$target_fdr)
    } else {
      sprintf("Recommended N = %.0f for %.0f%% power at FDR %.0f%%",
              ceiling(np$n), 100 * input$target_power, 100 * input$target_fdr)
    }
  })
}

shinyApp(ui, server)
