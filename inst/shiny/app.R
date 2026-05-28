# SCP Shiny app: simulation-based circadian power analysis (single-cohort MVP)
# Run locally: SCP::launchShiny()
# Deploy:      rsconnect::deployApp(system.file("shiny", package = "SCP"))

suppressPackageStartupMessages({
  library(shiny)
  library(SCP)
})

`%||%` <- function(a, b) if (!is.null(a) && length(a) > 0L) a else b

# ---------------------------------------------------------------------------
# UI
# ---------------------------------------------------------------------------

ui <- fluidPage(
  titlePanel("SCP: Circadian Study Sizer"),
  p(em("Pilot-calibrated power analysis for circadian transcriptomics studies."),
    "Pick a bundled pilot, set the target power and FDR, and read off the recommended sample size."),
  sidebarLayout(
    sidebarPanel(
      width = 4,
      h4("Pilot"),
      selectInput("species",   "Species",   choices = NULL),
      selectInput("dataset",   "Dataset",   choices = NULL),
      selectInput("tissue",    "Tissue",    choices = NULL),
      selectInput("condition", "Condition", choices = NULL),
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
    m <- tryCatch(SCP::scp_pilots(status = "ingested"),
                  error = function(e) NULL)
    if (is.null(m) || nrow(m) == 0L) return(NULL)
    m[m$sample_type %in% "in_vivo" | is.na(m$sample_type), , drop = FALSE]
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

  pilot <- reactive({
    req(input$species, input$dataset, input$tissue, input$condition)
    tryCatch(
      SCP::scp_load_pilot(input$species, input$dataset,
                          input$tissue, input$condition),
      error = function(e) NULL
    )
  })

  output$pilot_summary <- renderText({
    p <- pilot()
    if (is.null(p)) return("Pick a (species, dataset, tissue, condition) above.")
    amp <- p$amplitude   %||% NA_real_
    sig <- p$sigma_rhythmic %||% NA_real_
    r   <- if (length(amp) && length(sig) && length(amp) == length(sig))
             amp / sig else NA_real_
    r   <- r[is.finite(r) & r > 0]
    sprintf(
      paste0(
        "Pilot     : %s / %s / %s / %s\n",
        "Design    : %s\n",
        "Genes     : %s\n",
        "Prop. rhythmic (FDR 5%%) : %s\n",
        "Median r-tilde (A / sigma) : %s   IQR [%s, %s]"
      ),
      input$species, input$dataset, input$tissue, input$condition,
      manifest()[manifest()$species == input$species &
                 manifest()$dataset == input$dataset &
                 manifest()$tissue  == input$tissue   &
                 manifest()$condition == input$condition, "design"][1] %||% "unknown",
      format(p$ngenes %||% NA_integer_, big.mark = ","),
      if (is.null(p$prop_rhythmic)) "NA" else sprintf("%.1f%%", 100 * p$prop_rhythmic),
      if (length(r) == 0L) "NA" else sprintf("%.2f", stats::median(r)),
      if (length(r) == 0L) "NA" else sprintf("%.2f", stats::quantile(r, 0.25)),
      if (length(r) == 0L) "NA" else sprintf("%.2f", stats::quantile(r, 0.75))
    )
  })

  output$tod_plot <- renderPlot({
    p <- pilot()
    if (is.null(p)) {
      plot.new(); title(main = "No pilot selected")
      return(invisible())
    }
    times <- p$times %||% p$raw$times %||% numeric(0)
    if (length(times) == 0L) {
      plot.new(); title(main = "TOD metadata not available for this pilot")
      return(invisible())
    }
    t_mod <- times %% 24
    opar <- par(mar = c(3.2, 3.2, 1.6, 0.6), mgp = c(2, 0.6, 0))
    on.exit(par(opar))
    hist(t_mod, breaks = seq(0, 24, by = 1), col = "#a6cee3",
         border = "#1f78b4", xlim = c(0, 24),
         xlab = "Time of day (h, mod 24)", ylab = "Samples",
         main = sprintf("n = %d samples, %d distinct phase%s",
                        length(times), length(unique(round(t_mod, 1))),
                        if (length(unique(round(t_mod, 1))) == 1L) "" else "s"))
    rug(unique(round(t_mod, 2)), side = 3, col = "#e31a1c", lwd = 2)
  })

  sim_result <- eventReactive(input$run, {
    p <- pilot()
    req(p)
    withProgress(message = "Running simulation...", value = 0.1, {
      n_cores <- max(1L, min(4L, parallel::detectCores(logical = FALSE) - 1L))
      incProgress(0.2, detail = "Setting up design")

      design <- SCP::CircadianDesignOptions(
        sample_sizes = c(20, 40, 60, 80),
        nsims        = 20L,
        cts          = seq(0, 22, by = 4)
      )
      analysis <- SCP::CircadianAnalysisOptions(
        alpha    = input$target_fdr,
        K        = as.integer(input$K),
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
