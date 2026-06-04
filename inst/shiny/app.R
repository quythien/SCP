# SCP Shiny app: simulation-based circadian power analysis (single-cohort MVP)
# Run locally: SCP::launchShiny()
# Deploy:      rsconnect::deployApp(system.file("shiny", package = "SCP"))

suppressPackageStartupMessages({
  library(shiny)
  library(SCP)
})

# Lift the 5 MB default upload cap to 5 GB so a full gene matrix uploads as-is.
options(shiny.maxRequestSize = 5 * 1024^3)

`%||%` <- function(a, b) if (!is.null(a) && length(a) > 0L) a else b

# Capitalize the first letter for DISPLAY only; underlying values are unchanged.
.cap1 <- function(x) {
  x <- as.character(x)
  ifelse(nchar(x) > 0L, paste0(toupper(substr(x, 1, 1)), substr(x, 2, nchar(x))), x)
}
# Dropdowns read "Human"/"Active"; loader keys stay "human"/"active".
.cap_choices <- function(v) stats::setNames(as.character(v), .cap1(v))


# SCP_SHINY_DEMO=1 reads the slim 86-pilot subset from ./pilots/ (shinyapps.io
# demo) instead of the installed extdata/pilots.
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

# Worker count for mclapply: cap, leave a core free, handle detectCores()==NA,
# and fall back to serial on Windows (no fork()).
.safe_cores <- function(cap = 4L) {
  if (.Platform$OS.type == "windows") return(1L)
  nc <- tryCatch(parallel::detectCores(logical = FALSE), error = function(e) NA_integer_)
  if (is.na(nc) || nc < 1L)
    nc <- tryCatch(parallel::detectCores(), error = function(e) NA_integer_)
  if (is.na(nc) || nc < 1L) nc <- 1L
  max(1L, min(as.integer(cap), nc - 1L))
}

# Re-slice a pilot's rhythm_fit at a new threshold (mirrors SCP:::.reslice_pilot,
# inlined so it also works on raw readRDS in demo mode). rhythm_fit is p-sorted,
# so this is a prefix slice; the top-300 drive the effect-size distribution.
# Legacy pilots without rhythm_fit pass through unchanged.
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
  # Full gene count as denominator, immune to p$ngenes later holding the sim count.
  denom <- p$rhythm_denom %||% p$ngenes %||% nrow(rf)
  p$rhythm_denom   <- denom
  p$prop_rhythmic  <- nrow(cand) / denom
  if ("A2" %in% names(rf)) {
    # Two-harmonic pilot: keep the five fields gene-aligned.
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
  # paired_sigma=TRUE keeps (A, sigma) paired (narrow r-tilde, Panel A); FALSE
  # shuffles sigma for the wide stratified range (Panels B/C). Default: paired.
  if (!isTRUE(paired_sigma) && length(p$sigma_rhythmic) > 1L)
    p$sigma_rhythmic <- sample(p$sigma_rhythmic)
  p$paired_sigma <- isTRUE(paired_sigma)
  p
}

# Ensembl-style ID -> gene symbol per organism, used when an upload's IDs match
# `pat` and the annotation package is installed. Unmapped IDs keep their name.
.SPECIES_MAP <- list(
  human       = list(pkg = "org.Hs.eg.db",    key = "ENSEMBL", col = "SYMBOL",   pat = "^ENSG[0-9]+",       lab = "human (Ensembl)"),
  mouse       = list(pkg = "org.Mm.eg.db",    key = "ENSEMBL", col = "SYMBOL",   pat = "^ENSMUSG[0-9]+",    lab = "mouse (Ensembl)"),
  rat         = list(pkg = "org.Rn.eg.db",    key = "ENSEMBL", col = "SYMBOL",   pat = "^ENSRNOG[0-9]+",    lab = "rat (Ensembl)"),
  zebrafish   = list(pkg = "org.Dr.eg.db",    key = "ENSEMBL", col = "SYMBOL",   pat = "^ENSDARG[0-9]+",    lab = "zebrafish (Ensembl)"),
  fly         = list(pkg = "org.Dm.eg.db",    key = "FLYBASE", col = "SYMBOL",   pat = "^FBgn[0-9]+",       lab = "Drosophila (FlyBase)"),
  worm        = list(pkg = "org.Ce.eg.db",    key = "ENSEMBL", col = "SYMBOL",   pat = "^WBGene[0-9]+",     lab = "C. elegans (WormBase)"),
  yeast       = list(pkg = "org.Sc.sgd.db",   key = "ORF",     col = "GENENAME", pat = "^Y[A-P][LR][0-9]{3}[WC]", lab = "S. cerevisiae (SGD)"),
  malaria     = list(pkg = "org.Pf.plasmo.db", key = "ORF",    col = "SYMBOL",   pat = "^PF3D7_[0-9]",      lab = "P. falciparum (PlasmoDB)"),
  arabidopsis = list(pkg = "org.At.tair.db",  key = "TAIR",    col = "SYMBOL",   pat = "^AT[0-9CM]G[0-9]+", lab = "Arabidopsis (TAIR)")
)

.upload_symbols <- function(rf, do_map = TRUE, species = "auto") {
  if (is.null(rf) || !nrow(rf) || !("gene" %in% names(rf)))
    return(list(rf = rf, type = "unknown", detected = NA_character_, n_mapped = 0L))
  g <- as.character(rf$gene); sym <- g
  # Detect the dominant ID pattern regardless of the user's choice, so we can
  # suggest a species.
  hits     <- vapply(.SPECIES_MAP, function(s) mean(grepl(s$pat, g)), numeric(1))
  dom      <- which.max(hits)
  detected <- if (length(dom) && hits[dom] > 0.5) names(.SPECIES_MAP)[dom] else NA_character_
  map_with <- function(idx, s) {
    if (!any(idx) || !requireNamespace(s$pkg, quietly = TRUE)) return(invisible())
    db <- tryCatch(getExportedValue(s$pkg, s$pkg), error = function(e) NULL)
    if (is.null(db)) return(invisible())
    keys <- sub("\\..*", "", g[idx])
    m <- suppressMessages(tryCatch(
      AnnotationDbi::mapIds(db, keys = keys, column = s$col %||% "SYMBOL",
                            keytype = s$key, multiVals = "first"),
      error = function(e) NULL))
    if (!is.null(m)) { hit <- !is.na(m); sym[which(idx)[hit]] <<- m[hit] }
  }
  use <- if (identical(species, "auto")) detected else species
  if (!is.null(use) && use %in% names(.SPECIES_MAP)) {
    s <- .SPECIES_MAP[[use]]; type <- s$lab
    if (isTRUE(do_map)) map_with(rep(TRUE, length(g)), s)
  } else {
    type <- "gene symbols / names"
  }
  rf$symbol <- sym
  list(rf = rf, type = type, detected = detected,
       n_mapped = sum(sym != g))
}

# ---------------------------------------------------------------------------
# UI
# ---------------------------------------------------------------------------

ui <- fluidPage(
  tags$head(tags$style(HTML("
    .center-titles h1, .center-titles h2, .center-titles h3, .center-titles h4, .center-titles p {
      text-align: center;
    }
    /* center the figures in the right pane */
    #power_curve img { display: block; margin-left: auto; margin-right: auto; }
    #clock_panel img, #gene_cosinor img, #enrich_plot img, #tod_plot img { display: block; margin-left: auto; margin-right: auto; }
    /* compact gene table so the K=2 columns fit without wrapping */
    #gene_table table { font-size: 0.84em; margin-bottom: 0; }
    #gene_table th, #gene_table td { padding: 3px 7px; white-space: nowrap; }
  "))),
  div(class = "center-titles",
      titlePanel("Power Evaluation and Study Design for Circadian Biomarker Detection"),
      p(style = "font-size: 0.95em; color:#444; margin-bottom: 4px;",
        "Thien Quy Pham, Hao Wang, Marianne L. Seney, Colleen A. McClung, Guan Yu, and George C. Tseng"),
      p(em("Pilot-calibrated power analysis for circadian transcriptomics studies."),
        br(),
        "Pick a bundled pilot, set the target power and FDR, and read off the recommended sample size.")
  ),
  sidebarLayout(
    sidebarPanel(
      width = 3,
      h3("Circadian Power Study", style = "text-align:center; margin-top:2px; font-weight:700;"),
      hr(style = "margin-top:6px; margin-bottom:10px;"),
      h4("Pilot rhythmicity threshold"),
      radioButtons("rhy_stat", NULL,
                   choices = c("Adjusted FDR (BH)" = "q",
                               "Raw p-value"       = "p"),
                   selected = "q", inline = TRUE),
      selectInput("rhy_thresh", HTML("Threshold (&alpha;<sub>pilot</sub>)"),
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
        selectInput("design_filter", "Design",     choices = c("Active" = "active", "Passive" = "passive")),
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
        selectInput("upload_species", "Species of uploaded data",
                    choices = c("Auto-detect" = "auto", "Human" = "human",
                                "Mouse" = "mouse", "Rat" = "rat",
                                "Zebrafish" = "zebrafish", "Drosophila (fly)" = "fly",
                                "C. elegans (worm)" = "worm",
                                "S. cerevisiae (yeast)" = "yeast",
                                "P. falciparum" = "malaria",
                                "Arabidopsis" = "arabidopsis", "Other / none" = "none"),
                    selected = "auto"),
        checkboxInput("upload_map_sym",
                      "Map gene IDs to gene symbols (unmapped keep their original name)",
                      value = TRUE),
        helpText(em("Symbol mapping supported for: human, mouse, rat, zebrafish, ",
                    "Drosophila, C. elegans, S. cerevisiae, P. falciparum, Arabidopsis ",
                    "(the Metascape species with Bioconductor annotation). ",
                    "The app also auto-detects the species from your IDs and suggests one. ",
                    "Any other organism, gene symbols, or probe IDs are used as-is.")),
        helpText(em("Once both files load, the pilot is fit on the fly (5-15 sec). ",
                    "Upload limit 5 GB, so a full gene matrix can be used as-is. ",
                    "See the expected file layout on the right.")),
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
      sliderInput("n_min", HTML("N<sub>min</sub> (minimum sample size)"), min = 4, max = 100,
                  value = 20, step = 2),
      sliderInput("n_max", HTML("N<sub>max</sub> (maximum sample size)"), min = 20, max = 500,
                  value = 120, step = 10),
      sliderInput("n_step", "Step", min = 2, max = 60,
                  value = 20, step = 2),
      helpText(HTML("<em>Grid: N<sub>min</sub>, N<sub>min</sub> + step, ..., up to N<sub>max</sub>.</em>")),
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
      uiOutput("stale_banner"),
      hr(),
      actionButton("reset_app", "Reset app", class = "btn-block", width = "100%",
                   icon = icon("rotate-left")),
      helpText(em("Reload everything to the starting state.")),
      actionButton("capture_app", "Capture full screen (PNG)", class = "btn-block",
                   width = "100%", icon = icon("camera")),
      helpText(em("Save the entire app view as one image."))
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
          h4("Sampling design (TOD distribution)", style = "text-align:center;"),
          plotOutput("tod_plot", height = "auto", width = "auto")
        )
      ),
      # ---- Rhythmic gene explorer (pilot-level; shown before & after a run) ----
      hr(),
      div(id = "gene_explore_card",
          style = "background:#fbfcfe; border:1px solid #e3e9f0; border-radius:6px; padding:14px 16px; margin-top:6px;",
        div(style = "display:flex; justify-content:space-between; align-items:center;",
            h4("Circadian Biomarker Detection", style = "margin:0;"),
            actionButton("capture_explorer", "Capture this panel (PNG)",
                         icon = icon("camera"), class = "btn-sm")),
        helpText(em(paste0("Per-gene rhythmicity at the threshold set on the left (changing it ",
                           "updates the clock panel and table below). Curves are the fitted cosinor; ",
                           "the shaded band is ± 1.96σ (noise). Bundled pilots store fit ",
                           "estimates only, so they show the curve; uploaded data also overlays each ",
                           "sample's actual expression as points."))),
        # Banner echoing the current threshold + rhythmic count.
        div(style = "padding:6px 10px; background:#eef4fb; border-radius:4px; margin-bottom:8px;",
            uiOutput("explorer_thresh")),
        uiOutput("k2_status"),
        tags$b("Core clock genes"),
        plotOutput("clock_panel", height = "auto", width = "auto"),
        uiOutput("clock_note"),
        hr(),
        fluidRow(
          column(7,
            tags$b("Top rhythmic genes"),
            div(style = "max-height: 340px; overflow-y: auto; margin-top:6px;",
                tableOutput("gene_table")),
            div(style = "margin-top:6px;",
                helpText(em("Downloads include ALL fitted genes with full parameters (p, BH q, amplitude, sigma, r-tilde, peak, mesor).")),
                downloadButton("dl_genes",      "Download all genes (CSV)"),
                downloadButton("dl_genes_xlsx", "Download all genes (XLSX)"))
          ),
          column(5,
            tags$b("Gene detail"),
            selectizeInput("gene_pick", "Look up a gene (symbol or ID)",
                           choices = NULL, multiple = FALSE,
                           options = list(placeholder = "type a gene symbol...")),
            plotOutput("gene_cosinor", height = "260px"),
            verbatimTextOutput("gene_readout")
          )
        )
      ),
      # ---- Pathway enrichment (clusterProfiler hypergeometric + BH) ----
      div(style = "background:#fbfcfe; border:1px solid #e3e9f0; border-radius:6px; padding:14px 16px; margin-top:10px;",
        h4("Pathway enrichment", style = "margin-top:0;"),
        helpText(em(paste0("Over-representation of the rhythmic gene set (genes passing the threshold on ",
                           "the left) against the chosen ontology, via clusterProfiler (hypergeometric + ",
                           "Benjamini-Hochberg). Top significant terms shown as -log10(P). KEGG requires internet; GO/Reactome run offline."))),
        fluidRow(
          column(3,
            selectInput("enrich_db", "Ontology",
                        choices = c("KEGG pathways" = "kegg", "Reactome" = "reactome",
                                    "GO Biological Process" = "gobp",
                                    "GO Molecular Function" = "gomf",
                                    "GO Cellular Component" = "gocc"),
                        selected = "kegg")),
          column(3,
            selectInput("enrich_species", "Gene-set species",
                        choices = c("Auto (from pilot)" = "auto", "Human" = "human", "Mouse" = "mouse"),
                        selected = "auto")),
          column(3,
            radioButtons("enrich_bg", "Background",
                         choices = c("Tested genes (this pilot)" = "tested", "All genes (genome-wide)" = "all"),
                         selected = "tested")),
          column(3,
            selectInput("enrich_top", "Show top", choices = c(10, 15, 20), selected = 15))
        ),
        helpText(em(paste0("Background: 'Tested genes' restricts the universe to genes measured in this ",
                           "pilot (or your uploaded matrix), the rigorous choice used in the paper; ",
                           "'All genes' is the genome-wide default and gives smaller p-values."))),
        fluidRow(
          column(4,
            selectInput("enrich_sig", "Significance cutoff",
                        choices = c("Adjusted P < 0.05" = "adj_0.05", "Adjusted P < 0.01" = "adj_0.01",
                                    "P < 0.05" = "p_0.05", "P < 0.01" = "p_0.01"),
                        selected = "adj_0.05"))
        ),
        actionButton("run_enrich", "Run enrichment", class = "btn-primary"),
        plotOutput("enrich_plot", height = "auto", width = "auto"),
        uiOutput("enrich_note"),
        div(style = "margin-top:4px;",
            downloadButton("dl_enrich_top", "Download displayed terms (CSV)"),
            downloadButton("dl_enrich_all", "Download all significant (CSV)"))
      ),
      conditionalPanel(
        condition = "input.run > 0",
        hr(),
        h4("Power curve"),
        plotOutput("power_curve", height = "auto"),
        div(style = "margin-top: 6px;",
            downloadButton("dl_pdf", "Download figure (PDF)"),
            downloadButton("dl_png", "Download figure (PNG)")),
        hr(),
        div(style = "padding: 12px 14px; background:#f4f8fc; border-left: 4px solid #2c7fb8; border-radius: 0 4px 4px 0;",
            uiOutput("recommended_n_text"))
      ),
      conditionalPanel(
        condition = "input.run == 0",
        br(),
        helpText(em("Click 'Run simulation' to compute the power curve."))
      ),
      # Upload format guide, shown in the right panel before a run.
      conditionalPanel(
        condition = "input.pilot_source == 'upload' && input.run == 0",
        hr(),
        h4("Expected file layout"),
        fluidRow(
          column(6,
            tags$b("1. Expression CSV"),
            tags$div(tags$small("Genes in rows, samples in columns. First column = gene IDs, header row = sample IDs.")),
            tags$pre(style = "font-size:12px;",
                     ",Sample1,Sample2,Sample3\nGene1,4.21,3.88,5.10\nGene2,1.07,0.92,1.31\nGene3,2.55,2.71,2.40")
          ),
          column(6,
            tags$b("2. Times-of-day CSV"),
            tags$div(tags$small("One time-of-day (hours) per sample, no header. Must be in the same order as the expression columns.")),
            tags$pre(style = "font-size:12px;", "0.5\n6.0\n11.5")
          )
        ),
        helpText(em("A ready-made example ships with the package: ",
                    tags$code('system.file("extdata/example", package = "SCP")'), "."))
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
    updateSelectInput(session, "species", choices = .cap_choices(sp),
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
    updateSelectInput(session, "design_filter", choices = .cap_choices(designs),
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
    updateSelectInput(session, "condition", choices = .cap_choices(cd), selected = cd[1])
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
        # Full CircadianBioOptions (amplitude/sigma/phase + rhythm_fit), so the
        # upload drives both the simulation and the threshold knob.
        bio <- SCP::estCircadianParam(mat, tod, paired_sigma = TRUE, verbose = FALSE)
        incProgress(0.7)
        bio$times <- tod
        # Keep the matrix in-session so the gene plot can overlay sample points.
        bio$.raw_expr <- mat
        sm <- .upload_symbols(bio$rhythm_fit, do_map = isTRUE(input$upload_map_sym),
                              species = input$upload_species %||% "auto")
        bio$rhythm_fit <- sm$rf
        attr(bio, "id_type")  <- sm$type
        attr(bio, "detected") <- sm$detected
        attr(bio, "n_mapped") <- sm$n_mapped
        # K=2 variant so the two-harmonic detector has real (A2, phi2) truth to
        # simulate from. pilot() picks it up via $.pilot2 when K=2; NULL when the
        # fit can't be made (N < 6, or no gene passes the pre-screen).
        bio$.pilot2 <- tryCatch({
          b2 <- SCP::estCircadianParam2H(mat, tod, paired_sigma = TRUE, verbose = FALSE)
          b2$times     <- tod
          b2$.raw_expr <- mat
          b2$rhythm_fit <- .upload_symbols(b2$rhythm_fit,
                                           do_map  = isTRUE(input$upload_map_sym),
                                           species = input$upload_species %||% "auto")$rf
          b2
        }, error = function(e) NULL)
        incProgress(1)
        bio
      })
    }, error = function(e) e)

    if (inherits(res, "error")) {
      uploaded_pilot_state$pilot <- NULL
      uploaded_pilot_state$msg <- sprintf("Upload error: %s", conditionMessage(res))
    } else {
      uploaded_pilot_state$pilot <- res
      rt_up    <- (res$amplitude %||% NA_real_) / (res$sigma_rhythmic %||% NA_real_)
      n_mapped <- attr(res, "n_mapped") %||% 0L
      detected <- attr(res, "detected")
      sel      <- input$upload_species %||% "auto"
      labof    <- function(k) if (!is.na(k) && k %in% names(.SPECIES_MAP)) .SPECIES_MAP[[k]]$lab else NA
      # Suggest the auto-detected species, flagged when it differs from the choice.
      sugg <- if (is.na(detected)) "Auto-detect: IDs look like gene symbols/names already."
              else if (sel == "auto") sprintf("Auto-detected %s.", labof(detected))
              else if (sel == detected) sprintf("Matches detected species (%s).", labof(detected))
              else sprintf("Note: you selected %s, but the IDs look like %s. Switch the dropdown if that is wrong.",
                           .SPECIES_MAP[[sel]]$lab %||% sel, labof(detected))
      map_note <- if (!isTRUE(input$upload_map_sym)) "Mapping off (showing original IDs)."
                  else if (n_mapped > 0) sprintf("Mapped %d IDs to symbols.", n_mapped)
                  else "No IDs needed mapping (already symbols / no annotation match)."
      k2_note  <- if (!is.null(res$.pilot2))
                    "Two-harmonic (K=2) fit ready: switch the detector to K=2 to use it."
                  else "Two-harmonic (K=2) fit unavailable (need N>=6 and rhythmic genes); K=2 falls back to K=1."
      uploaded_pilot_state$msg <- sprintf(
        "OK: %d genes fit, %d rhythmic (%.0f%%). Median r-tilde = %.2f.\n%s\n%s\n%s",
        res$ngenes %||% NA_integer_,
        length(res$amplitude %||% numeric(0)),
        100 * (res$prop_rhythmic %||% NA_real_),
        stats::median(rt_up, na.rm = TRUE), sugg, map_note, k2_note)
    }
  })

  uploaded_pilot <- reactive({ uploaded_pilot_state$pilot })

  output$upload_status <- renderText({ uploaded_pilot_state$msg })

  # Map the canonical display tissue back to the raw manifest value. Don't filter
  # on condition here: it can lag during the dropdown cascade.
  .resolve_actual_tissue <- function() {
    m <- manifest()
    sub <- m[m$species  == input$species  &
             m$design   == input$design_filter &
             m$dataset  == input$dataset, , drop = FALSE]
    sub <- sub[.tissue_display(sub) == input$tissue, , drop = FALSE]
    if (nrow(sub) > 0L) sub$tissue[1] else input$tissue
  }

  # ---- bundled pilot ------------------------------------------------------
  pilot <- reactive({
    if (isTRUE(input$pilot_source == "upload")) {
      # Uploaded pilots carry a rhythm_fit table too, so honor the threshold
      # slider just like bundled pilots. For K=2, use the embedded two-harmonic
      # fit (built at upload) so the detector simulates from real (A2, phi2)
      # truth; fall back to the K=1 fit if the 2H fit was unavailable.
      up <- uploaded_pilot()
      want_k2 <- isTRUE(input$K == "2")
      base <- if (want_k2 && !is.null(up$.pilot2)) up$.pilot2 else up
      # Flag when K=2 was requested but the uploaded design can't identify it.
      if (want_k2 && is.null(up$.pilot2) && !is.null(base))
        attr(base, "k2_unavailable") <- TRUE
      return(.app_apply_threshold(base, input$rhy_stat,
                                  as.numeric(input$rhy_thresh)))
    }
    req(input$species, input$dataset, input$tissue, input$condition)
    actual_tissue <- .resolve_actual_tissue()
    K_req <- as.integer(input$K %||% 1L)
    # Load only when the full combo exists. The dropdown cascade briefly leaves
    # condition/dataset lagging tissue; wait it out rather than erroring.
    m0 <- manifest()
    combo_ok <- nrow(m0[m0$species == input$species & m0$dataset == input$dataset &
                        m0$tissue == actual_tissue & m0$condition == input$condition, ,
                        drop = FALSE]) > 0L
    req(combo_ok)
    p <- tryCatch(
      .read_pilot(input$species, input$dataset,
                  actual_tissue, input$condition, K = K_req),
      error = function(e) {
        if (K_req == 2L) {
          # No K=2 variant for this pilot: fall back to K=1 so the panel still
          # works; the k2_status note explains why instead of a blank screen.
          fb <- tryCatch(.read_pilot(input$species, input$dataset,
                                     actual_tissue, input$condition, K = 1L),
                         error = function(e2) NULL)
          if (!is.null(fb)) attr(fb, "k2_unavailable") <- TRUE
          fb
        } else {
          showNotification(
            sprintf("Could not load pilot [%s / %s / %s / %s]: %s",
                    input$species, input$dataset, actual_tissue, input$condition,
                    conditionMessage(e)),
            type = "error", duration = 12)
          NULL
        }
      }
    )
    if (!is.null(p) && !inherits(p, "CircadianBioOptions"))
      class(p) <- c("CircadianBioOptions", class(p))
    # Re-slice at the chosen threshold; on failure keep the loaded pilot rather
    # than blanking the panel.
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

  # Snapshot of the last run, to warn when settings drift afterward.
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
    # Uploaded pilots have no manifest row; summarise straight from the object.
    if (isTRUE(input$pilot_source == "upload")) {
      amp <- p$amplitude %||% numeric(0); sig <- p$sigma_rhythmic %||% numeric(0)
      k <- min(length(amp), length(sig))
      rv <- if (k > 0L) { v <- amp[seq_len(k)] / sig[seq_len(k)]; v[is.finite(v) & v > 0] } else numeric(0)
      nsamp <- length(p$times %||% p$raw$times %||% numeric(0))
      rlab  <- if (input$rhy_stat == "q") sprintf("FDR %.0f%%", 100*as.numeric(input$rhy_thresh))
               else sprintf("p < %.3g", as.numeric(input$rhy_thresh))
      return(sprintf(paste0(
        "Pilot     : User-uploaded\n",
        "Samples (n) : %s\nGenes     : %s\n",
        "Prop. rhythmic (%s) : %s\n",
        "Median r-tilde (A / sigma) : %s   IQR [%s, %s]"),
        if (nsamp > 0L) nsamp else "NA",
        format(p$ngenes %||% NA_integer_, big.mark = ","),
        rlab,
        if (is.finite(p$prop_rhythmic %||% NA_real_)) sprintf("%.1f%%", 100*p$prop_rhythmic) else "NA",
        if (length(rv) >= 3L) sprintf("%.2f", median(rv)) else "NA",
        if (length(rv) >= 3L) sprintf("%.2f", quantile(rv, .25)) else "NA",
        if (length(rv) >= 3L) sprintf("%.2f", quantile(rv, .75)) else "NA"))
    }
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
    # If no gene passed the threshold (weak tissues at FDR 5%), summarise r-tilde
    # from the full rhythm_fit pool so it isn't NA. prop_rhythmic stays honest.
    if (!is.finite(r_med) && !is.null(p$rhythm_fit) &&
        nrow(p$rhythm_fit) > 0L && !is.null(p$rhythm_fit$A)) {
      rfr <- p$rhythm_fit$A / p$rhythm_fit$sigma
      rfr <- rfr[is.finite(rfr) & rfr > 0]
      if (length(rfr) > 0L) {
        r_med <- median(rfr); r_q25 <- quantile(rfr, 0.25); r_q75 <- quantile(rfr, 0.75)
      }
    }
    # prop_rhythmic was already set at the chosen threshold in pilot(); for legacy
    # pilots without rhythm_fit, recompute from the raw p-vector.
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
        # Prefer the pilot's actual stored sample times; fall back to manifest n.
        nn <- length(p$times %||% p$raw$times %||% numeric(0))
        if (nn < 1L) {
          mrow <- manifest()[manifest()$species == input$species &
                             manifest()$dataset == input$dataset &
                             .tissue_display(manifest()) == input$tissue   &
                             manifest()$condition == input$condition, , drop = FALSE]
          nn <- if (nrow(mrow) > 0L && "n" %in% names(mrow)) mrow$n[1] else NA_integer_
        }
        if (is.null(nn) || is.na(nn) || nn < 1L) "NA" else format(nn, big.mark = ",")
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
  }, width = function() 420L, height = 200L)   # fixed size + centered via CSS

  sim_result <- eventReactive(input$run, {
    p <- pilot()
    req(p)
    withProgress(message = "Running simulation...", value = 0.1, {
      n_cores <- .safe_cores(4L)
      incProgress(0.2, detail = "Setting up design")

      # Coerce numeric fields (some pilots store them as character/list).
      .as_num <- function(x) suppressWarnings(as.numeric(x))
      for (fld in c("amplitude", "sigma_rhythmic", "phase",
                    "lBaselineExpr", "lOD")) {
        if (!is.null(p[[fld]])) p[[fld]] <- .as_num(p[[fld]])
      }
      # 2000 sim genes, as in the manuscript scripts; the runner resamples the
      # per-gene vectors to this length.
      p$ngenes <- 2000L
      if (!is.null(p$period)) p$period <- .as_num(p$period)
      # Redefine the rhythmic set at the sidebar threshold (top-300 by p among
      # passing genes), so the threshold drives the simulated effect-size pool.
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
      # Match amplitude/phase length to sigma_rhythmic (some pilots were built
      # with mismatched gene sets, e.g. amplitude 1122 vs sigma 300).
      if (!is.null(p$amplitude) && !is.null(p$sigma_rhythmic)) {
        n_sigma <- length(p$sigma_rhythmic)
        if (length(p$amplitude) != n_sigma && n_sigma > 0L) {
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
      # Build sample-size grid from user-chosen min/max/step. Guard against an
      # inverted range (n_min > n_max), which would make seq() throw "wrong sign
      # in 'by'": swap so the grid is always ascending, and force a positive step.
      n_lo   <- as.integer(input$n_min)
      n_hi   <- as.integer(input$n_max)
      n_by   <- max(1L, as.integer(input$n_step))
      if (n_lo > n_hi) { tmp <- n_lo; n_lo <- n_hi; n_hi <- tmp }
      n_grid <- seq(n_lo, n_hi, by = n_by)
      if (length(n_grid) < 2L) {
        n_grid <- n_lo + n_by * 0:3
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

      # Use the same seed as the manuscript figure scripts (GLOBAL_SEED = 2025)
      K_val <- as.integer(input$K)
      # Advance the bar per sample size (and per sim when sensitivity is on) so
      # it moves smoothly instead of sitting frozen during the long sim call.
      n_steps  <- length(n_grid) * (if (isTRUE(input$eff_sens)) 2L else 1L)
      step_amt <- 0.8 / max(n_steps, 1L)
      run_sim <- function(bio, label) tryCatch(
        { set.seed(2025)
          SCP::runSimsSingleCohort(
            bio.opts = bio, design.opts = design, analysis.opts = analysis,
            K = K_val, mc.cores = n_cores,
            progress = function(j, nj, n)
              incProgress(step_amt,
                          detail = sprintf("%s | n = %d  (%d of %d)", label, n, j, nj))) },
        error = function(e) {
          showNotification(sprintf("Simulation failed: %s", conditionMessage(e)),
                           type = "error", duration = 8)
          structure(list(error = conditionMessage(e)), class = "sim_error")
        }
      )
      # Panel A: paired (A, sigma) draw -> realistic marginal power (p is paired).
      res <- run_sim(p, if (isTRUE(input$eff_sens)) "Panel A" else "Simulating")
      # Effect-size sensitivity (Panels B/C): unpaired draw -> wide r-tilde, so
      # the stratified panels span the full effect-size range. Re-derived from
      # the same rhythm_fit table at the same threshold (no extra fitting).
      if (isTRUE(input$eff_sens) && !inherits(res, "sim_error")) {
        p_unpaired <- .app_apply_threshold(p, input$rhy_stat,
                                           as.numeric(input$rhy_thresh),
                                           paired_sigma = FALSE)
        res$.bc <- run_sim(p_unpaired, "Panels B & C")
      }
      incProgress(1, detail = "Done")
      last_run_state(isolate(current_state()))
      res
    })
  })

  # Shared drawing fn so the on-screen plot and the PDF/PNG downloads are identical.
  .draw_power <- function() {
    res <- sim_result()
    if (is.null(res)) { plot.new(); title(main = "Run a simulation first."); return(invisible()) }
    if (inherits(res, "sim_error")) { plot.new(); title(main = sprintf("Simulation error: %s", res$error %||% "unknown")); return(invisible()) }
    # A NaN power means zero detections in those replicates (0/0). Render it as 0
    # so plotting never hits a non-finite 'invalid graphics state' (which blanked
    # Panel A and crashed the device on weak/edge-case runs).
    .san <- function(x) {
      if (is.null(x) || inherits(x, "sim_error")) return(x)
      for (f in c("marginal_power","marginal_FDR","marginal_TD","marginal_FD",
                  "strat_power","strat_TD","strat_FD","strat_FDR"))
        if (!is.null(x[[f]])) x[[f]][!is.finite(x[[f]])] <- 0
      x
    }
    res <- .san(res); res$.bc <- .san(res$.bc)
    tfd <- as.numeric(input$target_fdr); tpw <- as.numeric(input$target_power)
    # Panel A always shows FDR 1/5/10/20% + the user's threshold (post-hoc, free);
    # the blue vline + recommended-N stay tied to the chosen FDR.
    fdr_lines <- sort(unique(c(0.01, 0.05, 0.10, 0.20, tfd)))
    bc <- res$.bc
    if (!is.null(bc) && !inherits(bc, "sim_error")) {
      # Vertical stack: each panel sized + styled like the single Panel A view,
      # so the three panels look identical in scale, just stacked.
      SCP::plotSingleCohortPower(bc, panel_a_res = res, panels = c("A", "B", "C"),
                                  vertical = TRUE,
                                  fdr_thresholds = fdr_lines, panel_fdr = tfd, vline_fdr = tfd,
                                  vline_power = tpw, r_max = NULL,
                                  cex_main = 1.9, cex_lab = 1.85, cex_axis = 1.55)
    } else {
      SCP::plotSingleCohortPower(res, panels = "A",
                                  fdr_thresholds = fdr_lines, panel_fdr = tfd, vline_fdr = tfd,
                                  vline_power = tpw, r_max = NULL,
                                  cex_main = 1.9, cex_lab = 1.85, cex_axis = 1.6)
    }
  }

  output$power_curve <- renderPlot({
    tryCatch(.draw_power(),
             error = function(e) { plot.new(); title(main = sprintf("Plot error: %s", conditionMessage(e))) })
  },
  # Both views are fixed-width and centered (via CSS) so they don't stretch
  # across the whole right pane on wide screens. The 3-panel view stacks the
  # panels vertically, each rendered large (Figure-3-like) rather than squished
  # side by side; single Panel A stays compact.
  width  = function() 680L,
  height = function() if (isTRUE(input$eff_sens)) 1450L else 470L)

  # Download the exact figure shown, as PDF or PNG.
  .fig_dims  <- function() if (isTRUE(input$eff_sens)) c(6.8, 14.5) else c(8, 5.5)
  .fig_fname <- function(ext) sprintf("SCP_power_%s.%s",
      gsub("[^A-Za-z0-9]+", "_", paste(input$species %||% "", input$tissue %||% "pilot", sep = "_")), ext)
  output$dl_pdf <- downloadHandler(
    filename = function() .fig_fname("pdf"),
    content  = function(file) {
      d <- .fig_dims(); grDevices::pdf(file, width = d[1], height = d[2]); on.exit(grDevices::dev.off())
      tryCatch(.draw_power(), error = function(e) { plot.new(); title(conditionMessage(e)) })
    })
  # Open a PNG device that works on headless servers. The default png() device
  # needs an X11 display (absent on shinyapps.io / most Linux deploys), so prefer
  # the cairo backend, then ragg, then fall back to the stock device.
  .open_png <- function(file, width, height, res = 120) {
    if (isTRUE(capabilities("cairo"))) {
      grDevices::png(file, width = width, height = height, res = res, type = "cairo")
    } else if (requireNamespace("ragg", quietly = TRUE)) {
      ragg::agg_png(file, width = width, height = height, units = "px", res = res)
    } else {
      grDevices::png(file, width = width, height = height, res = res)
    }
  }
  output$dl_png <- downloadHandler(
    filename = function() .fig_fname("png"),
    content  = function(file) {
      d <- .fig_dims(); .open_png(file, d[1]*110, d[2]*110, res = 120); on.exit(grDevices::dev.off())
      tryCatch(.draw_power(), error = function(e) { plot.new(); title(conditionMessage(e)) })
    })

  # Capture the ENTIRE app UI (inputs, summary, plot, everything you would
  # scroll to) into one PNG via html2canvas, so it can be emailed as a single
  # snapshot. This works regardless of viewport/scroll, unlike a browser print.
  observeEvent(input$capture_app, {
    if (!requireNamespace("shinyscreenshot", quietly = TRUE)) {
      showNotification(
        "Full-app capture needs the 'shinyscreenshot' package. Install it with install.packages('shinyscreenshot') and relaunch. (You can still use the PDF/PNG figure downloads.)",
        type = "warning", duration = 12)
      return(invisible(NULL))
    }
    shinyscreenshot::screenshot(filename = "SCP_app_view", timer = 0)
  })

  # Capture just the gene-explorer region (clock panel + table + cosinor) as one
  # PNG, mirroring the power-curve capture but scoped to this card.
  observeEvent(input$capture_explorer, {
    if (!requireNamespace("shinyscreenshot", quietly = TRUE)) {
      showNotification(
        "Capture needs the 'shinyscreenshot' package: install.packages('shinyscreenshot') and relaunch.",
        type = "warning", duration = 10)
      return(invisible(NULL))
    }
    shinyscreenshot::screenshot(id = "gene_explore_card",
                                filename = "SCP_gene_explorer", timer = 0)
  })


  output$recommended_n_text <- renderUI({
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
      return(div(style = "color:#a33; font-size: 0.95em;",
                 "Simulation did not return a marginal_power matrix."))
    }
    # Typographic hierarchy: a bold headline number, then a smaller caption,
    # so the box reads cleanly instead of one long oversized sentence.
    head_style <- "font-size: 1.55em; font-weight: 600; color:#2c7fb8; line-height: 1.2;"
    sub_style  <- "font-size: 0.92em; color:#556; margin-top: 3px; line-height: 1.3;"
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
      return(tagList(
        div(style = head_style, sprintf("Recommended N = %d", rec_n)),
        div(style = sub_style,
            sprintf("for %.0f%% power at FDR %.0f%%", 100 * tpw, 100 * tfd),
            tags$span(style = "color:#889; font-style: italic;", " (from simulation grid)"))
      ))
    }
    n_largest <- max(sizes)
    last_pow  <- mean_pow[which(sizes == n_largest)[1]]
    tagList(
      div(style = head_style, sprintf("N > %d needed", n_largest)),
      div(style = sub_style,
          sprintf("Largest N simulated (%d) reached %.0f%% power at FDR %.0f%%; a larger N is needed to hit the %.0f%% target.",
                  n_largest, 100 * last_pow, 100 * tfd, 100 * tpw))
    )
  })

  # ======================= Rhythmic gene explorer =========================
  # Core clock genes with common aliases (datasets use either the official HGNC
  # symbol or the common name, e.g. ARNTL vs BMAL1). Matching is alias-aware so a
  # gene is recognised whichever name the pilot/upload uses. All upper-case.
  .CLOCK_ALIASES <- list(
    ARNTL  = c("ARNTL","BMAL1","MOP3"),
    ARNTL2 = c("ARNTL2","BMAL2","MOP9"),
    CLOCK  = c("CLOCK"),
    NPAS2  = c("NPAS2","MOP4"),
    PER1   = c("PER1"), PER2 = c("PER2"), PER3 = c("PER3"),
    CRY1   = c("CRY1"), CRY2 = c("CRY2"),
    NR1D1  = c("NR1D1","REVERBA","REV-ERBA","NR1D1A","THRAL"),
    NR1D2  = c("NR1D2","REVERBB","REV-ERBB","RVR"),
    DBP    = c("DBP"), TEF = c("TEF"), HLF = c("HLF"),
    CIART  = c("CIART","CHRONO","GM129","C1ORF51"),
    NFIL3  = c("NFIL3","E4BP4"),
    RORA   = c("RORA"), RORB = c("RORB"), RORC = c("RORC")
  )
  # alias (upper) -> canonical name, for vectorised lookup
  .CLOCK_LOOKUP <- unlist(lapply(names(.CLOCK_ALIASES), function(cn)
    stats::setNames(rep(cn, length(.CLOCK_ALIASES[[cn]])), .CLOCK_ALIASES[[cn]])))
  .clock_canon <- function(syms) unname(.CLOCK_LOOKUP[toupper(syms)])   # NA if not a clock gene

  # Per-gene table from the loaded pilot's rhythm_fit: BH q-value, r-tilde, peak,
  # mesor (if stored). Sorted by p-value. Reacts to the left-panel threshold via
  # pilot() (which carries the full capped rhythm_fit regardless of threshold).
  gene_tbl <- reactive({
    p <- pilot(); req(p); rf <- p$rhythm_fit
    if (is.null(rf) || !nrow(rf)) return(NULL)
    n_total <- p$rhythm_denom %||% p$ngenes %||% nrow(rf)
    # K=2 pilots carry p_joint (joint two-harmonic F-test, the manuscript's K=2
    # detector). Option A: rank + threshold the explorer on p_joint for K=2; fall
    # back to the cosinor p (pvalue = p_K1) for K=1 (and old un-enriched K=2).
    is_k2 <- !is.null(rf$p_joint) && any(is.finite(rf$p_joint))
    pp    <- if (is_k2) rf$p_joint else rf$pvalue
    ord   <- order(pp)
    rf <- rf[ord, , drop = FALSE]; pp <- pp[ord]
    n  <- nrow(rf); i <- seq_len(n)
    q  <- pmin(rev(cummin(rev(pp * n_total / i))), 1)
    # Robust to pilots whose rhythm_fit predates the gene/symbol/mesor columns:
    # fall back to a row index for the ID, the ID as the gene name, NA mesor.
    gene_id <- if (!is.null(rf$gene))   as.character(rf$gene)   else as.character(i)
    gene_nm <- if (!is.null(rf$symbol)) as.character(rf$symbol) else gene_id
    mesor_v <- if (!is.null(rf$mesor))  rf$mesor                else rep(NA_real_, n)
    A2_v    <- if (!is.null(rf$A2))     rf$A2                   else rep(NA_real_, n)
    phi2_v  <- if (!is.null(rf$phi2))   rf$phi2                 else rep(NA_real_, n)
    p2h_v   <- if (!is.null(rf$p_2h))   rf$p_2h                 else rep(NA_real_, n)
    data.frame(
      Gene   = gene_nm,
      ID     = gene_id,
      p      = pp, q = q,
      p_K1   = rf$pvalue,                 # 1st-harmonic F-test (the pre-screen)
      p_2h   = p2h_v,                     # 2nd-harmonic-only F-test
      Amp    = rf$A, Sigma = rf$sigma, rtilde = rf$A / rf$sigma,
      Peak_h = round(rf$phi %% 24, 2),    # phi1 = 1st-harmonic acrophase
      Mesor  = mesor_v,
      A2     = A2_v, phi2 = phi2_v, K2 = is_k2,
      stringsAsFactors = FALSE)
  })

  # K=2 identifiability notice: shown only when the detector is K=2. Confirms a
  # working two-harmonic fit, or explains the K=1 fallback when the sampling
  # design has too few distinct time points to identify the 12h second harmonic.
  output$k2_status <- renderUI({
    if (!isTRUE(input$K == "2")) return(NULL)
    p <- pilot(); if (is.null(p)) return(NULL)
    rf <- p$rhythm_fit
    nt <- length(unique(round((p$times %||% p$cts %||% numeric(0)) %% 24, 3)))
    has_k2 <- !is.null(rf$A2) && any(is.finite(rf$A2)) && !isTRUE(attr(p, "k2_unavailable"))
    if (has_k2)
      div(style = "padding:7px 12px; background:#e9f5ea; border-left:4px solid #4a9c5d; border-radius:0 4px 4px 0; margin-bottom:8px; font-size:0.9em;",
          HTML(sprintf("<b>Two-harmonic (K=2) fit active.</b> The sampling design identifies the second harmonic (%d distinct time points).", nt)))
    else
      div(style = "padding:7px 12px; background:#fff3cd; border-left:4px solid #d68f00; border-radius:0 4px 4px 0; margin-bottom:8px; font-size:0.9em;",
          HTML(sprintf("<b>Two-harmonic (K=2) not identifiable for this dataset.</b> Its sampling design (%d distinct time point%s) cannot separately estimate the 12-hour second harmonic (the 5-parameter fit is rank-deficient). The panel and simulation use the single-harmonic (K=1) fit; switch the detector to K=1 to hide this notice.",
                       nt, if (nt == 1L) "" else "s")))
  })

  # Genes passing the chosen alpha_pilot (matches the summary's prop_rhythmic).
  gene_tbl_pass <- reactive({
    g <- gene_tbl(); req(g)
    thr  <- as.numeric(input$rhy_thresh)
    sval <- if (identical(input$rhy_stat, "q")) g$q else g$p
    g[sval < thr, , drop = FALSE]
  })

  # Human-readable description of the current alpha_pilot threshold.
  .thresh_label <- reactive({
    metric <- if (identical(input$rhy_stat, "q")) "BH-FDR" else "raw p"
    sprintf("%s < %s", metric, input$rhy_thresh %||% "0.01")
  })

  # Live banner so it is obvious the explorer reflects the left-panel threshold.
  output$explorer_thresh <- renderUI({
    g <- gene_tbl()
    if (is.null(g)) return(HTML("<em>Select a pilot to explore its rhythmic genes.</em>"))
    np <- nrow(gene_tbl_pass())
    HTML(sprintf("<b>Threshold (&alpha;<sub>pilot</sub>):</b> %s &nbsp;|&nbsp; <b>%s</b> rhythmic genes pass &nbsp;|&nbsp; clock panel + table below update with this setting.",
                 .thresh_label(), format(np, big.mark = ",")))
  })

  # Populate the gene lookup (server-side; rhythmic gene set can be large).
  observe({
    g <- gene_tbl()
    ch <- if (is.null(g) || !nrow(g)) character(0)
          else stats::setNames(g$ID, sprintf("%s  (p=%.1e)", g$Gene, g$p))
    updateSelectizeInput(session, "gene_pick", choices = ch,
                         selected = if (length(ch)) ch[[1]] else "", server = TRUE)
  })

  # Times of the (up to two) peaks of the fitted curve. For K=1 this is just the
  # cosinor acrophase; for K=2 it returns the two local maxima of the
  # two-harmonic waveform (the "two peaks" a bimodal gene actually shows).
  .twoharm_peaks <- function(A1, phi1, A2 = NA, phi2 = NA, period = 24) {
    w <- 2 * pi / period
    t <- seq(0, period, length.out = 289)[-289]            # 0..<24 at 5-min grid
    y <- A1 * cos(w * (t - phi1))
    if (!is.null(A2) && is.finite(A2) && !is.null(phi2) && is.finite(phi2))
      y <- y + A2 * cos(2 * w * (t - phi2))
    yl <- c(y[length(y)], y[-length(y)]); yr <- c(y[-1], y[1])
    ism <- y > yl & y >= yr                                # circular local maxima
    pk <- t[ism][order(-y[ism])]
    c(peak1 = if (length(pk) >= 1L) pk[1] else (phi1 %% period),
      peak2 = if (length(pk) >= 2L) pk[2] else NA_real_)
  }

  # Fitted cosinor curve: mesor-anchored if stored, else mean-centered; +/-1.96s
  # noise ribbon; peak(s) marked. `compact` (clock small-multiples) shows only the
  # gene title + peak; the full p/q/peak/(A/sigma) subtitle is reserved for the
  # large single-gene detail plot where there is room (avoids title collisions).
  # `pts` (optional list(t, y)) overlays each sample's actual expression on the
  # fit -- available for uploaded data (raw matrix in session); bundled pilots
  # store no raw data so they show the fitted curve + noise band only.
  .draw_cos <- function(row, period = 24, compact = FALSE, pts = NULL) {
    t <- seq(0, period, length.out = 240); w <- 2 * pi / period
    m0 <- if (!is.null(row$Mesor) && is.finite(row$Mesor)) row$Mesor else 0
    y  <- m0 + row$Amp * cos(w * (t - row$Peak_h))
    # K=2 pilots: add the second harmonic so the curve shows the real (possibly
    # bimodal / asymmetric) two-harmonic waveform, not a pure sinusoid.
    if (!is.null(row$A2) && is.finite(row$A2) && !is.null(row$phi2) && is.finite(row$phi2))
      y <- y + row$A2 * cos(2 * w * (t - row$phi2))
    band <- 1.96 * row$Sigma
    ylab <- if (m0 != 0) "Expression (log2)" else "Expression (centered)"
    yl   <- range(c(y + band, y - band, if (!is.null(pts)) pts$y))
    plot(t, y, type = "n", ylim = yl,
         xlab = "Time of day (h)", ylab = if (compact) "centered" else ylab,
         main = row$Gene, xaxt = "n", bty = "l",
         cex.lab = if (compact) 1.0 else 1.2,
         cex.main = if (compact) 1.2 else 1.45, font.main = 2)
    axis(1, at = seq(0, 24, if (compact) 6 else 4))
    polygon(c(t, rev(t)), c(y + band, rev(y - band)), col = "#2c7fb822", border = NA)
    if (!is.null(pts)) points(pts$t %% period, pts$y, pch = 19,
                              col = "#2c7fb899", cex = 1.1)   # actual sample points
    lines(t, y, lwd = 2.5, col = "#2c7fb8")
    # mark peak(s): one dashed line for K=1, two for a bimodal K=2 waveform.
    pk <- .twoharm_peaks(row$Amp, row$Peak_h, row$A2, row$phi2, period)
    abline(v = pk[["peak1"]], lty = 2, col = "grey50")
    if (is.finite(pk[["peak2"]])) abline(v = pk[["peak2"]], lty = 2, col = "grey50")
    # stats subtitle on every panel; compact (clock small-multiples) drops q and
    # uses a smaller font so the line fits a narrow panel without overlapping.
    pk_txt <- if (is.finite(pk[["peak2"]]))
      sprintf("peaks %.1f, %.1fh", pk[["peak1"]], pk[["peak2"]])
    else sprintf("peak %.1fh", pk[["peak1"]])
    sub <- if (compact)
      sprintf("p=%.0e  %s  A/s=%.1f", row$p, pk_txt, row$rtilde)
    else
      sprintf("p=%.1e   q=%.1e   %s   A/s=%.2f", row$p, row$q, pk_txt, row$rtilde)
    mtext(sub, side = 3, line = 0.35, cex = if (compact) 0.78 else 0.92, col = "grey25")
  }

  # Clock genes passing the threshold (alias-aware, deduped). Select the top 12
  # by significance (so the strongest show when capped), then display them in
  # ALPHABETICAL order for a predictable, consistent layout.
  clock_genes <- reactive({
    g <- gene_tbl_pass(); if (is.null(g) || !nrow(g)) return(NULL)
    canon <- .clock_canon(g$Gene)
    pr <- g[!is.na(canon), , drop = FALSE]; pc <- canon[!is.na(canon)]
    if (!nrow(pr)) return(pr)
    pr <- pr[!duplicated(pc), , drop = FALSE]
    pr <- pr[order(pr$p), , drop = FALSE]           # rank by significance for the cap
    pr <- utils::head(pr, 12L)
    pr[order(toupper(pr$Gene)), , drop = FALSE]      # plot alphabetically
  })
  # Layout: <=4 columns; rows wrap. FIXED per-panel size so 1-2 genes render at a
  # consistent aspect instead of stretching across the full width.
  .CLOCK_PW <- 210L; .CLOCK_PH <- 188L
  .clock_grid <- function() {
    n <- tryCatch(nrow(clock_genes()), error = function(e) 0L)
    n <- if (is.null(n) || is.na(n)) 0L else n
    if (n < 1L) return(list(n = 0L, nc = 1L, nr = 1L))
    # Balanced columns so no count leaves a lonely/ragged bottom row (max 12):
    # 1-3 -> 1 row; 4 -> 2x2 block; 5-6 -> 3 cols; 7-8 -> 4 cols; 9 -> 3x3;
    # 10-12 -> 4 cols.
    nc <- if (n <= 3L) n
          else if (n == 4L) 2L
          else if (n <= 6L) 3L
          else if (n <= 8L) 4L
          else if (n == 9L) 3L
          else 4L
    list(n = n, nc = nc, nr = ceiling(n / nc))
  }
  # layout() matrix at half-cell resolution (2*nc sub-columns) so a partial last
  # row is padded with blank cells on both sides -> its panels sit CENTERED
  # (e.g. n=10 -> bottom row [blank, P, P, blank]) instead of left-aligned.
  .clock_layout_mat <- function(n, nc) {
    nr <- ceiling(n / nc); drawn <- 0L; rows <- vector("list", nr)
    for (r in seq_len(nr)) {
      k <- min(nc, n - drawn); ids <- seq.int(drawn + 1L, drawn + k)
      row <- if (k == nc) rep(ids, each = 2L)
             else { pad <- nc - k; c(rep(0L, pad), rep(ids, each = 2L), rep(0L, pad)) }
      rows[[r]] <- row; drawn <- drawn + k
    }
    do.call(rbind, rows)
  }

  output$clock_panel <- renderPlot({
    pr <- clock_genes()
    if (is.null(pr) || !nrow(pr)) {
      # No clock gene passes the threshold (common for a small uploaded gene
      # subset). Draw the message with zero plot margins: the fallback device
      # is short (see height fn) and R's DEFAULT margins exceed it, which would
      # otherwise throw "figure margins too large" instead of showing the note.
      op <- par(mar = c(0, 0, 0, 0)); on.exit(par(op))
      plot.new()
      text(0.5, 0.5, sprintf("No core clock genes pass %s in this pilot.",
                             .thresh_label()), cex = 1.1)
      return()
    }
    gd <- .clock_grid()
    layout(.clock_layout_mat(gd$n, gd$nc))
    op <- par(mar = c(3.3, 3.4, 3.4, 0.7), mgp = c(1.9, 0.6, 0))
    on.exit({ par(op); layout(1) })
    for (i in seq_len(gd$n)) .draw_cos(pr[i, ], compact = TRUE)
  },
  # device sized to the grid -> each panel keeps a fixed ratio for any gene count
  width  = function() { gd <- .clock_grid(); if (gd$n < 1L) 460L else gd$nc * .CLOCK_PW },
  height = function() { gd <- .clock_grid(); if (gd$n < 1L) 160L else gd$nr * .CLOCK_PH })

  output$clock_note <- renderUI({
    g <- gene_tbl_pass(); if (is.null(g)) return(NULL)
    present <- unique(.clock_canon(g$Gene)); present <- present[!is.na(present)]
    cg <- clock_genes()
    shown <- if (is.null(cg) || !nrow(cg)) character(0) else {
      s <- unique(.clock_canon(cg$Gene)); s[!is.na(s)] }
    absent <- setdiff(names(.CLOCK_ALIASES), present)   # not passing the threshold
    extra  <- setdiff(present, shown)                   # passing but beyond the 12-panel cap
    tagList(
      if (length(extra))
        helpText(em(sprintf("Other clock genes passing %s but not shown (panel caps at 12): %s",
                            .thresh_label(), paste(sort(extra), collapse = ", ")))),
      if (length(absent))
        helpText(em(sprintf("Clock genes not passing %s here: %s",
                            .thresh_label(), paste(absent, collapse = ", "))))
    )
  })

  output$gene_table <- renderTable({
    g <- gene_tbl_pass(); req(g)
    if (!nrow(g)) return(data.frame(Note = "No genes pass the chosen threshold."))
    out <- head(g, 100L)
    is_k2 <- isTRUE(any(out$K2)) && any(is.finite(out$p_2h))
    # Compact headers under K=2 so all columns fit one row (the K=2 context is
    # already shown by the banner above, so the headers can be short).
    if (is_k2) {
      pks <- mapply(.twoharm_peaks, out$Amp, out$Peak_h, out$A2, out$phi2)
      df <- data.frame(Gene = out$Gene,
                       a = formatC(out$p,   format = "e", digits = 1),
                       b = formatC(out$q,   format = "e", digits = 1),
                       c = formatC(out$p_2h, format = "e", digits = 1),
                       d = round(out$rtilde, 2),
                       e = round(pks["peak1", ], 1),
                       f = round(pks["peak2", ], 1),
                       check.names = FALSE, stringsAsFactors = FALSE)
      names(df) <- c("Gene", "p (joint)", "q", "p (2nd)", "r̃", "Peak1", "Peak2")
    } else {
      df <- data.frame(Gene = out$Gene,
                       a = formatC(out$p, format = "e", digits = 2),
                       b = formatC(out$q, format = "e", digits = 2),
                       c = round(out$rtilde, 2),
                       d = out$Peak_h,
                       check.names = FALSE, stringsAsFactors = FALSE)
      names(df) <- c("Gene", "p-value", "q-value", "r̃ (A/σ)", "Peak (h)")
    }
    df
  }, striped = TRUE, hover = TRUE, width = "100%", digits = 3)

  # Full export: ALL fitted genes (the pilot's complete rhythm_fit, p-sorted)
  # with every parameter and a clean column order/labels.
  gene_export <- reactive({
    g <- gene_tbl(); req(g)
    is_k2 <- isTRUE(any(g$K2)) && any(is.finite(g$p_2h))
    if (!is_k2) {
      # K=1: one harmonic, one peak (= the acrophase).
      data.frame(
        gene       = g$Gene, gene_id = g$ID,
        pvalue     = g$p, qvalue_BH = g$q,
        amplitude  = g$Amp, sigma = g$Sigma, r_tilde = g$rtilde,
        peak_hours = g$Peak_h, mesor = g$Mesor,
        stringsAsFactors = FALSE)
    } else {
      # K=2: ONE mesor; TWO harmonic amplitudes (A1,A2) + acrophases (phi1,phi2);
      # TWO derived waveform peaks (peak1,peak2 = the actual maxima). pvalue is
      # the joint F(4,N-5); p_1st/2nd_harmonic are the component tests.
      pks <- mapply(.twoharm_peaks, g$Amp, g$Peak_h, g$A2, g$phi2)
      data.frame(
        gene = g$Gene, gene_id = g$ID,
        pvalue_joint = g$p, qvalue_BH = g$q,
        p_1st_harmonic = g$p_K1, p_2nd_harmonic = g$p_2h,
        mesor = g$Mesor,
        A1 = g$Amp,  phi1_hours = g$Peak_h,
        A2 = g$A2,   phi2_hours = g$phi2,
        sigma = g$Sigma, r_tilde = g$rtilde,
        peak1_hours = round(pks["peak1", ], 3),
        peak2_hours = round(pks["peak2", ], 3),
        stringsAsFactors = FALSE)
    }
  })
  .gene_fname <- function() gsub("[^A-Za-z0-9]+", "_",
      paste(if (isTRUE(input$pilot_source == "upload")) "upload" else input$species %||% "pilot",
            input$dataset %||% "", input$tissue %||% ""))

  output$dl_genes <- downloadHandler(
    filename = function() sprintf("SCP_genes_%s.csv", .gene_fname()),
    content  = function(f) utils::write.csv(gene_export(), f, row.names = FALSE))

  output$dl_genes_xlsx <- downloadHandler(
    filename = function() sprintf("SCP_genes_%s.xlsx", .gene_fname()),
    content  = function(f) {
      g <- gene_export()
      if (requireNamespace("writexl", quietly = TRUE)) {
        writexl::write_xlsx(g, f)
      } else if (requireNamespace("openxlsx", quietly = TRUE)) {
        openxlsx::write.xlsx(g, f)
      } else {
        # graceful fallback: still deliver the data as CSV-in-xlsx-name
        utils::write.csv(g, f, row.names = FALSE)
        showNotification("XLSX writer not installed; wrote CSV content instead.",
                         type = "warning", duration = 6)
      }
    })

  output$gene_cosinor <- renderPlot({
    g <- gene_tbl(); req(g, input$gene_pick)
    row <- g[g$ID == input$gene_pick, , drop = FALSE]
    if (!nrow(row)) { plot.new(); title("Gene not rhythmic in this pilot (p >= 0.2)."); return() }
    # Overlay actual sample points when raw data is available (uploaded pilots).
    p <- pilot(); pts <- NULL
    rx <- p$.raw_expr
    if (!is.null(rx) && input$gene_pick %in% rownames(rx) && !is.null(p$times))
      pts <- list(t = as.numeric(p$times), y = as.numeric(rx[input$gene_pick, ]))
    .draw_cos(row[1, ], pts = pts)
  })

  output$gene_readout <- renderText({
    g <- gene_tbl(); req(g, input$gene_pick)
    row <- g[g$ID == input$gene_pick, , drop = FALSE]
    if (!nrow(row)) return("Gene not in the rhythmic set (p >= 0.2).")
    r <- row[1, ]
    is_k2 <- isTRUE(r$K2) && is.finite(r$A2)
    mesor_line <- if (is.finite(r$Mesor)) sprintf("%.3f\n", r$Mesor)
                  else "(not stored for this pilot)\n"
    if (is_k2) {
      # Full two-harmonic readout: both harmonics' (A, phi), both waveform peaks,
      # and all three p-values; r-tilde is the 1st-harmonic effect size A1/sigma.
      pk <- .twoharm_peaks(r$Amp, r$Peak_h, r$A2, r$phi2)
      paste0(
        sprintf("Gene       : %s  (%s)\n", r$Gene, r$ID),
        sprintf("p (joint)  : %.3e   <- two-harmonic detector\n", r$p),
        sprintf("q (BH)     : %.3e\n", r$q),
        sprintf("p (1st h.) : %.3e\n", r$p_K1),
        sprintf("p (2nd h.) : %.3e\n", r$p_2h),
        sprintf("A1 / phi1  : %.3f / %.2f h   (1st harmonic)\n", r$Amp, r$Peak_h),
        sprintf("A2 / phi2  : %.3f / %.2f h   (2nd harmonic)\n", r$A2, r$phi2),
        sprintf("Sigma      : %.3f\n", r$Sigma),
        sprintf("r-tilde    : %.2f   (A1/sigma)\n", r$rtilde),
        sprintf("Peak 1     : %.2f h\n", pk[["peak1"]]),
        if (is.finite(pk[["peak2"]])) sprintf("Peak 2     : %.2f h\n", pk[["peak2"]]) else "",
        sprintf("Mesor      : %s", mesor_line))
    } else {
      paste0(
        sprintf("Gene      : %s  (%s)\n", r$Gene, r$ID),
        sprintf("p-value   : %.3e\n", r$p),
        sprintf("q (BH)    : %.3e\n", r$q),
        sprintf("Amplitude : %.3f\n", r$Amp),
        sprintf("Sigma     : %.3f\n", r$Sigma),
        sprintf("r-tilde   : %.2f   (A/sigma)\n", r$rtilde),
        sprintf("Peak      : %.2f h\n", r$Peak_h),
        sprintf("Mesor     : %s", mesor_line))
    }
  })

  # Reset: reload the session to a clean blank state from anywhere (recovers a
  # stuck/frozen UI without restarting the R process).
  observeEvent(input$reset_app, { session$reload() })

  # ======================= Pathway enrichment (Enrichr) ===================
  # Human-readable label for the chosen ontology.
  .enrich_lib_label <- function(db, sp) {
    org <- if (identical(sp, "mouse")) "mouse" else "human"
    switch(db, kegg = sprintf("KEGG (%s)", org), reactome = sprintf("Reactome (%s)", org),
      gobp = "GO Biological Process", gomf = "GO Molecular Function",
      gocc = "GO Cellular Component", "KEGG")
  }
  # Resolve the gene-set species: explicit choice, else inferred from the pilot.
  .enrich_species <- function() {
    s <- input$enrich_species %||% "auto"
    if (s != "auto") return(s)
    if (isTRUE(input$pilot_source == "upload")) {
      up <- uploaded_pilot()
      d  <- if (!is.null(up)) attr(up, "detected") else NA
      if (!is.na(d) && d %in% c("human", "mouse")) return(d)
      if (identical(input$upload_species, "mouse")) return("mouse")
      return("human")
    }
    if (identical(input$species, "mouse")) "mouse" else "human"
  }
  # The tested-gene universe (background) for the loaded pilot, in its native ID
  # space: bundled pilots store $tested_genes; uploads use the raw matrix genes.
  .enrich_universe <- function() {
    p <- pilot(); if (is.null(p)) return(NULL)
    p$tested_genes %||% (if (!is.null(p$.raw_expr)) rownames(p$.raw_expr) else NULL)
  }

  # clusterProfiler over-representation (matches the paper: enrichKEGG/enrichGO/
  # enrichPathway with a custom universe). Runs on the button.
  enrich_res <- eventReactive(input$run_enrich, {
    g <- gene_tbl_pass()
    validate(need(!is.null(g) && nrow(g) > 0,
                  "No rhythmic genes at the current threshold to enrich."))
    sp    <- .enrich_species()
    orgdb <- if (sp == "mouse") "org.Mm.eg.db" else "org.Hs.eg.db"
    for (pk in c("clusterProfiler", "AnnotationDbi", orgdb))
      if (!requireNamespace(pk, quietly = TRUE))
        return(list(err = sprintf("Package '%s' is not installed (needed for enrichment).", pk)))
    org <- getExportedValue(orgdb, orgdb)
    to_entrez <- function(ids) {
      ids <- unique(ids[!is.na(ids) & nzchar(ids)])
      kt  <- if (mean(grepl("^ENS", ids)) > 0.5) "ENSEMBL" else "SYMBOL"
      keys <- sub("\\..*", "", ids)
      m <- suppressMessages(tryCatch(AnnotationDbi::mapIds(org, keys = keys, column = "ENTREZID",
                                     keytype = kt, multiVals = "first"), error = function(e) NULL))
      unique(m[!is.na(m)])
    }
    # background
    want_tested <- identical(input$enrich_bg %||% "tested", "tested")
    uni_ids <- if (want_tested) .enrich_universe() else NULL
    bg_used <- if (want_tested && !is.null(uni_ids)) "tested" else "all"
    fg  <- to_entrez(unique(g$Gene))
    if (!length(fg)) return(list(err = "No rhythmic genes could be mapped to Entrez IDs."))
    uni <- if (identical(bg_used, "tested")) to_entrez(uni_ids) else NULL
    db  <- input$enrich_db %||% "kegg"
    organism <- if (sp == "mouse") "mmu" else "hsa"
    ek <- withProgress(message = "Running enrichment...", value = 0.4, {
      tryCatch({
        r <- switch(db,
          kegg     = clusterProfiler::enrichKEGG(fg, organism = organism, universe = uni,
                                                 pvalueCutoff = 1, qvalueCutoff = 1, minGSSize = 3),
          reactome = ReactomePA::enrichPathway(fg, organism = if (sp == "mouse") "mouse" else "human",
                                               universe = uni, pvalueCutoff = 1, qvalueCutoff = 1,
                                               minGSSize = 3, readable = TRUE),
          gobp     = clusterProfiler::enrichGO(fg, OrgDb = org, ont = "BP", universe = uni,
                                               pvalueCutoff = 1, qvalueCutoff = 1, minGSSize = 3, readable = TRUE),
          gomf     = clusterProfiler::enrichGO(fg, OrgDb = org, ont = "MF", universe = uni,
                                               pvalueCutoff = 1, qvalueCutoff = 1, minGSSize = 3, readable = TRUE),
          gocc     = clusterProfiler::enrichGO(fg, OrgDb = org, ont = "CC", universe = uni,
                                               pvalueCutoff = 1, qvalueCutoff = 1, minGSSize = 3, readable = TRUE),
          clusterProfiler::enrichKEGG(fg, organism = organism, universe = uni, pvalueCutoff = 1))
        incProgress(1); r
      }, error = function(e) e)
    })
    if (inherits(ek, "error"))
      return(list(err = sprintf("Enrichment failed%s: %s",
                                if (db == "kegg") " (KEGG needs internet)" else "", conditionMessage(ek))))
    rd <- if (is.null(ek)) data.frame() else as.data.frame(ek)
    df <- if (!nrow(rd)) data.frame() else data.frame(
      Term = rd$Description, P.value = rd$pvalue, Adjusted.P.value = rd$p.adjust,
      Overlap = paste0(rd$Count, "/", sub(".*/", "", rd$GeneRatio)),
      Genes = rd$geneID, stringsAsFactors = FALSE)
    list(res = df, lib = .enrich_lib_label(db, sp), species = sp,
         n_genes = length(fg), bg = bg_used, n_uni = length(uni))
  })

  # All terms passing the user's significance cutoff (>= 3 genes), p-sorted.
  enrich_sig <- reactive({
    e <- enrich_res(); req(e)
    if (!is.null(e$err) || is.null(e$res) || !nrow(e$res)) return(NULL)
    r <- e$res
    parts  <- strsplit(input$enrich_sig %||% "adj_0.05", "_")[[1]]
    metric <- parts[1]; cut <- as.numeric(parts[2])
    pv <- if (metric == "adj") r$Adjusted.P.value else r$P.value
    k  <- as.integer(sub("/.*", "", r$Overlap))
    r  <- r[is.finite(pv) & k >= 3 & pv < cut, , drop = FALSE]
    if (!nrow(r)) return(r[0, , drop = FALSE])
    r[order(r$P.value), , drop = FALSE]
  })
  enrich_top <- reactive({
    r <- enrich_sig(); if (is.null(r) || !nrow(r)) return(r)
    head(r, as.integer(input$enrich_top %||% 15L))
  })

  output$enrich_plot <- renderPlot({
    e <- enrich_res()
    if (is.null(e)) { plot.new(); title("Pick an ontology and click 'Run enrichment'."); return() }
    if (!is.null(e$err)) { plot.new(); title(e$err); return() }
    r <- enrich_top()
    if (is.null(r) || !nrow(r)) {
      plot.new(); title(sprintf("No terms pass %s with >=3 genes.",
                                gsub("_", " ", input$enrich_sig %||% "adj 0.05"))); return()
    }
    term <- sub(" \\(GO:.*\\)$", "", r$Term)              # strip GO id suffix
    term <- ifelse(nchar(term) > 46, paste0(substr(term, 1, 43), "..."), term)
    val  <- -log10(r$P.value)
    ord  <- order(val)                                    # most significant on top
    val <- val[ord]; term <- term[ord]
    pal <- grDevices::colorRampPalette(c("#fee0b6", "#f1a340", "#b35806"))(length(val))
    # Term names as y-axis labels (las=1) with a wide left margin; no in-bar text.
    op <- par(mar = c(4.2, 15, 2.6, 1.2), mgp = c(2.4, 0.6, 0)); on.exit(par(op))
    barplot(val, horiz = TRUE, names.arg = term, las = 1, col = pal, border = NA,
            xlab = "-log10(P)", cex.names = 0.82, cex.lab = 1.15,
            main = sprintf("Top %d enriched terms (%s)", nrow(r), e$lib),
            cex.main = 1.2, xlim = c(0, max(val) * 1.02))
  },
  # fixed width (centered via CSS) so it doesn't over-stretch; height grows with
  # the number of bars so they stay a readable thickness.
  width  = function() 760L,
  height = function() {
    r <- tryCatch(enrich_top(), error = function(e) NULL)
    n <- if (is.null(r) || !nrow(r)) 1L else nrow(r)
    as.integer(120L + n * 26L)
  })

  output$enrich_note <- renderUI({
    e <- enrich_res()
    if (is.null(e)) return(helpText(em("Pick an ontology and background, set the cutoff, and click 'Run enrichment'.")))
    if (!is.null(e$err)) return(helpText(em(e$err)))
    nsig <- if (is.null(enrich_sig())) 0L else nrow(enrich_sig())
    bg <- if (identical(e$bg, "tested"))
            sprintf("tested-gene background (%s genes)", format(e$n_uni, big.mark = ","))
          else "genome-wide background"
    helpText(em(sprintf("clusterProfiler %s, %s. %d rhythmic genes submitted. %d terms significant at %s; showing top %d. x-axis = -log10(P).",
                        e$lib, bg, e$n_genes, nsig,
                        gsub("_", " ", input$enrich_sig %||% "adj 0.05"),
                        min(nsig, as.integer(input$enrich_top %||% 15L)))))
  })

  output$dl_enrich_top <- downloadHandler(
    filename = function() sprintf("SCP_enrichment_%s_top_%s.csv", input$enrich_db %||% "kegg", .gene_fname()),
    content  = function(f) {
      r <- enrich_top()
      utils::write.csv(if (is.null(r) || !nrow(r)) data.frame(note = "no significant terms") else r,
                       f, row.names = FALSE)
    })
  output$dl_enrich_all <- downloadHandler(
    filename = function() sprintf("SCP_enrichment_%s_allsig_%s.csv", input$enrich_db %||% "kegg", .gene_fname()),
    content  = function(f) {
      r <- enrich_sig()
      utils::write.csv(if (is.null(r) || !nrow(r)) data.frame(note = "no significant terms") else r,
                       f, row.names = FALSE)
    })
}

shinyApp(ui, server)
