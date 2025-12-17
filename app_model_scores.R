# Lightweight app for server-side files:  read single models table + protein abundance file -> generate scores & compare models
# - Expects R/generate_scores.R and R/model_comparison.R to be provided in R/
# - Uses data. table for fast IO and processing
#
# Expectations for generate_scores: 
#   generate_scores(models_dt, proteins_dt, selected_models, standardize = FALSE, platform = NULL)
#   Should return a data.table/data.frame (or list with $scores) containing 'individual_id' and columns 'score_<modelname>'.
#
# Key UI/behavior:
# - User provides server-side paths for models and proteins. 
# - After loading models, user sees a preview and can select which columns define grouping (model structure).
# - If no grouping columns chosen, app treats entire file as a single model.
# - Detects id-like columns in proteins; user selects which is the primary id.
# - Platform selector (3k / HT) passed to generate_scores.
# - Text input for server-side path to save scores (auto-save after generation + manual save button).
#
# Run:  shiny::runApp("app_model_scores. R", host="0.0.0.0", port=3838)

library(shiny)
library(shinydashboard)
library(DT)
library(ggplot2)
library(data.table)
library(tools)

# Source optional helpers provided by the user
if (file.exists("R/model_comparison.R")) {
  tryCatch(source("R/model_comparison. R"), error = function(e) message("Error sourcing R/model_comparison.R:  ", e$message))
}
if (file.exists("R/generate_scores.R")) {
  tryCatch(source("R/generate_scores.R"), error = function(e) message("Error sourcing R/generate_scores. R: ", e$message))
}

# ---------------------------------------------------------------------
# Fast reader using data.table
# ---------------------------------------------------------------------
read_local_table_dt <- function(path) {
  if (! file.exists(path)) stop("File not found:  ", path)
  dt <- tryCatch({
    data.table::fread(path)
  }, error = function(e) {
    stop("fread failed for ", path, ": ", e$message)
  })
  return(dt)
}

model_name_from_path <- function(path) {
  basename(tools::file_path_sans_ext(path))
}

DEFAULT_GROUP_CANDIDATES <- c("cohort", "phenotype", "method", "component", "compa", "compb", "compc", "compd")

# ---------------------------------------------------------------------
# UI
# ---------------------------------------------------------------------
ui <- dashboardPage(
  dashboardHeader(title = "Model -> Scores (server files, data. table)"),
  dashboardSidebar(width = 340,
    sidebarMenu(
      menuItem("Data", tabName = "data", icon = icon("database")),
      menuItem("Generate Scores", tabName = "scores", icon = icon("calculator")),
      menuItem("Compare Models", tabName = "compare", icon = icon("balance-scale"))
    )
  ),
  dashboardBody(
    tabItems(
      tabItem(tabName = "data",
        fluidRow(
          box(title = "Models (path on instance)", width = 8, status = "primary", solidHeader = TRUE,
              p("Provide a single models table path on the server.  The table must include columns: protein_name, coefficient.  You will select which columns define models after preview."),
              textInput("models_path", "Models file path:", placeholder = "/data/models_table.csv"),
              actionButton("load_models_btn", "Load & preview models file", class = "btn-primary")
          ),
          box(title = "Proteins (path on instance)", width = 4, status = "primary", solidHeader = TRUE,
              p("Protein abundance file (CSV or TSV) on the server. Required columns: an id-like column, protein_name, protein_level"),
              textInput("proteins_path", "Proteins file path:", placeholder = "/data/proteins_testcohort.csv"),
              actionButton("load_proteins_btn", "Load proteins"),
              hr(),
              uiOutput("proteins_id_ui")
          )
        ),
        fluidRow(
          box(title = "Loaded models (preview & grouping)", width = 12, 
              DT::dataTableOutput("models_preview"),
              hr(),
              uiOutput("models_grouping_ui")  # dynamic:  column picker + apply button
          )
        ),
        fluidRow(
          box(title = "Loaded proteins (preview)", width = 6, DT:: dataTableOutput("proteins_preview")),
          box(title = "Model filters (available after grouping)", width = 6, uiOutput("models_filters_ui"))
        ),
        fluidRow(
          box(width = 12, verbatimTextOutput("load_message"))
        )
      ),
      
      tabItem(tabName = "scores",
        fluidRow(
          box(title = "Score generation options", width = 4, status = "warning", solidHeader = TRUE,
              uiOutput("select_models_ui"),
              fluidRow(
                column(6, actionButton("select_all_models", "Select all", width = "100%")),
                column(6, actionButton("clear_models", "Clear", width = "100%"))
              ),
              checkboxInput("standardize_scores", "Standardize scores (z-score)", value = FALSE),
              selectInput("platform", "Platform:", choices = c("3k","HT"), selected = "3k"),
              br(),
              actionButton("generate_scores", "Generate scores", class = "btn-primary"),
              br(), br(),
              textInput("output_scores_path", "Output scores path (server):", placeholder = "/data/output/scores.csv"),
              actionButton("save_scores_now", "Save scores now", class = "btn-default")
          ),
          box(title = "Score distribution controls", width = 8, status = "info", solidHeader = TRUE,
              radioButtons("score_plot_mode", "Plot mode", choices = c("Density"="density","Binned mean (equal-sized bins)"="binned"), inline = TRUE),
              conditionalPanel(condition = "input.score_plot_mode == 'binned'",
                               sliderInput("score_nbins", "Number of bins", min=2, max=50, value=10)
              )
          )
        ),
        fluidRow(
          box(title = "Scores (preview)", width = 12, DT::dataTableOutput("scores_table"))
        ),
        fluidRow(
          box(title = "Score distribution", width = 8, plotOutput("scores_dist_plot", height = "420px")),
          box(title = "Per-model proteins used/missing", width = 4, DT::dataTableOutput("proteins_used_table"))
        )
      ),
      
      tabItem(tabName = "compare",
        fluidRow(
          box(title = "Model overlap & score correlations", width = 12, status = "primary", solidHeader = TRUE,
              verbatimTextOutput("overlap_summary"),
              DT::dataTableOutput("overlap_table"),
              br(),
              conditionalPanel(condition = "output.multipleModels == true",
                               h4("Score correlation matrix"),
                               plotOutput("score_corr_heatmap", height = "420px")
              )
          )
        )
      )
    )
  )
)

# ---------------------------------------------------------------------
# Server
# ---------------------------------------------------------------------
server <- function(input, output, session) {
  rv <- reactiveValues(
    models_dt_raw = NULL,    # raw models table before grouping
    models_dt = NULL,        # final models table with model_name after grouping
    proteins_dt_raw = NULL,  # original loaded proteins
    proteins_dt = NULL,      # normalized proteins with chosen individual_id
    scores_dt = NULL,
    load_msgs = character(),
    detected_groups = character(),
    detected_id_cols = character()
  )
  
  # helper functions
  detect_id_like_cols <- function(dt) {
    if (is.null(dt)) return(character(0))
    cols <- names(dt)
    cols[grepl("id", cols, ignore.case = TRUE)]
  }
  find_col_name <- function(dt, target_lc) {
    nm <- names(dt); idx <- which(tolower(nm) == target_lc)
    if (length(idx) == 1) return(nm[idx])
    return(NULL)
  }
  
  # ---------- Load models (step 1: load raw, show preview, let user define grouping) ----------
  observeEvent(input$load_models_btn, {
    req(input$models_path)
    p <- trimws(input$models_path)
    if (! nzchar(p) || !file.exists(p)) {
      rv$models_dt_raw <- NULL
      rv$models_dt <- NULL
      rv$load_msgs <- c(rv$load_msgs, paste0("Models file not found: ", p))
      return()
    }
    tryCatch({
      dt <- read_local_table_dt(p); setDT(dt)
      # validate required columns (case-insensitive)
      pn_col <- find_col_name(dt, "protein_name"); coef_col <- find_col_name(dt, "coefficient")
      if (is.null(pn_col) || is.null(coef_col)) stop("Models file must contain columns 'protein_name' and 'coefficient' (case-insensitive)")
      # normalize those columns to standard names
      if (pn_col != "protein_name") setnames(dt, pn_col, "protein_name")
      if (coef_col != "coefficient") setnames(dt, coef_col, "coefficient")
      dt[, protein_name := as.character(protein_name)]
      dt[, coefficient := as. numeric(coefficient)]
      
      # store raw loaded table (before constructing model_name)
      rv$models_dt_raw <- copy(dt)
      rv$models_dt <- NULL  # clear final models_dt until user applies grouping
      rv$detected_groups <- character()
      
      # detect candidate grouping columns (exclude protein_name, coefficient, model_name if present)
      exclude <- c("protein_name", "coefficient", "model_name")
      all_cols <- names(dt)
      candidate_groups <- setdiff(all_cols, exclude)
      # optionally filter to known candidates (you can remove this to allow all columns)
      suggested <- intersect(candidate_groups, DEFAULT_GROUP_CANDIDATES)
      
      # Render grouping UI:  multi-select + Apply button
      output$models_grouping_ui <- renderUI({
        tagList(
          p("Select which column(s) define models (e.g., cohort, phenotype, method, component). If left empty, the app will treat the entire file as a single model."),
          selectizeInput("grouping_columns", "Grouping columns:", 
                         choices = candidate_groups, 
                         selected = suggested,  # pre-select suggested candidates
                         multiple = TRUE, 
                         options = list(placeholder = "Pick columns or leave empty for single model")),
          actionButton("apply_grouping_btn", "Apply grouping & build model names", class = "btn-success")
        )
      })
      
      rv$load_msgs <- c(rv$load_msgs, paste0("Models file loaded (raw): ", basename(p), " (", nrow(dt), " rows). Review preview and select grouping columns, then click 'Apply grouping'. "))
    }, error = function(e) {
      rv$models_dt_raw <- NULL
      rv$models_dt <- NULL
      rv$load_msgs <- c(rv$load_msgs, paste0("Error loading models:  ", e$message))
    })
  })
  
  # Show preview of raw models table
  output$models_preview <- DT:: renderDataTable({ 
    req(rv$models_dt_raw); 
    DT::datatable(head(rv$models_dt_raw, 200), options = list(scrollX = TRUE, pageLength = 10)) 
  })
  
  # ---------- Apply grouping:  build model_name based on user selection ----------
  observeEvent(input$apply_grouping_btn, {
    req(rv$models_dt_raw)
    dt <- copy(rv$models_dt_raw)
    chosen_groups <- input$grouping_columns %||% character(0)
    
    # If user chose grouping columns, construct model_name by pasting them
    if (length(chosen_groups) > 0) {
      # validate chosen columns exist
      missing <- setdiff(chosen_groups, names(dt))
      if (length(missing) > 0) {
        showNotification(paste("Selected grouping columns not found in table:", paste(missing, collapse=", ")), type = "error")
        return()
      }
      dt[, model_name := do.call(paste, c(. SD, sep = "_")), .SDcols = chosen_groups]
      rv$detected_groups <- chosen_groups
      rv$load_msgs <- c(rv$load_msgs, paste0("Grouping applied: model_name built from [", paste(chosen_groups, collapse = ", "), "]. ", length(unique(dt$model_name)), " model(s) created."))
    } else {
      # No grouping:  treat as single model
      dt[, model_name := "single_model"]
      rv$detected_groups <- character(0)
      rv$load_msgs <- c(rv$load_msgs, "No grouping columns selected. Treated entire table as a single model named 'single_model'.")
    }
    
    # keep necessary columns:  model_name, protein_name, coefficient, plus any grouping columns
    keep_cols <- unique(c("model_name", "protein_name", "coefficient", chosen_groups))
    dt2 <- dt[, ..keep_cols]
    dt2[, model_name := as.character(model_name)]
    
    # store final models table
    rv$models_dt <- dt2
    showNotification(paste("Model names constructed:", length(unique(dt2$model_name)), "model(s)"), type = "message")
  })
  
  # render filters UI (only after models_dt exists after Apply grouping)
  output$models_filters_ui <- renderUI({
    req(rv$models_dt)
    present_cols <- rv$detected_groups
    if (length(present_cols) == 0) return(helpText("No grouping columns selected; filters not applicable."))
    inputs <- lapply(present_cols, function(col) {
      vals <- sort(unique(na.omit(rv$models_dt[[col]])))
      selectInput(inputId = paste0("filter_", col), label = col, choices = c("All", vals), selected = "All", multiple = FALSE, selectize = TRUE)
    })
    tagList(p("Use filters to narrow available models. "), do.call(fluidRow, lapply(inputs, function(x) column(width = floor(12/length(inputs)), x))))
  })
  
  # ---------- Load proteins and detect id-like columns ----------
  observeEvent(input$load_proteins_btn, {
    req(input$proteins_path)
    p <- trimws(input$proteins_path)
    if (!nzchar(p) || !file.exists(p)) {
      rv$proteins_dt_raw <- NULL
      rv$load_msgs <- c(rv$load_msgs, paste0("Proteins file not found: ", p))
      return()
    }
    tryCatch({
      dt <- read_local_table_dt(p); setDT(dt)
      rv$proteins_dt_raw <- dt
      rv$detected_id_cols <- detect_id_like_cols(dt)
      rv$load_msgs <- c(rv$load_msgs, paste0("Proteins loaded: ", basename(p), " (", nrow(dt), " rows). Detected id-like columns: ", ifelse(length(rv$detected_id_cols)>0, paste(rv$detected_id_cols, collapse=", "), "<none>")))
    }, error = function(e) {
      rv$proteins_dt_raw <- NULL
      rv$load_msgs <- c(rv$load_msgs, paste0("Error loading proteins: ", e$message))
    })
  })
  
  output$proteins_preview <- DT::renderDataTable({ req(rv$proteins_dt_raw); DT::datatable(head(rv$proteins_dt_raw,200), options = list(scrollX = TRUE, pageLength = 10)) })
  output$load_message <- renderText({ paste(rv$load_msgs, collapse = "\n") })
  
  # ---------- dynamic UI for id selection ----------
  output$proteins_id_ui <- renderUI({
    if (is.null(rv$proteins_dt_raw)) return(helpText("Load proteins file to detect id columns. "))
    id_cols <- rv$detected_id_cols
    if (length(id_cols) == 0) {
      tagList(helpText("No id-like columns detected.  Enter the column name to use as individual id: "), textInput("proteins_id_manual", "Individual id column name:", value = "individual_id"))
    } else {
      tagList(helpText("Pick which detected id column to use as the primary individual id:"), selectInput("proteins_id_choice", "Individual id column:", choices = id_cols, selected = id_cols[1], multiple = FALSE, selectize = TRUE))
    }
  })
  
  # Helper:  normalize proteins table (rename chosen id to 'individual_id')
  normalize_proteins_dt <- reactive({
    req(rv$proteins_dt_raw)
    dt <- copy(rv$proteins_dt_raw)
    chosen <- NULL
    if (! is.null(input$proteins_id_choice) && nzchar(input$proteins_id_choice)) chosen <- input$proteins_id_choice
    if (is.null(chosen) || ! nzchar(chosen)) {
      if (!is.null(input$proteins_id_manual) && nzchar(input$proteins_id_manual)) chosen <- input$proteins_id_manual
    }
    if (is.null(chosen) || !nzchar(chosen)) {
      if ("individual_id" %in% names(dt)) chosen <- "individual_id"
    }
    if (is.null(chosen) || !nzchar(chosen)) stop("No individual id column selected or found.")
    if (!(chosen %in% names(dt))) stop(paste0("Selected id column '", chosen, "' not found in proteins table"))
    # detect protein_name and protein_level (case-insensitive)
    if (! ("protein_name" %in% names(dt))) {
      nm <- names(dt)[tolower(names(dt)) == "protein_name"]
      if (length(nm)==1) setnames(dt, nm, "protein_name")
    }
    if (!("protein_level" %in% names(dt))) {
      nm <- names(dt)[tolower(names(dt)) == "protein_level"]
      if (length(nm)==1) setnames(dt, nm, "protein_level")
    }
    if (! all(c("protein_name","protein_level") %in% names(dt))) stop("Proteins table must contain 'protein_name' and 'protein_level' columns (case-insensitive)")
    if (chosen != "individual_id") setnames(dt, chosen, "individual_id")
    dt[, individual_id := as.character(individual_id)]
    dt[, protein_name := as.character(protein_name)]
    dt[, protein_level := as. numeric(protein_level)]
    dt
  })
  
  # ---------- filtered model names reactive ----------
  filtered_model_names <- reactive({
    req(rv$models_dt)
    dt <- copy(rv$models_dt)
    filter_inputs <- paste0("filter_", rv$detected_groups)
    if (length(filter_inputs) == 0) return(unique(dt$model_name))
    for (fn in filter_inputs) {
      if (! is.null(input[[fn]]) && nzchar(input[[fn]]) && input[[fn]] != "All") {
        colname <- sub("^filter_", "", fn)
        real_col <- names(dt)[tolower(names(dt)) == tolower(colname)][1]
        if (! is.null(real_col) && real_col != "") dt <- dt[get(real_col) %in% input[[fn]]]
      }
    }
    unique(dt$model_name)
  })
  
  # ---------- Model selection UI ----------
  output$select_models_ui <- renderUI({
    choices <- filtered_model_names()
    if (is.null(choices) || length(choices) == 0) helpText("No models available for selection.  Apply grouping first or check filters.")
    else selectizeInput("models_selected", "Models to score:", choices = sort(choices), multiple = TRUE, selected = choices[1:min(2, length(choices))], options = list(placeholder = "Search/select models", maxOptions = 500, allowEmptyOption = TRUE))
  })
  
  observeEvent(input$select_all_models, { choices <- filtered_model_names(); if (! is.null(choices) && length(choices)>0) updateSelectizeInput(session, "models_selected", selected = choices) })
  observeEvent(input$clear_models, { updateSelectizeInput(session, "models_selected", selected = character(0)) })
  
  observe({
    choices <- filtered_model_names()
    if (! is.null(choices)) {
      current_sel <- isolate(input$models_selected)
      new_sel <- current_sel[current_sel %in% choices]
      if (length(new_sel) == 0 && length(choices) > 0) new_sel <- choices[1:min(2,length(choices))]
      tryCatch(updateSelectizeInput(session, "models_selected", choices = sort(choices), selected = new_sel, server = TRUE), error = function(e) {})
    }
  })
  
  # ---------- Generate scores:  call user-supplied generate_scores (with platform) ----------
  observeEvent(input$generate_scores, {
    req(rv$models_dt, rv$proteins_dt_raw, input$models_selected)
    proteins_for_scoring <- tryCatch(normalize_proteins_dt(), error = function(e) { showNotification(paste("Proteins normalization error:", e$message), type = "error"); NULL })
    if (is.null(proteins_for_scoring)) return()
    selected <- input$models_selected; standardize <- isTRUE(input$standardize_scores); platform_choice <- input$platform %||% "3k"
    if (! exists("generate_scores", mode = "function")) { showNotification("generate_scores(... ) function not found.  Please provide R/generate_scores.R", type = "error"); return() }
    gen_res <- NULL
    tryCatch({
      gen_res <- generate_scores(models_dt = copy(rv$models_dt), proteins_dt = copy(proteins_for_scoring), selected_models = selected, standardize = standardize, platform = platform_choice)
    }, error = function(e) {
      msg <- e$message %||% ""
      if (grepl("unused argument|formal argument|unused", msg, ignore.case = TRUE)) {
        tryCatch({ gen_res <<- generate_scores(models_dt = copy(rv$models_dt), proteins_dt = copy(proteins_for_scoring), selected_models = selected, standardize = standardize) }, error = function(e2) { showNotification(paste("generate_scores retry error:", e2$message), type = "error"); gen_res <<- NULL })
      } else {
        showNotification(paste("generate_scores error:", e$message), type = "error"); gen_res <<- NULL
      }
    })
    if (is.null(gen_res)) return()
    if (is.list(gen_res) && !is.null(gen_res$scores)) scores_dt <- as.data.table(gen_res$scores) else scores_dt <- as.data.table(gen_res)
    if (! ("individual_id" %in% names(scores_dt))) { showNotification("generate_scores did not return 'individual_id' column", type = "error"); return() }
    expected_cols <- paste0("score_", selected); missing_cols <- setdiff(expected_cols, names(scores_dt))
    if (length(missing_cols) > 0) { showNotification(paste("Missing score columns from generate_scores:", paste(missing_cols, collapse = ", ")), type = "error"); return() }
    rv$scores_dt <- scores_dt
    showNotification("Scores generated and stored in memory", type = "message")
    # Auto-save if output path provided
    outp <- trimws(input$output_scores_path %||% "")
    if (nzchar(outp)) {
      tryCatch({
        dir.create(dirname(outp), recursive = TRUE, showWarnings = FALSE)
        data.table:: fwrite(rv$scores_dt, outp)
        showNotification(paste("Scores written to", outp), type = "message")
      }, error = function(e) {
        showNotification(paste("Error writing scores to path:", e$message), type = "error")
      })
    }
  })
  
  # Manual save button
  observeEvent(input$save_scores_now, {
    req(rv$scores_dt)
    outp <- trimws(input$output_scores_path %||% "")
    if (!nzchar(outp)) { showNotification("Please provide an output path in the input field first.", type = "warning"); return() }
    tryCatch({
      dir.create(dirname(outp), recursive = TRUE, showWarnings = FALSE)
      data.table::fwrite(rv$scores_dt, outp)
      showNotification(paste("Scores written to", outp), type = "message")
    }, error = function(e) {
      showNotification(paste("Error writing scores to path:", e$message), type = "error")
    })
  })
  
  output$scores_table <- DT::renderDataTable({ req(rv$scores_dt); DT::datatable(head(rv$scores_dt,200), options = list(scrollX = TRUE, pageLength = 10)) })
  
  # ---------- Proteins used / missing per model ----------
  proteins_used_table_dt <- reactive({
    req(rv$models_dt, rv$proteins_dt_raw, input$models_selected)
    sel <- input$models_selected
    pn_col <- names(rv$proteins_dt_raw)[tolower(names(rv$proteins_dt_raw)) == "protein_name"]
    present_prots <- if (length(pn_col)==1) unique(rv$proteins_dt_raw[[pn_col]]) else unique(rv$proteins_dt_raw[[which(tolower(names(rv$proteins_dt_raw))=="protein_name")]])
    out_list <- vector("list", length(sel))
    for (i in seq_along(sel)) {
      m <- sel[i]
      reqs <- unique(rv$models_dt[model_name == m, protein_name])
      pres <- intersect(reqs, present_prots)
      miss <- setdiff(reqs, pres)
      out_list[[i]] <- data.table(model = m, n_required = length(reqs), n_present = length(pres), n_missing = length(miss), missing_proteins = if (length(miss)>0) paste(miss, collapse = "; ") else "")
    }
    rbindlist(out_list)
  })
  output$proteins_used_table <- DT::renderDataTable({ req(proteins_used_table_dt()); DT::datatable(proteins_used_table_dt(), options = list(scrollX = TRUE, pageLength = 10)) })
  
  # ---------- Score distribution ----------
  output$scores_dist_plot <- renderPlot({
    req(rv$scores_dt, input$models_selected)
    dt <- copy(rv$scores_dt); sel <- input$models_selected; score_cols <- paste0("score_", sel); available <- intersect(score_cols, names(dt))
    if (length(available) == 0) { plot. new(); text(0.5,0.5,"No scores available.  Generate scores first."); return() }
    long_dt <- melt(dt, id.vars = "individual_id", measure.vars = available, variable.name = "score_name", value.name = "score", variable.factor = FALSE)
    long_dt[, score_name := sub("^score_", "", score_name)]
    if (input$score_plot_mode == "density") {
      ggplot(long_dt, aes(x = score, color = score_name, fill = score_name)) + geom_density(alpha = 0.25, na.rm = TRUE) + theme_minimal() + xlab("Score") + ylab("Density")
    } else {
      nb <- as.integer(input$score_nbins %||% 10)
      long_dt[, bin := as.integer(ntile(score, nb)), by = .(score_name)]
      summarize_dt <- long_dt[! is.na(score), .(n = .N, mean = mean(score, na.rm = TRUE), sd = sd(score, na.rm = TRUE), se = ifelse(. N>0, sd(score, na. rm = TRUE)/sqrt(.N), NA_real_)), by = .(score_name, bin)]
      summarize_dt[, ci_low := mean - 1.96 * se]; summarize_dt[, ci_high := mean + 1.96 * se]
      ggplot(summarize_dt, aes(x = factor(bin), y = mean, group = 1)) + geom_point() + geom_errorbar(aes(ymin = ci_low, ymax = ci_high), width = 0.2) + facet_wrap(~ score_name, scales = "free_y") + xlab("Score bin") + ylab("Mean score ± 95% CI") + theme_minimal()
    }
  })
  
  # ---------- Overlap & correlations ----------
  overlap_info_dt <- reactive({
    req(rv$models_dt, input$models_selected)
    sel <- input$models_selected
    prot_sets <- lapply(sel, function(m) unique(rv$models_dt[model_name == m, protein_name]))
    names(prot_sets) <- sel
    all_prots <- unique(unlist(prot_sets))
    presence_mat <- sapply(prot_sets, function(s) all_prots %in% s)
    presence_dt <- as.data.table(cbind(protein_name = all_prots, as.data.table(presence_mat)))
    overlap_count <- length(Reduce(intersect, prot_sets))
    unique_counts <- sapply(sel, function(m) length(setdiff(prot_sets[[m]], unlist(prot_sets[names(prot_sets) != m]))))
    list(presence_dt = presence_dt, overlap_count = overlap_count, unique_counts = unique_counts, prot_sets = prot_sets)
  })
  output$overlap_summary <- renderText({ req(overlap_info_dt()); oi <- overlap_info_dt(); paste0("Overlapping proteins across selected models: ", oi$overlap_count, "\nUnique proteins per model:\n", paste(names(oi$unique_counts), oi$unique_counts, sep=":  ", collapse = "; ")) })
  output$overlap_table <- DT:: renderDataTable({ req(overlap_info_dt()); DT::datatable(overlap_info_dt()$presence_dt, options = list(scrollX = TRUE, pageLength = 10)) })
  outputOptions(output, "multipleModels", suspendWhenHidden = FALSE)
  output$multipleModels <- reactive({ length(input$models_selected %||% character(0)) >= 2 })
  output$score_corr_heatmap <- renderPlot({
    req(rv$scores_dt, input$models_selected)
    sel <- input$models_selected
    score_cols <- paste0("score_", sel)
    present <- intersect(score_cols, names(rv$scores_dt))
    if (length(present) < 2) { plot.new(); text(0.5,0.5,"At least two generated scores required for correlation plot. "); return() }
    mat_dt <- rv$scores_dt[, ..present]; mat <- as.matrix(mat_dt); colnames(mat) <- sub("^score_", "", present)
    corm <- cor(mat, use = "pairwise.complete.obs"); corm_melt <- as.data.table(as.table(corm)); setnames(corm_melt, c("ModelX","ModelY","Correlation"))
    ggplot(corm_melt, aes(x = ModelX, y = ModelY, fill = Correlation)) + geom_tile() + scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) + geom_text(aes(label = round(Correlation, 2))) + theme_minimal() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
  })
}

shinyApp(ui, server)
