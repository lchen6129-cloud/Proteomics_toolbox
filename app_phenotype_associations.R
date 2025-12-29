# Phenotype Association Analysis App
# 
# Takes score file (from scoring app) + phenotype file -> performs association analysis
# 
# Expected inputs:
# - Score file: rgcpid, score_model1, score_model2, ...   
# - Phenotype file: rgcpid (or RGCpID), phenotype columns, covariates
#
# Expected functions in R/phenotype_association.R:  
# - plot_continuous_distribution(pheno_vector, pheno_name)
# - get_continuous_summary_stats(pheno_vector)
# - rint_transform(pheno_vector)  # NEW: RINTing function
# - run_linear_regression(outcome, predictors, covariates, data, stratify_by = NULL)
# - plot_linear_regression_results(results, n_bins = 10)  # UPDATED: added n_bins
# - plot_binary_distribution(pheno_vector, pheno_name)
# - get_binary_summary_stats(pheno_vector)
# - run_logistic_regression(outcome, predictors, covariates, data, stratify_by = NULL)
# - plot_logistic_regression_results(results, n_bins = 10)  # UPDATED: added n_bins
# - plot_survival_curves(time, event, pheno_name)
# - get_survival_summary_stats(time, event)
# - run_cox_regression(time, event, predictors, covariates, data, stratify_by = NULL)
# - plot_cox_regression_results(results, n_bins = 10, high_risk_cutoff = 0.5)  # UPDATED: added n_bins and high_risk_cutoff
# - get_high_risk_summary(data, score_col, cutoff_percentile)  # NEW: for high-risk group analysis
#
# Run:    shiny::runApp("app_phenotype_association.R", host="0.0.0.0", port=3838)

library(shiny)
library(shinydashboard)
library(DT)
library(ggplot2)
library(data.table)
library(tools)

# Source association functions
if (file.exists("R/phenotype_association.R")) {
  tryCatch(source("R/phenotype_association. R"), error = function(e) message("Error sourcing R/phenotype_association.R:    ", e$message))
} else {
  message("R/phenotype_association.R not found - you must provide this file with analysis functions")
}

# ---------------------------------------------------------------------
# Helper functions
# ---------------------------------------------------------------------
read_local_table_dt <- function(path) {
  if (!file.exists(path)) stop("File not found: ", path)
  dt <- tryCatch({
    data.table::fread(path, nThread = 4, showProgress = TRUE)
  }, error = function(e) {
    stop("fread failed for ", path, ": ", e$message)
  })
  return(dt)
}

detect_id_like_cols <- function(dt) {
  if (is.null(dt)) return(character(0))
  cols <- names(dt)
  id_like <- cols[grepl("rgcpid|id", cols, ignore.case = TRUE)]
  return(id_like)
}

detect_score_cols <- function(dt) {
  if (is.null(dt)) return(character(0))
  cols <- names(dt)
  score_cols <- cols[grepl("^score_", cols, ignore.case = TRUE)]
  return(score_cols)
}

# ---------------------------------------------------------------------
# UI
# ---------------------------------------------------------------------
ui <- dashboardPage(
  dashboardHeader(title = "Phenotype Association Analysis"),
  dashboardSidebar(width = 280,
    sidebarMenu(
      menuItem("Data", tabName = "data", icon = icon("database")),
      menuItem("Continuous Phenotypes", tabName = "continuous", icon = icon("chart-line")),
      menuItem("Binary Phenotypes", tabName = "binary", icon = icon("circle-half-stroke")),
      menuItem("Time-to-Event", tabName = "survival", icon = icon("clock"))
    )
  ),
  dashboardBody(
    tabItems(
      # ==================== DATA TAB ====================
      tabItem(tabName = "data",
        fluidRow(
          box(title = "Scores (from scoring app)", width = 6, status = "primary", solidHeader = TRUE,
              p("Provide path to score file.    Expected columns: rgcpid, score_model1, score_model2, ..."),
              textInput("scores_path", "Scores file path:", placeholder = "/data/scores_output.csv"),
              actionButton("load_scores_btn", "Load scores", class = "btn-primary"),
              hr(),
              uiOutput("scores_id_ui")
          ),
          box(title = "Phenotypes", width = 6, status = "primary", solidHeader = TRUE,
              p("Provide path to phenotype file.  Expected:   rgcpid/RGCpID column + phenotype/covariate columns"),
              textInput("phenotypes_path", "Phenotypes file path:", placeholder = "/data/phenotypes. csv"),
              actionButton("load_phenotypes_btn", "Load phenotypes", class = "btn-primary"),
              hr(),
              uiOutput("phenotypes_id_ui")
          )
        ),
        fluidRow(
          box(title = "Loaded scores (preview)", width = 6, 
              DT::dataTableOutput("scores_preview"),
              verbatimTextOutput("scores_summary")
          ),
          box(title = "Loaded phenotypes (preview)", width = 6, 
              DT::dataTableOutput("phenotypes_preview"),
              verbatimTextOutput("phenotypes_summary")
          )
        ),
        fluidRow(
          box(width = 12, title = "Merged data status", status = "info",
              verbatimTextOutput("merge_status"),
              actionButton("merge_data_btn", "Merge scores + phenotypes", class = "btn-success")
          )
        ),
        fluidRow(
          box(width = 12, verbatimTextOutput("load_messages"))
        )
      ),
      
      # ==================== CONTINUOUS PHENOTYPES TAB ====================
      tabItem(tabName = "continuous",
        fluidRow(
          box(title = "Phenotype Selection", width = 3, status = "primary", solidHeader = TRUE,
              uiOutput("continuous_pheno_ui"),
              checkboxInput("apply_rint_cont", "Apply RINT transformation", value = FALSE),
              actionButton("preview_continuous_btn", "Preview Distribution", class = "btn-info")
          ),
          box(title = "Distribution & Summary", width = 9, status = "info", solidHeader = TRUE,
              plotOutput("continuous_dist_plot", height = "300px"),
              verbatimTextOutput("continuous_summary_stats")
          )
        ),
        fluidRow(
          box(title = "Analysis Setup", width = 12, status = "warning", solidHeader = TRUE,
              fluidRow(
                column(2, uiOutput("scores_select_cont_ui")),
                column(2, uiOutput("covariates_select_cont_ui")),
                column(2, uiOutput("stratify_select_cont_ui")),
                column(2, 
                       br(),
                       sliderInput("n_bins_cont", "Plot bins:", min = 5, max = 50, value = 10, step = 5)
                ),
                column(4, 
                       br(),
                       actionButton("run_continuous_btn", "Run Linear Regression", class = "btn-primary", style = "width:  100%;"),
                       br(), br(),
                       textInput("continuous_output_path", "Output path (optional):", placeholder = "/data/output/continuous_results.csv"),
                       actionButton("save_continuous_btn", "Save results", class = "btn-default", style = "width: 100%;")
                )
              )
          )
        ),
        fluidRow(
          box(title = "Linear Regression Results", width = 12, 
              DT::dataTableOutput("continuous_results_table")
          )
        ),
        fluidRow(
          box(title = "Regression Plots", width = 12,
              plotOutput("continuous_regression_plot", height = "500px")
          )
        )
      ),
      
      # ==================== BINARY PHENOTYPES TAB ====================
      tabItem(tabName = "binary",
        fluidRow(
          box(title = "Phenotype Selection", width = 3, status = "primary", solidHeader = TRUE,
              uiOutput("binary_pheno_ui"),
              actionButton("preview_binary_btn", "Preview Distribution", class = "btn-info")
          ),
          box(title = "Distribution & Summary", width = 9, status = "info", solidHeader = TRUE,
              plotOutput("binary_dist_plot", height = "300px"),
              verbatimTextOutput("binary_summary_stats")
          )
        ),
        fluidRow(
          box(title = "Analysis Setup", width = 12, status = "warning", solidHeader = TRUE,
              fluidRow(
                column(2, uiOutput("scores_select_binary_ui")),
                column(2, uiOutput("covariates_select_binary_ui")),
                column(2, uiOutput("stratify_select_binary_ui")),
                column(2,
                       br(),
                       sliderInput("n_bins_binary", "Plot bins:", min = 5, max = 50, value = 10, step = 5)
                ),
                column(4,
                       br(),
                       actionButton("run_binary_btn", "Run Logistic Regression", class = "btn-primary", style = "width: 100%;"),
                       br(), br(),
                       textInput("binary_output_path", "Output path (optional):", placeholder = "/data/output/binary_results.csv"),
                       actionButton("save_binary_btn", "Save results", class = "btn-default", style = "width: 100%;")
                )
              )
          )
        ),
        fluidRow(
          box(title = "Logistic Regression Results", width = 12, 
              DT::dataTableOutput("binary_results_table")
          )
        ),
        fluidRow(
          box(title = "Regression Plots (OR/CI)", width = 12,
              plotOutput("binary_regression_plot", height = "500px")
          )
        )
      ),
      
      # ==================== TIME-TO-EVENT TAB ====================
      tabItem(tabName = "survival",
        fluidRow(
          box(title = "Phenotype Selection", width = 3, status = "primary", solidHeader = TRUE,
              uiOutput("survival_event_ui"),
              uiOutput("survival_time_ui"),
              sliderInput("censor_cutoff_years", "Censor cutoff (years):", 
                          min = 1, max = 30, value = 10, step = 1),
              actionButton("preview_survival_btn", "Preview Survival Curves", class = "btn-info")
          ),
          box(title = "Survival Curves & Summary", width = 9, status = "info", solidHeader = TRUE,
              plotOutput("survival_curves_plot", height = "300px"),
              verbatimTextOutput("survival_summary_stats")
          )
        ),
        fluidRow(
          box(title = "Analysis Setup", width = 12, status = "warning", solidHeader = TRUE,
              fluidRow(
                column(2, uiOutput("scores_select_survival_ui")),
                column(2, uiOutput("covariates_select_survival_ui")),
                column(2, uiOutput("stratify_select_survival_ui")),
                column(2,
                       br(),
                       sliderInput("n_bins_survival", "Plot bins:", min = 5, max = 50, value = 10, step = 5)
                ),
                column(4,
                       br(),
                       actionButton("run_survival_btn", "Run Cox Regression", class = "btn-primary", style = "width: 100%;"),
                       br(), br(),
                       textInput("survival_output_path", "Output path (optional):", placeholder = "/data/output/survival_results. csv"),
                       actionButton("save_survival_btn", "Save results", class = "btn-default", style = "width: 100%;")
                )
              )
          )
        ),
        fluidRow(
          box(title = "Cox Regression Results", width = 12, 
              DT::dataTableOutput("survival_results_table")
          )
        ),
        fluidRow(
          box(title = "Hazard Ratio Plots", width = 12,
              plotOutput("survival_regression_plot", height = "500px")
          )
        ),
        fluidRow(
          box(title = "High-Risk Group Analysis", width = 12, status = "success", solidHeader = TRUE,
              sliderInput("high_risk_cutoff", "High-risk cutoff (top percentile):", 
                          min = 0, max = 100, value = 50, step = 5, post = "%"),
              helpText("Define high-risk group as individuals above this percentile of the score distribution."),
              verbatimTextOutput("high_risk_summary"),
              plotOutput("high_risk_plot", height = "400px")
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
    scores_dt = NULL,
    phenotypes_dt = NULL,
    merged_dt = NULL,
    scores_id_col = NULL,
    phenotypes_id_col = NULL,
    detected_score_cols = character(),
    continuous_results = NULL,
    binary_results = NULL,
    survival_results = NULL,
    survival_analysis_data = NULL,  # Store for high-risk analysis
    load_msgs = character()
  )
  
  # ==================== LOAD SCORES ====================
  observeEvent(input$load_scores_btn, {
    req(input$scores_path)
    p <- trimws(input$scores_path)
    if (! nzchar(p) || !file.exists(p)) {
      rv$scores_dt <- NULL
      rv$load_msgs <- c(rv$load_msgs, paste0("Scores file not found: ", p))
      return()
    }
    
    showNotification("Loading scores...", id = "load_scores", duration = NULL)
    
    tryCatch({
      dt <- read_local_table_dt(p)
      setDT(dt)
      rv$scores_dt <- dt
      rv$detected_score_cols <- detect_score_cols(dt)
      rv$load_msgs <- c(rv$load_msgs, paste0("Scores loaded:    ", basename(p), " (", nrow(dt), " rows, ", length(rv$detected_score_cols), " score columns detected)"))
      removeNotification("load_scores")
      showNotification(paste("Loaded", nrow(dt), "individuals with", length(rv$detected_score_cols), "scores"), type = "message", duration = 3)
    }, error = function(e) {
      rv$scores_dt <- NULL
      rv$load_msgs <- c(rv$load_msgs, paste0("Error loading scores: ", e$message))
      removeNotification("load_scores")
      showNotification(paste("Error:", e$message), type = "error")
    })
  })
  
  output$scores_id_ui <- renderUI({
    if (is. null(rv$scores_dt)) return(helpText("Load scores file first.  "))
    id_cols <- detect_id_like_cols(rv$scores_dt)
    if (length(id_cols) == 0) {
      tagList(helpText("No ID column detected.   Enter manually:  "), 
              textInput("scores_id_manual", "ID column name:", value = "rgcpid"))
    } else {
      tagList(helpText("Select ID column: "), 
              selectInput("scores_id_choice", "ID column:", choices = id_cols, selected = id_cols[1]))
    }
  })
  
  output$scores_preview <- DT::renderDataTable({
    req(rv$scores_dt)
    dt_preview <- head(rv$scores_dt, 100)
    DT::datatable(dt_preview, options = list(scrollX = TRUE, pageLength = 10),
                  caption = paste0("Showing first 100 rows (total:  ", nrow(rv$scores_dt), " rows)"))
  })
  
  output$scores_summary <- renderText({
    req(rv$scores_dt)
    paste0("Rows: ", nrow(rv$scores_dt), "\n",
           "Columns: ", ncol(rv$scores_dt), "\n",
           "Score columns detected: ", length(rv$detected_score_cols), "\n",
           "Score names: ", paste(sub("^score_", "", rv$detected_score_cols), collapse = ", "))
  })
  
  # ==================== LOAD PHENOTYPES ====================
  observeEvent(input$load_phenotypes_btn, {
    req(input$phenotypes_path)
    p <- trimws(input$phenotypes_path)
    if (!nzchar(p) || !file.exists(p)) {
      rv$phenotypes_dt <- NULL
      rv$load_msgs <- c(rv$load_msgs, paste0("Phenotypes file not found: ", p))
      return()
    }
    
    showNotification("Loading phenotypes...", id = "load_pheno", duration = NULL)
    
    tryCatch({
      dt <- read_local_table_dt(p)
      setDT(dt)
      rv$phenotypes_dt <- dt
      rv$load_msgs <- c(rv$load_msgs, paste0("Phenotypes loaded:  ", basename(p), " (", nrow(dt), " rows, ", ncol(dt), " columns)"))
      removeNotification("load_pheno")
      showNotification(paste("Loaded", nrow(dt), "individuals with", ncol(dt), "phenotype columns"), type = "message", duration = 3)
    }, error = function(e) {
      rv$phenotypes_dt <- NULL
      rv$load_msgs <- c(rv$load_msgs, paste0("Error loading phenotypes:    ", e$message))
      removeNotification("load_pheno")
      showNotification(paste("Error:", e$message), type = "error")
    })
  })
  
  output$phenotypes_id_ui <- renderUI({
    if (is. null(rv$phenotypes_dt)) return(helpText("Load phenotypes file first. "))
    id_cols <- detect_id_like_cols(rv$phenotypes_dt)
    if (length(id_cols) == 0) {
      tagList(helpText("No ID column detected.  Enter manually:"), 
              textInput("phenotypes_id_manual", "ID column name:", value = "rgcpid"))
    } else {
      tagList(helpText("Select ID column:"), 
              selectInput("phenotypes_id_choice", "ID column:", choices = id_cols, selected = id_cols[1]))
    }
  })
  
  output$phenotypes_preview <- DT::renderDataTable({
    req(rv$phenotypes_dt)
    preview_cols <- min(20, ncol(rv$phenotypes_dt))
    dt_preview <- head(rv$phenotypes_dt[, 1:preview_cols], 100)
    DT::datatable(dt_preview, options = list(scrollX = TRUE, pageLength = 10),
                  caption = paste0("Showing first 100 rows × first ", preview_cols, " columns (total: ", 
                                   nrow(rv$phenotypes_dt), " × ", ncol(rv$phenotypes_dt), ")"))
  })
  
  output$phenotypes_summary <- renderText({
    req(rv$phenotypes_dt)
    paste0("Rows: ", nrow(rv$phenotypes_dt), "\n",
           "Columns: ", ncol(rv$phenotypes_dt))
  })
  
  output$load_messages <- renderText({
    paste(rv$load_msgs, collapse = "\n")
  })
  
  # ==================== MERGE DATA ====================
  observeEvent(input$merge_data_btn, {
    req(rv$scores_dt, rv$phenotypes_dt)
    
    scores_id <- input$scores_id_choice %||% input$scores_id_manual %||% "rgcpid"
    pheno_id <- input$phenotypes_id_choice %||% input$phenotypes_id_manual %||% "rgcpid"
    
    if (!(scores_id %in% names(rv$scores_dt))) {
      showNotification(paste("Scores ID column not found:", scores_id), type = "error")
      return()
    }
    if (!(pheno_id %in% names(rv$phenotypes_dt))) {
      showNotification(paste("Phenotypes ID column not found:", pheno_id), type = "error")
      return()
    }
    
    scores_copy <- copy(rv$scores_dt)
    pheno_copy <- copy(rv$phenotypes_dt)
    
    setnames(scores_copy, scores_id, "id")
    setnames(pheno_copy, pheno_id, "id")
    
    merged <- merge(scores_copy, pheno_copy, by = "id", all = FALSE, suffixes = c("_score", "_pheno"))
    
    rv$merged_dt <- merged
    rv$scores_id_col <- scores_id
    rv$phenotypes_id_col <- pheno_id
    
    showNotification(paste("Merged:", nrow(merged), "individuals with complete data"), type = "message")
    rv$load_msgs <- c(rv$load_msgs, paste0("Merged:    ", nrow(merged), " individuals (", 
                                            nrow(rv$scores_dt), " in scores, ", 
                                            nrow(rv$phenotypes_dt), " in phenotypes)"))
  })
  
  output$merge_status <- renderText({
    if (is.null(rv$merged_dt)) {
      "No merged data yet.    Load scores and phenotypes, then click 'Merge'."
    } else {
      paste0("✓ Merged data ready: ", nrow(rv$merged_dt), " individuals\n",
             "Columns: ", ncol(rv$merged_dt), "\n",
             "Score columns: ", length(rv$detected_score_cols), "\n",
             "Available for analysis")
    }
  })
  
  # ==================== HELPER:   Available phenotype columns ====================
  available_pheno_cols <- reactive({
    req(rv$merged_dt)
    exclude <- c("id", rv$detected_score_cols)
    setdiff(names(rv$merged_dt), exclude)
  })
  
  # ==================== CONTINUOUS PHENOTYPES ====================
  
  output$continuous_pheno_ui <- renderUI({
    req(rv$merged_dt)
    choices <- available_pheno_cols()
    selectInput("continuous_pheno", "Select continuous phenotype:", choices = choices)
  })
  
  output$scores_select_cont_ui <- renderUI({
    req(rv$merged_dt)
    selectInput("scores_cont", "Select score(s):", choices = rv$detected_score_cols, multiple = TRUE, selected = rv$detected_score_cols[1])
  })
  
  output$covariates_select_cont_ui <- renderUI({
    req(rv$merged_dt)
    choices <- available_pheno_cols()
    selectInput("covariates_cont", "Select covariates (optional):", choices = choices, multiple = TRUE)
  })
  
  output$stratify_select_cont_ui <- renderUI({
    req(rv$merged_dt)
    choices <- c("None", available_pheno_cols())
    selectInput("stratify_cont", "Stratify by (optional):", choices = choices, selected = "None")
  })
  
  # Preview distribution button
  observeEvent(input$preview_continuous_btn, {
    req(rv$merged_dt, input$continuous_pheno)
    
    pheno_name <- input$continuous_pheno
    pheno_vec <- rv$merged_dt[[pheno_name]]
    
    # Apply RINT if requested
    if (input$apply_rint_cont && exists("rint_transform")) {
      pheno_vec <- rint_transform(pheno_vec)
      pheno_name <- paste0(pheno_name, " (RINTed)")
    }
    
    output$continuous_dist_plot <- renderPlot({
      req(exists("plot_continuous_distribution"))
      plot_continuous_distribution(pheno_vec, pheno_name)
    })
    
    output$continuous_summary_stats <- renderText({
      req(exists("get_continuous_summary_stats"))
      stats <- get_continuous_summary_stats(pheno_vec)
      paste(capture.output(print(stats)), collapse = "\n")
    })
  })
  
  # Run regression
  observeEvent(input$run_continuous_btn, {
    req(rv$merged_dt, input$continuous_pheno, input$scores_cont)
    
    pheno_name <- input$continuous_pheno
    scores <- input$scores_cont
    covariates <- input$covariates_cont %||% character(0)
    stratify_var <- if (input$stratify_cont == "None") NULL else input$stratify_cont
    n_bins <- input$n_bins_cont %||% 10
    
    # Prepare data with optional RINTing
    analysis_data <- copy(rv$merged_dt)
    outcome_col <- pheno_name
    
    if (input$apply_rint_cont && exists("rint_transform")) {
      analysis_data[[paste0(pheno_name, "_rint")]] <- rint_transform(analysis_data[[pheno_name]])
      outcome_col <- paste0(pheno_name, "_rint")
      showNotification("Applied RINT transformation to outcome", type = "message", duration = 3)
    }
    
    if (exists("run_linear_regression")) {
      results <- run_linear_regression(
        outcome = outcome_col,
        predictors = scores,
        covariates = covariates,
        data = analysis_data,
        stratify_by = stratify_var
      )
      rv$continuous_results <- results
      
      output$continuous_results_table <- DT::renderDataTable({
        if (is.data.frame(results) || is.data.table(results)) {
          DT::datatable(results, options = list(scrollX = TRUE, pageLength = 25))
        } else if (is.list(results) && ! is.null(results$summary)) {
          DT::datatable(results$summary, options = list(scrollX = TRUE, pageLength = 25))
        } else {
          DT::datatable(data.frame(Message = "Results format not recognized"))
        }
      })
      
      output$continuous_regression_plot <- renderPlot({
        req(exists("plot_linear_regression_results"))
        plot_linear_regression_results(results, n_bins = n_bins)
      })
      
      showNotification("Continuous analysis complete", type = "message")
    } else {
      showNotification("run_linear_regression() not found in R/phenotype_association.R", type = "error")
    }
  })
  
  observeEvent(input$save_continuous_btn, {
    req(rv$continuous_results)
    outp <- trimws(input$continuous_output_path %||% "")
    if (!nzchar(outp)) {
      showNotification("Provide output path first.", type = "warning")
      return()
    }
    tryCatch({
      dir.create(dirname(outp), recursive = TRUE, showWarnings = FALSE)
      if (is.data.frame(rv$continuous_results) || is.data.table(rv$continuous_results)) {
        fwrite(rv$continuous_results, outp)
      } else if (is.list(rv$continuous_results) && !is.null(rv$continuous_results$summary)) {
        fwrite(rv$continuous_results$summary, outp)
      }
      showNotification(paste("Saved:", outp), type = "message")
    }, error = function(e) {
      showNotification(paste("Save error:", e$message), type = "error")
    })
  })
  
  # ==================== BINARY PHENOTYPES ====================
  
  output$binary_pheno_ui <- renderUI({
    req(rv$merged_dt)
    choices <- available_pheno_cols()
    selectInput("binary_pheno", "Select binary phenotype:", choices = choices)
  })
  
  output$scores_select_binary_ui <- renderUI({
    req(rv$merged_dt)
    selectInput("scores_binary", "Select score(s):", choices = rv$detected_score_cols, multiple = TRUE, selected = rv$detected_score_cols[1])
  })
  
  output$covariates_select_binary_ui <- renderUI({
    req(rv$merged_dt)
    choices <- available_pheno_cols()
    selectInput("covariates_binary", "Select covariates (optional):", choices = choices, multiple = TRUE)
  })
  
  output$stratify_select_binary_ui <- renderUI({
    req(rv$merged_dt)
    choices <- c("None", available_pheno_cols())
    selectInput("stratify_binary", "Stratify by (optional):", choices = choices, selected = "None")
  })
  
  observeEvent(input$preview_binary_btn, {
    req(rv$merged_dt, input$binary_pheno)
    
    pheno_name <- input$binary_pheno
    pheno_vec <- rv$merged_dt[[pheno_name]]
    
    output$binary_dist_plot <- renderPlot({
      req(exists("plot_binary_distribution"))
      plot_binary_distribution(pheno_vec, pheno_name)
    })
    
    output$binary_summary_stats <- renderText({
      req(exists("get_binary_summary_stats"))
      stats <- get_binary_summary_stats(pheno_vec)
      paste(capture. output(print(stats)), collapse = "\n")
    })
  })
  
  observeEvent(input$run_binary_btn, {
    req(rv$merged_dt, input$binary_pheno, input$scores_binary)
    
    pheno_name <- input$binary_pheno
    scores <- input$scores_binary
    covariates <- input$covariates_binary %||% character(0)
    stratify_var <- if (input$stratify_binary == "None") NULL else input$stratify_binary
    n_bins <- input$n_bins_binary %||% 10
    
    if (exists("run_logistic_regression")) {
      results <- run_logistic_regression(
        outcome = pheno_name,
        predictors = scores,
        covariates = covariates,
        data = rv$merged_dt,
        stratify_by = stratify_var
      )
      rv$binary_results <- results
      
      output$binary_results_table <- DT::renderDataTable({
        if (is. data.frame(results) || is.data.table(results)) {
          DT::datatable(results, options = list(scrollX = TRUE, pageLength = 25))
        } else if (is.list(results) && !is.null(results$summary)) {
          DT::datatable(results$summary, options = list(scrollX = TRUE, pageLength = 25))
        } else {
          DT::datatable(data.frame(Message = "Results format not recognized"))
        }
      })
      
      output$binary_regression_plot <- renderPlot({
        req(exists("plot_logistic_regression_results"))
        plot_logistic_regression_results(results, n_bins = n_bins)
      })
      
      showNotification("Binary analysis complete", type = "message")
    } else {
      showNotification("run_logistic_regression() not found in R/phenotype_association.R", type = "error")
    }
  })
  
  observeEvent(input$save_binary_btn, {
    req(rv$binary_results)
    outp <- trimws(input$binary_output_path %||% "")
    if (!nzchar(outp)) {
      showNotification("Provide output path first.", type = "warning")
      return()
    }
    tryCatch({
      dir.create(dirname(outp), recursive = TRUE, showWarnings = FALSE)
      if (is.data.frame(rv$binary_results) || is.data.table(rv$binary_results)) {
        fwrite(rv$binary_results, outp)
      } else if (is.list(rv$binary_results) && !is.null(rv$binary_results$summary)) {
        fwrite(rv$binary_results$summary, outp)
      }
      showNotification(paste("Saved:", outp), type = "message")
    }, error = function(e) {
      showNotification(paste("Save error:", e$message), type = "error")
    })
  })
  
  # ==================== TIME-TO-EVENT (SURVIVAL) ====================
  
  output$survival_event_ui <- renderUI({
    req(rv$merged_dt)
    choices <- available_pheno_cols()
    selectInput("survival_event", "Select event indicator (0/1):", choices = choices)
  })
  
  output$survival_time_ui <- renderUI({
    req(rv$merged_dt)
    choices <- available_pheno_cols()
    selectInput("survival_time", "Select time variable:", choices = choices)
  })
  
  output$scores_select_survival_ui <- renderUI({
    req(rv$merged_dt)
    selectInput("scores_survival", "Select score(s):", choices = rv$detected_score_cols, multiple = TRUE, selected = rv$detected_score_cols[1])
  })
  
  output$covariates_select_survival_ui <- renderUI({
    req(rv$merged_dt)
    choices <- available_pheno_cols()
    selectInput("covariates_survival", "Select covariates (optional):", choices = choices, multiple = TRUE)
  })
  
  output$stratify_select_survival_ui <- renderUI({
    req(rv$merged_dt)
    choices <- c("None", available_pheno_cols())
    selectInput("stratify_survival", "Stratify by (optional):", choices = choices, selected = "None")
  })
  
  observeEvent(input$preview_survival_btn, {
    req(rv$merged_dt, input$survival_event, input$survival_time)
    
    event_name <- input$survival_event
    time_name <- input$survival_time
    censor_cutoff <- input$censor_cutoff_years
    
    time_vec <- rv$merged_dt[[time_name]]
    event_vec <- rv$merged_dt[[event_name]]
    
    time_censored <- pmin(time_vec, censor_cutoff)
    event_censored <- ifelse(time_vec > censor_cutoff, 0, event_vec)
    
    output$survival_curves_plot <- renderPlot({
      req(exists("plot_survival_curves"))
      plot_survival_curves(time_censored, event_censored, event_name)
    })
    
    output$survival_summary_stats <- renderText({
      req(exists("get_survival_summary_stats"))
      stats <- get_survival_summary_stats(time_censored, event_censored)
      paste(capture.output(print(stats)), collapse = "\n")
    })
  })
  
  observeEvent(input$run_survival_btn, {
    req(rv$merged_dt, input$survival_event, input$survival_time, input$scores_survival)
    
    event_name <- input$survival_event
    time_name <- input$survival_time
    scores <- input$scores_survival
    covariates <- input$covariates_survival %||% character(0)
    stratify_var <- if (input$stratify_survival == "None") NULL else input$stratify_survival
    censor_cutoff <- input$censor_cutoff_years
    n_bins <- input$n_bins_survival %||% 10
    
    time_vec <- rv$merged_dt[[time_name]]
    event_vec <- rv$merged_dt[[event_name]]
    
    time_censored <- pmin(time_vec, censor_cutoff)
    event_censored <- ifelse(time_vec > censor_cutoff, 0, event_vec)
    
    if (exists("run_cox_regression")) {
      analysis_data <- copy(rv$merged_dt)
      analysis_data[, time_censored := time_censored]
      analysis_data[, event_censored := event_censored]
      
      # Store for high-risk analysis
      rv$survival_analysis_data <- analysis_data
      
      results <- run_cox_regression(
        time = "time_censored",
        event = "event_censored",
        predictors = scores,
        covariates = covariates,
        data = analysis_data,
        stratify_by = stratify_var
      )
      rv$survival_results <- results
      
      output$survival_results_table <- DT::renderDataTable({
        if (is.data. frame(results) || is.data.table(results)) {
          DT::datatable(results, options = list(scrollX = TRUE, pageLength = 25))
        } else if (is.list(results) && !is.null(results$summary)) {
          DT::datatable(results$summary, options = list(scrollX = TRUE, pageLength = 25))
        } else {
          DT::datatable(data. frame(Message = "Results format not recognized"))
        }
      })
      
      output$survival_regression_plot <- renderPlot({
        req(exists("plot_cox_regression_results"))
        plot_cox_regression_results(results, n_bins = n_bins)
      })
      
      showNotification("Survival analysis complete", type = "message")
    } else {
      showNotification("run_cox_regression() not found in R/phenotype_association.R", type = "error")
    }
  })
  
  # ==================== HIGH-RISK GROUP ANALYSIS ====================
  observe({
    req(rv$survival_analysis_data, rv$survival_results, input$scores_survival)
    
    cutoff_pct <- input$high_risk_cutoff / 100
    data <- rv$survival_analysis_data
    
    # Use first selected score for high-risk analysis
    score_col <- input$scores_survival[1]
    
    if (exists("get_high_risk_summary")) {
      output$high_risk_summary <- renderText({
        summary_stats <- get_high_risk_summary(data, score_col, cutoff_pct)
        paste(capture.output(print(summary_stats)), collapse = "\n")
      })
    } else {
      # Fallback:  simple summary
      output$high_risk_summary <- renderText({
        score_vec <- data[[score_col]]
        threshold <- quantile(score_vec, cutoff_pct, na.rm = TRUE)
        high_risk <- score_vec >= threshold
        
        paste0("High-risk cutoff (", input$high_risk_cutoff, "th percentile): ", round(threshold, 3), "\n",
               "High-risk group: ", sum(high_risk, na.rm = TRUE), " individuals (", 
               round(100 * mean(high_risk, na.rm = TRUE), 1), "%)\n",
               "Low-risk group: ", sum(! high_risk, na.rm = TRUE), " individuals (", 
               round(100 * mean(! high_risk, na.rm = TRUE), 1), "%)\n",
               "\nHigh-risk group events: ", sum(data$event_censored[high_risk], na.rm = TRUE), "/", sum(high_risk, na.rm = TRUE),
               "\nLow-risk group events:  ", sum(data$event_censored[!high_risk], na.rm = TRUE), "/", sum(!high_risk, na.rm = TRUE))
      })
    }
    
    # High-risk KM plot
    output$high_risk_plot <- renderPlot({
      req(exists("plot_survival_curves"))
      score_vec <- data[[score_col]]
      threshold <- quantile(score_vec, cutoff_pct, na. rm = TRUE)
      data[, risk_group := ifelse(get(score_col) >= threshold, "High Risk", "Low Risk")]
      
      # Plot stratified by risk group
      library(survival)
      library(survminer)
      fit <- survfit(Surv(time_censored, event_censored) ~ risk_group, data = data)
      
      # Use ggsurvplot if available, otherwise base plot
      if (requireNamespace("survminer", quietly = TRUE)) {
        survminer::ggsurvplot(fit, data = data, 
                              risk.table = TRUE, 
                              pval = TRUE,
                              title = paste0("Survival by Risk Group (cutoff:  ", input$high_risk_cutoff, "th percentile)"),
                              xlab = "Time (years)",
                              ylab = "Survival Probability")
      } else {
        plot(fit, col = c("blue", "red"), lwd = 2, 
             xlab = "Time (years)", ylab = "Survival Probability",
             main = paste0("Survival by Risk Group (cutoff: ", input$high_risk_cutoff, "th percentile)"))
        legend("topright", legend = c("High Risk", "Low Risk"), col = c("red", "blue"), lwd = 2)
      }
    })
  })
  
  observeEvent(input$save_survival_btn, {
    req(rv$survival_results)
    outp <- trimws(input$survival_output_path %||% "")
    if (!nzchar(outp)) {
      showNotification("Provide output path first.", type = "warning")
      return()
    }
    tryCatch({
      dir. create(dirname(outp), recursive = TRUE, showWarnings = FALSE)
      if (is.data.frame(rv$survival_results) || is.data.table(rv$survival_results)) {
        fwrite(rv$survival_results, outp)
      } else if (is.list(rv$survival_results) && !is.null(rv$survival_results$summary)) {
        fwrite(rv$survival_results$summary, outp)
      }
      showNotification(paste("Saved:", outp), type = "message")
    }, error = function(e) {
      showNotification(paste("Save error:", e$message), type = "error")
    })
  })
}

shinyApp(ui, server)
