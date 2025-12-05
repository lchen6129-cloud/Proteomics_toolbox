library(shiny)
library(shinydashboard)
library(DT)
library(readr)
library(readxl)
library(dplyr)
library(tidyr)

# Source utility functions
source("R/utils.R")
source("R/model_comparison.R")
source("R/score_generation.R")
source("R/phenotype_association.R")

# UI
ui <- dashboardPage(
  dashboardHeader(title = "Protein Model Analysis"),
  dashboardSidebar(
    sidebarMenu(
      menuItem("Data Loading", tabName = "data_loading", icon = icon("folder-open")),
      menuItem("Model Comparison", tabName = "comparison", icon = icon("exchange-alt")),
      menuItem("Score Generation", tabName = "scores", icon = icon("calculator")),
      menuItem("Phenotype Association", tabName = "phenotype", icon = icon("microscope"))
    )
  ),
  dashboardBody(
    tabItems(
      # Data Loading Tab
      tabItem(tabName = "data_loading",
        fluidPage(
          h2("Data Loading"),
          p("Provide file paths to your data files (CSV or Excel format)"),
          br(),
          box(
            title = "File Paths",
            status = "primary",
            solidHeader = TRUE,
            width = 12,
            textInput("path_models", "Path to Models Table (Table 1):", 
                     placeholder = "/path/to/models_table.csv"),
            textInput("path_proteins", "Path to Protein Levels Table (Table 2):", 
                     placeholder = "/path/to/protein_levels.csv"),
            textInput("path_phenotype", "Path to Phenotype Data (Table 3):", 
                     placeholder = "/path/to/phenotype_data.csv"),
            br(),
            actionButton("load_data_btn", "Load Data", class = "btn-primary", icon = icon("play")),
            br(), br(),
            uiOutput("data_load_message")
          ),
          br(),
          fluidRow(
            column(6,
              box(
                title = "Models Table Preview",
                status = "info",
                solidHeader = TRUE,
                width = 12,
                DT::dataTableOutput("preview_models")
              )
            ),
            column(6,
              box(
                title = "Protein Levels Table Preview",
                status = "info",
                solidHeader = TRUE,
                width = 12,
                DT::dataTableOutput("preview_proteins")
              )
            )
          ),
          fluidRow(
            column(12,
              box(
                title = "Phenotype Data Preview",
                status = "info",
                solidHeader = TRUE,
                width = 12,
                DT::dataTableOutput("preview_phenotype")
              )
            )
          )
        )
      ),
      
      # Model Comparison Tab
      tabItem(tabName = "comparison",
        fluidPage(
          h2("Model Comparison"),
          br(),
          box(
            title = "Select Models to Compare",
            status = "primary",
            solidHeader = TRUE,
            width = 12,
            fluidRow(
              column(4, uiOutput("select_model1")),
              column(4, uiOutput("select_model2")),
              column(4, 
                br(),
                actionButton("compare_btn", "Compare", class = "btn-info", icon = icon("balance-scale"))
              )
            )
          ),
          br(),
          box(
            title = "Comparison Results",
            status = "success",
            solidHeader = TRUE,
            width = 12,
            tabsetPanel(
              tabPanel("Summary",
                br(),
                fluidRow(
                  column(4, 
                    valueBoxOutput("overlap_count")
                  ),
                  column(4,
                    valueBoxOutput("unique_model1")
                  ),
                  column(4,
                    valueBoxOutput("unique_model2")
                  )
                ),
                br(),
                h4("Overlapping Proteins:"),
                DT::dataTableOutput("overlap_proteins")
              ),
              tabPanel("Detailed Comparison",
                br(),
                DT::dataTableOutput("detailed_comparison")
              )
            )
          ),
          br(),
          box(
            title = "Save Results",
            status = "warning",
            solidHeader = TRUE,
            width = 12,
            fluidRow(
              column(10,
                textInput("output_comparison", "Output File Path:", 
                         placeholder = "/path/to/output/comparison_results.csv")
              ),
              column(2,
                br(),
                actionButton("save_comparison", "Save", class = "btn-warning", icon = icon("save"))
              )
            ),
            br(),
            uiOutput("comparison_save_message")
          )
        )
      ),
      
      # Score Generation Tab
      tabItem(tabName = "scores",
        fluidPage(
          h2("Score Generation"),
          br(),
          box(
            title = "Select Models and Generate Scores",
            status = "primary",
            solidHeader = TRUE,
            width = 12,
            fluidRow(
              column(8, uiOutput("select_models_score")),
              column(4,
                br(),
                actionButton("generate_scores_btn", "Generate Scores", class = "btn-info", icon = icon("cogs"))
              )
            )
          ),
          br(),
          box(
            title = "Score Results",
            status = "success",
            solidHeader = TRUE,
            width = 12,
            tabsetPanel(
              tabPanel("Scores",
                br(),
                DT::dataTableOutput("scores_table")
              ),
              tabPanel("Missing Variables Report",
                br(),
                DT::dataTableOutput("missing_vars_table")
              ),
              tabPanel("Summary",
                br(),
                verbatimTextOutput("scores_summary")
              )
            )
          ),
          br(),
          box(
            title = "Save Results",
            status = "warning",
            solidHeader = TRUE,
            width = 12,
            fluidRow(
              column(10,
                textInput("output_scores", "Output File Path:", 
                         placeholder = "/path/to/output/individual_scores.csv")
              ),
              column(2,
                br(),
                actionButton("save_scores", "Save", class = "btn-warning", icon = icon("save"))
              )
            ),
            br(),
            uiOutput("scores_save_message")
          )
        )
      ),
      
      # Phenotype Association Tab
      tabItem(tabName = "phenotype",
        fluidPage(
          h2("Phenotype Association Analysis"),
          br(),
          box(
            title = "Select Analysis Parameters",
            status = "primary",
            solidHeader = TRUE,
            width = 12,
            fluidRow(
              column(3, uiOutput("select_model_pheno")),
              column(3, 
                selectInput("regression_type", "Regression Type:",
                           choices = c("Linear" = "linear", 
                                      "Logistic" = "logistic", 
                                      "Cox" = "cox"))
              ),
              column(3, uiOutput("select_phenotypes")),
              column(3,
                br(),
                actionButton("run_regression", "Run Analysis", class = "btn-info", icon = icon("play"))
              )
            )
          ),
          br(),
          box(
            title = "Regression Results",
            status = "success",
            solidHeader = TRUE,
            width = 12,
            tabsetPanel(
              tabPanel("Results Table",
                br(),
                DT::dataTableOutput("regression_results")
              ),
              tabPanel("Model Statistics",
                br(),
                verbatimTextOutput("model_stats")
              )
            )
          ),
          br(),
          box(
            title = "Save Results",
            status = "warning",
            solidHeader = TRUE,
            width = 12,
            fluidRow(
              column(10,
                textInput("output_phenotype", "Output File Path:", 
                         placeholder = "/path/to/output/phenotype_association_results.csv")
              ),
              column(2,
                br(),
                actionButton("save_phenotype", "Save", class = "btn-warning", icon = icon("save"))
              )
            ),
            br(),
            uiOutput("phenotype_save_message")
          )
        )
      )
    )
  )
)

# Server
server <- function(input, output, session) {
  
  # Reactive values to store data
  data <- reactiveValues(
    models = NULL,
    proteins = NULL,
    phenotype = NULL,
    data_loaded = FALSE
  )
  
  # Data loading
  observeEvent(input$load_data_btn, {
    tryCatch({
      # Load models table
      data$models <- load_file(input$path_models)
      
      # Load protein levels table
      data$proteins <- load_file(input$path_proteins)
      
      # Load phenotype data
      data$phenotype <- load_file(input$path_phenotype)
      
      data$data_loaded <- TRUE
      
      output$data_load_message <- renderUI({
        div(class = "alert alert-success",
          h4("✓ Data loaded successfully!"),
          p(sprintf("Models: %d rows, %d columns", nrow(data$models), ncol(data$models))),
          p(sprintf("Proteins: %d rows, %d columns", nrow(data$proteins), ncol(data$proteins))),
          p(sprintf("Phenotype: %d rows, %d columns", nrow(data$phenotype), ncol(data$phenotype)))
        )
      })
    }, error = function(e) {
      output$data_load_message <- renderUI({
        div(class = "alert alert-danger",
          h4("✗ Error loading data:"),
          p(as.character(e$message))
        )
      })
    })
  })
  
  # Preview tables
  output$preview_models <- DT::renderDataTable({
    if(is.null(data$models)) return(NULL)
    DT::datatable(head(data$models, 10))
  })
  
  output$preview_proteins <- DT::renderDataTable({
    if(is.null(data$proteins)) return(NULL)
    DT::datatable(head(data$proteins, 10))
  })
  
  output$preview_phenotype <- DT::renderDataTable({
    if(is.null(data$phenotype)) return(NULL)
    DT::datatable(head(data$phenotype, 10))
  })
  
  # Model Comparison
  output$select_model1 <- renderUI({
    if(! data$data_loaded) return(NULL)
    selectInput("model1", "Model 1:", 
               choices = unique(data$models$model_name))
  })
  
  output$select_model2 <- renderUI({
    if(!data$data_loaded) return(NULL)
    selectInput("model2", "Model 2:", 
               choices = unique(data$models$model_name))
  })
  
  comparison_results <- reactiveValues(data = NULL)
  
  observeEvent(input$compare_btn, {
    tryCatch({
      results <- compare_models(data$models, input$model1, input$model2)
      comparison_results$data <- results
      
      output$overlap_count <- renderValueBox({
        valueBox(
          length(results$overlap_proteins),
          "Overlapping Proteins",
          icon = icon("circle"),
          color = "blue"
        )
      })
      
      output$unique_model1 <- renderValueBox({
        valueBox(
          results$unique_count_model1,
          paste("Unique in", input$model1),
          icon = icon("circle"),
          color = "green"
        )
      })
      
      output$unique_model2 <- renderValueBox({
        valueBox(
          results$unique_count_model2,
          paste("Unique in", input$model2),
          icon = icon("circle"),
          color = "red"
        )
      })
      
      output$overlap_proteins <- DT::renderDataTable({
        DT::datatable(results$overlap_data)
      })
      
      output$detailed_comparison <- DT::renderDataTable({
        DT::datatable(results$detailed_comparison)
      })
    }, error = function(e) {
      showNotification(paste("Error:", e$message), type = "error")
    })
  })
  
  # Save comparison results
  observeEvent(input$save_comparison, {
    if(is.null(comparison_results$data)) {
      output$comparison_save_message <- renderUI({
        div(class = "alert alert-danger", "No results to save.  Please run comparison first.")
      })
      return()
    }
    
    tryCatch({
      save_path <- input$output_comparison
      validate_output_path(save_path)
      
      # Write comparison summary and also save detailed table
      summary_df <- data.frame(
        Metric = c("Overlapping Proteins", 
                   paste("Unique in", input$model1),
                   paste("Unique in", input$model2)),
        Count = c(length(comparison_results$data$overlap_proteins),
                 comparison_results$data$unique_count_model1,
                 comparison_results$data$unique_count_model2)
      )
      
      # Create a basic combined CSV: write summary, then append overlap and detailed as separate CSVs
      # Simpler: write three separate files (summary, overlap, detailed) with suffixes
      base_path <- save_path
      if (tools::file_ext(base_path) == "") base_path <- paste0(base_path, ".csv")
      dir <- dirname(base_path)
      name_base <- tools::file_path_sans_ext(basename(base_path))
      
      summary_path <- file.path(dir, paste0(name_base, "_summary.csv"))
      overlap_path <- file.path(dir, paste0(name_base, "_overlap.csv"))
      detailed_path <- file.path(dir, paste0(name_base, "_detailed.csv"))
      
      write.csv(summary_df, summary_path, row.names = FALSE)
      write.csv(comparison_results$data$overlap_data, overlap_path, row.names = FALSE)
      write.csv(comparison_results$data$detailed_comparison, detailed_path, row.names = FALSE)
      
      output$comparison_save_message <- renderUI({
        div(class = "alert alert-success",
          h4("✓ Results saved successfully!"),
          p(paste("Files saved to:", summary_path)),
          p(paste("Overlap:", overlap_path)),
          p(paste("Detailed:", detailed_path))
        )
      })
    }, error = function(e) {
      output$comparison_save_message <- renderUI({
        div(class = "alert alert-danger",
          h4("✗ Error saving results:"),
          p(as.character(e$message))
        )
      })
    })
  })
  
  # Score Generation
  output$select_models_score <- renderUI({
    if(! data$data_loaded) return(NULL)
    selectInput("models_score", "Select Models:", 
               choices = unique(data$models$model_name),
               multiple = TRUE)
  })
  
  scores_results <- reactiveValues(data = NULL, warnings = NULL)
  
  observeEvent(input$generate_scores_btn, {
    tryCatch({
      results <- generate_scores(
        data$models, 
        data$proteins, 
        input$models_score
      )
      
      scores_results$data <- results$scores
      scores_results$warnings <- results$warnings
      
      output$scores_table <- DT::renderDataTable({
        DT::datatable(results$scores)
      })
      
      output$missing_vars_table <- DT::renderDataTable({
        DT::datatable(results$warnings)
      })
      
      output$scores_summary <- renderPrint({
        cat("Score Generation Summary\n")
        cat("========================\n")
        cat(sprintf("Number of individuals scored: %d\n", nrow(results$scores)))
        cat(sprintf("Number of models used: %d\n", length(input$models_score)))
        cat(sprintf("Individuals with warnings: %d\n", ifelse(nrow(results$warnings)>0, length(unique(results$warnings$individual_id)), 0)))
      })
    }, error = function(e) {
      showNotification(paste("Error:", e$message), type = "error")
    })
  })
  
  # Save scores results
  observeEvent(input$save_scores, {
    if(is.null(scores_results$data)) {
      output$scores_save_message <- renderUI({
        div(class = "alert alert-danger", "No results to save. Please generate scores first.")
      })
      return()
    }
    
    tryCatch({
      save_path <- input$output_scores
      validate_output_path(save_path)
      
      base_path <- save_path
      if (tools::file_ext(base_path) == "") base_path <- paste0(base_path, ".csv")
      dir <- dirname(base_path)
      name_base <- tools::file_path_sans_ext(basename(base_path))
      
      scores_path <- file.path(dir, paste0(name_base, "_scores.csv"))
      warnings_path <- file.path(dir, paste0(name_base, "_warnings.csv"))
      
      write.csv(scores_results$data, scores_path, row.names = FALSE)
      write.csv(scores_results$warnings, warnings_path, row.names = FALSE)
      
      output$scores_save_message <- renderUI({
        div(class = "alert alert-success",
          h4("✓ Results saved successfully!"),
          p(paste("Scores:", scores_path)),
          p(paste("Warnings:", warnings_path))
        )
      })
    }, error = function(e) {
      output$scores_save_message <- renderUI({
        div(class = "alert alert-danger",
          h4("✗ Error saving results:"),
          p(as.character(e$message))
        )
      })
    })
  })
  
  # Phenotype Association
  output$select_model_pheno <- renderUI({
    if(!data$data_loaded) return(NULL)
    selectInput("model_pheno", "Select Model:", 
               choices = unique(data$models$model_name))
  })
  
  output$select_phenotypes <- renderUI({
    if(is.null(data$phenotype)) return(NULL)
    # Exclude common ID/non-phenotype columns
    pheno_cols <- setdiff(names(data$phenotype), c("individual_id", "id"))
    selectInput("phenotypes", "Select Phenotypes:", 
               choices = pheno_cols,
               multiple = TRUE)
  })
  
  regression_results <- reactiveValues(data = NULL, stats = NULL)
  
  observeEvent(input$run_regression, {
    tryCatch({
      results <- run_phenotype_association(
        data$models,
        data$proteins,
        data$phenotype,
        input$model_pheno,
        input$phenotypes,
        input$regression_type
      )
      
      regression_results$data <- results$results_table
      regression_results$stats <- results$model_stats
      
      output$regression_results <- DT::renderDataTable({
        DT::datatable(results$results_table, options = list(scrollX = TRUE))
      })
      
      output$model_stats <- renderPrint({
        cat(results$model_stats)
      })
    }, error = function(e) {
      showNotification(paste("Error:", e$message), type = "error")
    })
  })
  
  # Save phenotype results
  observeEvent(input$save_phenotype, {
    if(is.null(regression_results$data)) {
      output$phenotype_save_message <- renderUI({
        div(class = "alert alert-danger", "No results to save. Please run analysis first.")
      })
      return()
    }
    
    tryCatch({
      save_path <- input$output_phenotype
      validate_output_path(save_path)
      
      base_path <- save_path
      if (tools::file_ext(base_path) == "") base_path <- paste0(base_path, ".csv")
      dir <- dirname(base_path)
      name_base <- tools::file_path_sans_ext(basename(base_path))
      
      results_path <- file.path(dir, paste0(name_base, "_regression_results.csv"))
      stats_path <- file.path(dir, paste0(name_base, "_model_stats.txt"))
      
      write.csv(regression_results$data, results_path, row.names = FALSE)
      writeLines(regression_results$stats, con = stats_path)
      
      output$phenotype_save_message <- renderUI({
        div(class = "alert alert-success",
          h4("✓ Results saved successfully!"),
          p(paste("Results:", results_path)),
          p(paste("Model summary:", stats_path))
        )
      })
    }, error = function(e) {
      output$phenotype_save_message <- renderUI({
        div(class = "alert alert-danger",
          h4("✗ Error saving results:"),
          p(as.character(e$message))
        )
      })
    })
  })
}

shinyApp(ui, server)
