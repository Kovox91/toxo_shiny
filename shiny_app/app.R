library(shiny)
library(tidyverse)
library(plotly)
library(DT)

# Source helper functions
source("R/data_loading.R")
source("R/plot_module.R")

# Load all data once at startup
data_env <- load_all_data(base_path = ".")
all_data <- data_env$data
better_names_ref <- data_env$better_names_ref
stats_data <- data_env$stats_data

# Get unique values for filters
experiments <- sort(unique(all_data$experiment))
all_genes <- sort(unique(all_data$Geneid))

# ============================================================================
# UI
# ============================================================================

ui <- fluidPage(
  titlePanel("Toxoplasma Poly(A) Tail Analysis"),
  sidebarLayout(
    sidebarPanel(
      width = 3,

      # Data source selector
      h4("1. Select Datasets"),
      uiOutput("dataset_checkboxes"),
      hr(),

      # Gene filter
      h4("2. Select Genes"),
      selectizeInput(
        inputId = "gene_select",
        label = "Genes:",
        choices = NULL, # Will be populated reactively
        multiple = TRUE,
        options = list(
          maxItems = 6,
          plugins = list("remove_button")
        )
      ),
      p(
        "Select 1 to 6 genes. Empty selection is not allowed.",
        class = "text-muted"
      ),
      hr(),

      # Condition filter
      h4("3. Filter Conditions"),
      uiOutput("condition_checkboxes"),
      hr(),

      # Summary info
      h5("Dataset Info:"),
      uiOutput("data_summary")
    ),
    mainPanel(
      width = 9,

      # Tabs for different outputs
      tabsetPanel(
        tabPanel(
          "Line Plot",
          uiOutput("plot_container"),
          br(),
          p(
            "Dashed lines: Cryo EM derived Poly(A) tail lengths according to Wang et al. 2024 (teal) and Shikha et al. 2025 (black)",
            class = "text-muted"
          )
        ),
        tabPanel(
          "Ridge Plot",
          br(),
          plotlyOutput("ridge_plot", height = "900px")
        ),
        tabPanel(
          "Statistics",
          br(),
          fluidRow(
            column(
              width = 4,
              uiOutput("stats_cond1_checkboxes")
            ),
            column(
              width = 4,
              uiOutput("stats_cond2_checkboxes")
            ),
            column(
              width = 4,
              numericInput(
                "stats_fdr_max",
                "Max FDR:",
                value = 0.05
              )
            )
          ),
          fluidRow(
            column(
              width = 4,
              numericInput(
                "stats_wasserstein_min",
                "Min Wasserstein:",
                value = 1
              )
            ),
            column(
              width = 4,
              numericInput(
                "stats_wnorm_min",
                "Min W_Norm:",
                value = 1.5
              )
            ),
            column(
              width = 4,
              p("Statistics table uses currently selected datasets from the sidebar.",
                class = "text-muted"
              )
            )
          ),
          DT::dataTableOutput("stats_table"),
          br(),
          wellPanel(
            h5("About the statistics"),
            p(HTML("<b>Wasserstein distance</b> (Earth Mover's Distance) measures how much the poly(A) tail length distribution of one condition must be \"moved\" to match another. Larger values indicate greater distributional differences between conditions.")),
            p(HTML("<b>FDR</b> is computed from a permutation test: reads are randomly redistributed between two conditions 1,000 times to generate a null distribution of Wasserstein distances, and the observed value is ranked within this null to compute a p-value, which is then FDR-corrected across all gene–condition pair tests.")),
            p(HTML("<b>p95_W</b> is the 95th percentile of the Wasserstein distance between technical replicates within a condition for the respective gene. It represents the expected technical background noise.")),
            p(HTML("<b>W_Norm</b> is the fold increase of the Wasserstein distance between conditions over this technical background (p95_W). Values greater than 1 indicate the biological difference exceeds the technical noise floor."))
          )
        )
      )
    )
  )
)

# ============================================================================
# SERVER
# ============================================================================

server <- function(input, output, session) {
  # Reactive: Generate dataset checkboxes
  output$dataset_checkboxes <- renderUI({
    checkboxGroupInput(
      inputId = "selected_datasets",
      label = "Datasets:",
      choices = experiments,
      selected = experiments, # Default: all selected
      inline = FALSE
    )
  })

  # Reactive: Get selected datasets
  selected_datasets <- reactive({
    req(input$selected_datasets)
    input$selected_datasets
  })

  # Reactive: Filter data by selected datasets
  filtered_by_source <- reactive({
    all_data %>%
      filter(experiment %in% selected_datasets())
  })

  # Reactive: Update gene choices only when dataset selection changes.
  # This avoids input-update feedback loops that can cause repeated re-rendering.
  observeEvent(selected_datasets(),
    {
      genes_in_source <- sort(unique(filtered_by_source()$Geneid))

      # Default genes to select
      default_genes <- c("LSUC", "RNA33_1", "RNA34_F", "RNA42")

      # Isolate avoids making this observer depend on gene_select itself.
      current_selection <- isolate(input$gene_select)

      if (is.null(current_selection)) {
        available_defaults <- intersect(default_genes, genes_in_source)
        if (length(available_defaults) > 0) {
          persistent_selection <- available_defaults
        } else {
          persistent_selection <- genes_in_source[1:min(3, length(genes_in_source))]
        }
      } else {
        # Keep only genes that still exist in current dataset selection.
        # Empty remains empty (invalid selection handled by validate() in outputs).
        persistent_selection <- intersect(current_selection, genes_in_source)
      }

      if (length(persistent_selection) > 6) {
        persistent_selection <- persistent_selection[1:6]
      }

      updateSelectizeInput(
        session,
        inputId = "gene_select",
        choices = genes_in_source,
        selected = persistent_selection
      )
    },
    ignoreInit = FALSE
  )

  # Reactive: Filter data by selected genes
  filtered_by_gene <- reactive({
    selected_genes <- input$gene_select

    if (is.null(selected_genes) || length(selected_genes) == 0) {
      # Empty selection is intentionally invalid for plotting.
      filtered_by_source() %>%
        slice(0)
    } else {
      filtered_by_source() %>%
        filter(Geneid %in% selected_genes)
    }
  })

  # Reactive: Get available conditions for this source & genes
  available_conditions <- reactive({
    sort(unique(na.omit(filtered_by_gene()$con_unified)))
  })

  # Reactive: Get selected conditions
  selected_conditions <- reactive({
    if (is.null(input$conditions) || length(input$conditions) == 0) {
      available_conditions()
    } else {
      input$conditions
    }
  })

  # Reactive: Generate condition checkboxes
  output$condition_checkboxes <- renderUI({
    conditions <- available_conditions()

    # Create checkboxes, pre-select all
    checkboxGroupInput(
      inputId = "conditions",
      label = "Conditions:",
      choices = conditions,
      selected = conditions,
      inline = FALSE
    )
  })

  # Reactive: Final filtered data
  plot_data <- reactive({
    filtered_by_gene() %>%
      filter(con_unified %in% selected_conditions())
  })

  # Reactive: Selected genes for plotting (for reference lines)
  selected_genes_for_plot <- reactive({
    if (is.null(input$gene_select) || length(input$gene_select) == 0) {
      character(0)
    } else {
      input$gene_select
    }
  })

  # Reactive: Calculate plot height based on number of genes and datasets
  plot_height <- reactive({
    n_genes <- length(unique(plot_data()$Geneid))
    n_datasets <- length(unique(plot_data()$experiment))

    if (n_datasets > 1) {
      # Stack experiments: ~150px per row per gene, ~100px per dataset
      n_rows <- ceiling(n_genes / 3)
      height_px <- n_rows * 300 + n_datasets * 150
    } else {
      # Single dataset mode
      n_rows <- ceiling(n_genes / 3)
      height_px <- n_rows * 250
    }

    # Ensure minimum height
    max(600, height_px)
  })

  # Output: Plot container with dynamic height
  output$plot_container <- renderUI({
    plotlyOutput("polyA_plot", height = paste0(plot_height(), "px"))
  })

  # Output: Main plot
  output$polyA_plot <- renderPlotly({
    validate(
      need(
        !is.null(input$gene_select) && length(input$gene_select) > 0,
        "Please select at least 1 gene."
      ),
      need(
        length(input$gene_select) <= 6,
        "Too many genes selected. Please select at most 6 genes."
      )
    )

    make_polyA_plot(
      data = plot_data(),
      better_names_ref = better_names_ref,
      selected_genes = selected_genes_for_plot(),
      compare_all = length(unique(plot_data()$experiment)) > 1,
      interactive = TRUE
    )
  })

  # Output: Ridge plot
  output$ridge_plot <- renderPlotly({
    validate(
      need(
        !is.null(input$gene_select) && length(input$gene_select) > 0,
        "Please select at least 1 gene."
      ),
      need(
        length(input$gene_select) <= 6,
        "Too many genes selected. Please select at most 6 genes."
      )
    )

    make_ridge_plot(
      data = plot_data(),
      interactive = TRUE
    )
  })

  # Output: Dataset info box
  output$data_summary <- renderUI({
    data <- plot_data()

    HTML(paste0(
      "<b>Genes shown:</b> ", length(unique(data$Geneid)), "<br>",
      "<b>Conditions:</b> ", paste(selected_conditions(), collapse = ", "), "<br>",
      "<b>Data points:</b> ", nrow(data)
    ))
  })

  # Statistics tab: constrain to selected datasets from sidebar
  stats_by_dataset <- reactive({
    req(stats_data)
    stats_data %>%
      filter(Experiment %in% selected_datasets())
  })

  # Render Condition1 checkboxes
  output$stats_cond1_checkboxes <- renderUI({
    sdata <- stats_by_dataset()
    cond1_vals <- sort(unique(sdata$Condition1))

    checkboxGroupInput(
      "stats_cond1_filter",
      "Filter Condition1:",
      choices = cond1_vals,
      selected = cond1_vals,
      inline = FALSE
    )
  })

  # Render Condition2 checkboxes
  output$stats_cond2_checkboxes <- renderUI({
    sdata <- stats_by_dataset()
    cond2_vals <- sort(unique(sdata$Condition2))

    checkboxGroupInput(
      "stats_cond2_filter",
      "Filter Condition2:",
      choices = cond2_vals,
      selected = cond2_vals,
      inline = FALSE
    )
  })

  stats_filtered <- reactive({
    sdata <- stats_by_dataset()

    # Convert numeric inputs, using defaults if NULL or invalid
    max_fdr <- as.numeric(input$stats_fdr_max)
    if (is.na(max_fdr)) max_fdr <- 0.05

    min_wasserstein <- as.numeric(input$stats_wasserstein_min)
    if (is.na(min_wasserstein)) min_wasserstein <- 1

    min_wnorm <- as.numeric(input$stats_wnorm_min)
    if (is.na(min_wnorm)) min_wnorm <- 1.5

    if (!is.null(input$stats_cond1_filter) && length(input$stats_cond1_filter) > 0) {
      sdata <- sdata %>% filter(Condition1 %in% input$stats_cond1_filter)
    }

    if (!is.null(input$stats_cond2_filter) && length(input$stats_cond2_filter) > 0) {
      sdata <- sdata %>% filter(Condition2 %in% input$stats_cond2_filter)
    }

    sdata %>%
      filter(
        Wasserstein >= min_wasserstein,
        FDR <= max_fdr,
        W_Norm >= min_wnorm
      ) %>%
      arrange(FDR, desc(W_Norm))
  })

  output$stats_table <- DT::renderDataTable({
    DT::datatable(
      stats_filtered(),
      options = list(
        pageLength = 20,
        scrollX = TRUE,
        autoWidth = TRUE
      ),
      rownames = FALSE
    )
  })
}

# ============================================================================
# RUN APP
# ============================================================================

shinyApp(ui = ui, server = server)
