# Toxoplasma Poly(A) Tail Analysis Shiny App

A Shiny web application for interactive visualization and exploration of poly(A) tail length analysis data from your dissertation.

## Features

✓ **Data source selection**: Choose between Helicase, HPR, Shikha et al., or Wang 2024 datasets  
✓ **Dynamic gene filtering**: Multi-select up to 55 genes with live preview  
✓ **Condition filtering**: Toggle conditions on/off with instant plot updates  
✓ **Interactive plots**: Zoom, hover, and download plots with plotly  
✓ **Summary statistics**: View aggregated stats by gene/condition/dataset  
✓ **Reference lines**: Cryoem structural reference data overlaid (Shikha & Wang)  

## Running the App

### Prerequisites

Make sure you have the required R packages installed:

```r
install.packages(c("shiny", "tidyverse", "plotly", "DT"))
```

### Launch

From the `shiny_app/` directory:

```r
shiny::runApp()
```

Or from the project root:

```r
shiny::runApp("shiny_app/")
```

The app will open in your default browser at `http://localhost:3838` (or similar).

## File Structure

```
shiny_app/
├── app.R                 # Main Shiny application
├── R/
│   ├── data_loading.R    # Data loading & preparation
│   ├── plot_module.R     # Plotting functions
│   └── ...               # Add more modules as you expand
├── www/                  # Static web assets (CSS, images, etc.)
└── README.md            # This file
```

## How It Works

1. **Sidebar controls** (left panel):
   - Select a data source (experiment)
   - Choose genes to display
   - Filter by experimental conditions

2. **Main panel** (right side):
   - **Plot tab**: Interactive line plot of polyA length vs read frequency
   - **Summary tab**: Aggregated statistics
   - **Data tab**: Full raw data table

3. **Instant reactivity**: All plots and tables update as soon as you change any filter

## Next Steps for Expansion

To add more analyses (as you expand into a full "dissertation data browser"):

1. Create new R modules in `R/` folder (e.g., `R/coverage_module.R`)
2. Add corresponding data loading functions
3. Create new `tabPanel()` tabs in the UI
4. Add reactive filters as needed

Example structure for a new analysis:
```
tabPanel(
  "New Analysis",
  plotlyOutput("new_plot"),
  DT::dataTableOutput("new_table")
)
```

## Deployment

To deploy to shinyapps.io:

```r
# Install deployment tools
install.packages("rsconnect")

# Deploy (from shiny_app/ directory)
rsconnect::deployApp()
```

## Notes

- Data is loaded once at startup for performance
- All reactions are reactive, not observeEvent() — provides instant updates
- Reference lines (red = Shikha2025, black = Wang2024) update based on selected genes
