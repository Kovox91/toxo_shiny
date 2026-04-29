# Plot generation module
library(ggplot2)
library(plotly)
library(dplyr)

theme_shiny <- function() {
  theme_bw(base_family = "Liberation Sans") +
    theme(
      text = element_text(family = "Liberation Sans", colour = "black"),
      
      axis.title = element_text(size = 11),
      axis.text = element_text(size = 10),
      
      legend.title = element_text(size = 10, face = "bold"),
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.text = element_text(size = 9),
      legend.box.spacing = unit(3, "mm"),
      legend.key.width = unit(6, "mm"),
      
      plot.title = element_text(size = 12, face = "bold"),
      plot.margin = margin(3, 3, 3, 3, "mm"),
      
      panel.border = element_blank(),
      panel.grid.major = element_line(color = "grey90"),
      panel.grid.minor = element_blank(),
      
      strip.text = element_text(size = 10, margin = margin(b = 2)),
      strip.background = element_rect(
        fill = "grey95",
        colour = NA
      )
    )
}

#' Generate polyA length plot
#'
#' @param data Filtered dataset from load_all_data
#' @param better_names_ref Reference dataframe with cryoem info
#' @param selected_genes Vector of genes to plot
#' @param compare_all If TRUE, facet by experiment; if FALSE, single dataset
#' @param interactive If TRUE, returns plotly object; if FALSE, returns ggplot
#'
#' @return ggplot or plotly object
#'
make_polyA_plot <- function(data, better_names_ref, selected_genes, compare_all = FALSE, interactive = TRUE) {
  
  if (nrow(data) == 0) {
    return(
      ggplot() +
        annotate("text", x = 0.5, y = 0.5, label = "No data to display",
                 size = 6) +
        theme_void()
    )
  }
  
  # Create base plot
  p <- ggplot(data,
              aes(
                x = polyA,
                y = n_mean,
                color = con_unified
              )) +
    geom_line(size = 0.5) +
    theme_shiny() +
    labs(
      x = "Poly(A) tail length (nt)",
      y = "Relative Read Frequency (%)",
      color = "Condition"
    )
  
    p <- p +
      facet_grid(experiment ~ Geneid,
                 scales = "free")
  
  p <- p +
    scale_x_continuous(
      breaks = function(x) seq(0, ceiling(x[2]), by = 2),
      labels = function(x) ifelse(x %% 2 == 0, x, "")
    )
  
  # Add reference lines only for genes that are actually plotted
  genes_plotted <- unique(data$Geneid)
  facet_keys <- data %>%
    distinct(experiment, Geneid)
  
  # Add reference lines for Shikha2025
  shikha_refs <- better_names_ref %>%
    filter(gtf_rRNA_id %in% genes_plotted & cryoem_source == "Shikha2025") %>%
    transmute(
      Geneid = gtf_rRNA_id,
      xintercept = cryoem_bases_trimmed - 0.05
    ) %>%
    distinct() %>%
    inner_join(facet_keys, by = "Geneid")
  
  if (nrow(shikha_refs) > 0) {
    p <- p +
      geom_vline(
        data = shikha_refs,
        aes(xintercept = xintercept),
        color = "black",
        linetype = "dashed",
        size = 0.3,
        inherit.aes = FALSE
      )
  }
  
  # Add reference lines for Wang2024
  wang_refs <- better_names_ref %>%
    filter(gtf_rRNA_id %in% genes_plotted & cryoem_source == "Wang2024") %>%
    transmute(
      Geneid = gtf_rRNA_id,
      xintercept = cryoem_bases_trimmed + 0.05
    ) %>%
    distinct() %>%
    inner_join(facet_keys, by = "Geneid")
  
  if (nrow(wang_refs) > 0) {
    p <- p +
      geom_vline(
        data = wang_refs,
        aes(xintercept = xintercept),
        color = "#009e73ff",
        linetype = "dashed",
        size = 0.3,
        inherit.aes = FALSE
      )
  }
  
  if (interactive) {
    p <- ggplotly(p, tooltip = c("x", "y", "colour")) %>%
      layout(
        hovermode = "closest",
        font = list(family = "Liberation Sans")
      )
  }
  
  return(p)
}

#' Create summary table for filtered data
#'
#' @param data Filtered dataset
#'
#' @return DT table
#'
make_summary_table <- function(data) {
  summary <- data %>%
    group_by(Geneid, con_unified, experiment) %>%
    summarise(
      n_measurements = n(),
      mean_polyA = round(mean(polyA, na.rm = TRUE), 2),
      max_freq = round(max(n_mean, na.rm = TRUE), 2),
      .groups = "drop"
    ) %>%
    arrange(Geneid, experiment)
  
  return(summary)
}

#' Generate ridge plot showing polyA length distributions
#'
#' @param data Filtered dataset
#' @param interactive If TRUE, returns plotly object; if FALSE, returns ggplot
#'
#' @return ggplot or plotly object
#'
make_ridge_plot <- function(data, interactive = TRUE) {
  
  if (nrow(data) == 0) {
    return(
      ggplot() +
        annotate("text", x = 0.5, y = 0.5, label = "No data to display",
                 size = 6) +
        theme_void()
    )
  }

  ridge_source <- data %>%
    filter(!is.na(con_unified), !is.na(polyA), !is.na(n_mean)) %>%
    filter(
      (experiment %in% c("rmRHel1", "rmRSiC") & !con_unified %in% c("PolyA-KO", "Untreated")) |
      (experiment == "Shikha et al. 2025" & con_unified %in% c("PolyA-KO", "Untreated")) |
      (!experiment %in% c("rmRHel1", "rmRSiC", "Shikha et al. 2025"))
    )

  if (nrow(ridge_source) == 0) {
    return(
      ggplot() +
        annotate("text", x = 0.5, y = 0.5, label = "No data after condition filtering",
                 size = 6) +
        theme_void()
    )
  }

  # Build offset lines: each condition gets a baseline offset, then n_mean is
  # scaled and added on top to preserve the original line shape.
  condition_levels <- ridge_source %>%
    pull(con_unified) %>%
    as.character() %>%
    unique() %>%
    sort()

  if (length(condition_levels) == 0) {
    return(
      ggplot() +
        annotate("text", x = 0.5, y = 0.5, label = "No condition data to display",
                 size = 6) +
        theme_void()
    )
  }

  ridge_data <- ridge_source %>%
    group_by(experiment, Geneid, con_unified, polyA) %>%
    summarise(
      n_mean = mean(n_mean, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(experiment_index = dense_rank(experiment) - 1) %>%
    group_by(experiment) %>%
    mutate(
      condition_index_local = dense_rank(con_unified) - 1
    ) %>%
    ungroup() %>%
    mutate(
      con_unified = factor(con_unified, levels = condition_levels),
      condition_index = experiment_index * 100 + condition_index_local
    ) %>%
    group_by(experiment, Geneid) %>%
    mutate(
      panel_max = max(n_mean, na.rm = TRUE),
      panel_max = ifelse(panel_max <= 0, 1, panel_max),
      y_offset = condition_index + (n_mean / panel_max) * 3,
      hover_text = paste0(
        "Gene: ", Geneid,
        "<br>Dataset: ", experiment,
        "<br>Condition: ", as.character(con_unified),
        "<br>Poly(A): ", round(polyA, 2),
        "<br>Relative freq (%): ", round(n_mean, 2)
      )
    ) %>%
    arrange(experiment, Geneid, con_unified, polyA) %>%
    ungroup()

  y_breaks <- ridge_data %>%
    distinct(experiment, condition_index, con_unified) %>%
    mutate(y_label = as.character(con_unified))

  y_label_map <- y_breaks %>%
    distinct(condition_index, y_label) %>%
    arrange(condition_index)

  y_break_values <- y_label_map$condition_index
  names(y_break_values) <- y_label_map$y_label

  p <- ggplot(ridge_data,
              aes(
                x = polyA,
                y = y_offset,
                color = con_unified,
                fill = con_unified,
                group = interaction(experiment, Geneid, con_unified),
                text = hover_text
              )) +
    geom_hline(
      data = y_breaks,
      aes(yintercept = condition_index),
      color = "grey88",
      linewidth = 0.25,
      inherit.aes = FALSE
    ) +
    geom_ribbon(
      aes(ymin = condition_index, ymax = y_offset),
      alpha = 0.2,
      color = NA,
      show.legend = FALSE
    ) +
    geom_line(linewidth = 0.35) +
    theme_shiny() +
    scale_y_continuous(
      breaks = y_break_values,
      labels = names(y_break_values),
      expand = expansion(mult = c(0.02, 0.06))
    ) +
    scale_x_continuous(
      breaks = function(x) seq(0, ceiling(x[2]), by = 2),
      labels = function(x) ifelse(x %% 2 == 0, x, "")
    ) +
    labs(
      x = "Poly(A) tail length (nt)",
      y = "Condition (offset baselines)",
      color = "Condition"
    ) +
    theme(
      legend.position = "none",
      axis.text.y = element_text(size = 8)
    )

    p <- p + facet_grid(experiment ~ Geneid, scales = "free")

  
  if (interactive) {
    p <- ggplotly(p, tooltip = "text") %>%
      layout(
        hovermode = "closest",
        font = list(family = "Liberation Sans")
      )
  }
  
  return(p)
}
