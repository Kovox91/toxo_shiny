# Data loading and preparation for Shiny app
library(tidyverse)

load_all_data <- function(base_path = ".") {
  
  # Load datasets
  helicase <- readRDS(file.path(base_path, "input/mrhel1_polyA.RDS")) %>%
    select(everything(), -n) %>%
    mutate(experiment = "rmRHel1")
  
  hpr <- readRDS(file.path(base_path, "input/mrsic_polyA.RDS")) %>%
    select(everything(), "con" = condition) %>%
    mutate(experiment = "rmRSiC")

  # Load pre-curated statistics data
  stats_file <- file.path(base_path, "input/polyA_stats_curated.rds")
  if (!file.exists(stats_file)) {
    stop(
      "Missing curated stats file: input/polyA_stats_curated.rds\n",
      "Run scripts/build_curated_stats.R once to generate it."
    )
  }
  stats_data <- readRDS(stats_file)
  
  # Load references
  names_raw <- read_csv(
    file.path(base_path, "reference/naming_help.csv"),
    show_col_types = FALSE
  )
  
  polyA_shikha_raw <- read.delim2(
    file.path(base_path, "reference/shikha_polyA.tsv")
  )
  
  polyA_wang_raw <- read_csv(
    file.path(base_path, "reference/wang_polyA_lengths.csv"),
    skip = 2,
    show_col_types = FALSE
  )
  
  better_names_raw <- read.delim2(
    file.path(base_path, "reference/cryoem_to_gtf_mapping.tsv")
  )
  
  # Prepare names reference
  names <- better_names_raw %>%
    select(
      "name" = cryoem_fragment,
      "Geneid" = gtf_rRNA_id,
      "ID_long" = gtf_rRNA_id,
      "Source" = cryoem_source
    ) %>%
    mutate(name = base::gsub("-", "", name)) %>%
    filter(!Geneid %in% c("No_GTF_Match", "SSU-X"))
  
  # Prepare Shikha data
  polyA_shikha <- polyA_shikha_raw %>%
    select(
      "name" = Nomenclature,
      type,
      "polyA" = polya_length,
      "n_mean" = polya_pct,
      "polyA_struct" = polyA
    ) %>%
    mutate(
      experiment = "Shikha et al. 2025",
      name = base::gsub("-", "", name),
      n_mean = as.numeric(n_mean)
    ) %>%
    left_join(filter(names, Source == "Shikha2025"), by = "name",
              relationship = "many-to-many") %>%
    select("Geneid" = ID_long, "con" = type, polyA, n_mean, experiment, polyA_struct, name)
  
  # Prepare Wang data
  polyA_wang <- polyA_wang_raw %>%
    select("name" = Nomenclature,
           seq = `rRNA Sequence (5'-3')`) %>%
    mutate(
      name = base::gsub("-", "", name),
      polya_length = nchar(sub(".*?(A+)$", "\\1", seq)),
      polya_length = ifelse(grepl("A+$", seq), polya_length, 0)
    ) %>%
    select(-seq) %>%
    left_join(filter(names, Source == "Wang2024"), by = "name",
              relationship = "many-to-many") %>%
    mutate(
      ID_long = ifelse(is.na(ID_long), name, ID_long),
      ID_long = ifelse(ID_long == "RNA37/38_F", "RNA37-38_F", Geneid)
    )
  
  # Combine all datasets
  all_data <- rbind(helicase, hpr) %>%
    mutate(n_mean = 100 * n_mean) %>%
    rbind(polyA_shikha) %>%
    mutate(
      con_unified = case_when(
        # controls
        con == "Untreated" ~ "Untreated",
        con == "PolyA-KO" ~ "PolyA-KO",
        
        # parental (mR*)
        con %in% c("mRHel-FLAG -ATc", "mRSiC-FLAG -ATc") ~ "parental 0d",
        con %in% c("mRHel-FLAG 3d ATc") ~ "parental 3d",
        
        # promoter replaced (rmR*) — separate days
        con %in% c("rmRHel-FLAG -ATc", "rmRHel-FLAG 0d ATc",
                   "rmRSiC-FLAG -ATc") ~ "PR 0d",
        
        con %in% c("rmRHel-FLAG 1d ATc", "rmRSiC-FLAG 1d ATc") ~ "PR 1d",
        con %in% c("rmRHel-FLAG 2d ATc", "rmRSiC-FLAG 2d ATc") ~ "PR 2d",
        con %in% c("rmRHel-FLAG 3d ATc", "rmRSiC-FLAG 3d ATc") ~ "PR 3d",
        
        TRUE ~ NA_character_
      ),
      Geneid = ifelse(is.na(Geneid), name, Geneid)
    )
  
  return(list(
    data = all_data,
    better_names_ref = better_names_raw,
    stats_data = stats_data
  ))
}
