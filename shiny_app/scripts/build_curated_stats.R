library(tidyverse)

base_path <- "."

normalize_condition <- function(x) {
  case_when(
    str_to_lower(x) %in% c("par", "par0") ~ "Parental 0",
    str_to_lower(x) == "par3" ~ "Parental 3",
    str_to_lower(x) %in% c("pr_0", "pr0") ~ "Promoter Replaced 0",
    str_detect(str_to_lower(x), "^pr_?[0-9]+$") ~
      paste("Promoter Replaced", str_replace(str_to_lower(x), "^pr_?", "")),
    TRUE ~ x
  )
}

# Function to assign cardinality order to conditions
get_condition_order <- function(x) {
  case_when(
    x == "Parental 0" ~ 1,
    x == "Promoter Replaced 0" ~ 2,
    x == "Parental 3" ~ 3,
    x == "Promoter Replaced 1" ~ 4,
    x == "Promoter Replaced 2" ~ 5,
    x == "Promoter Replaced 3" ~ 6,
    TRUE ~ 999
  )
}

hel_stats <- readRDS(file.path(base_path, "input/mrhel_polyA_stats.rds")) %>%
  mutate(experiment = "rmRHel1")

sic_stats <- readRDS(file.path(base_path, "input/mrsic_polyA_stats.rds")) %>%
  mutate(experiment = "rmRSiC")

stats_data <- bind_rows(hel_stats, sic_stats) %>%
  transmute(
    Experiment = experiment,
    Gene = Geneid,
    Condition1_temp = normalize_condition(cond1),
    Condition2_temp = normalize_condition(cond2),
    Wasserstein = round(as.numeric(wasserstein), 2),
    FDR = round(as.numeric(p_perm_fdr), 3),
    p95_W = round(as.numeric(p95_w), 2),
    W_Norm = round(as.numeric(W_norm), 1)
  ) %>%
  rowwise() %>%
  mutate(
    # Determine which condition should be first based on cardinality order
    order1 = get_condition_order(Condition1_temp),
    order2 = get_condition_order(Condition2_temp),
    Condition1 = if_else(order1 <= order2, Condition1_temp, Condition2_temp),
    Condition2 = if_else(order1 <= order2, Condition2_temp, Condition1_temp)
  ) %>%
  ungroup() %>%
  select(Experiment, Gene, Condition1, Condition2, Wasserstein, FDR, p95_W, W_Norm)

out_file <- file.path(base_path, "input/polyA_stats_curated.rds")
saveRDS(stats_data, out_file)

cat("Wrote curated stats to:", out_file, "\n")
cat("Rows:", nrow(stats_data), " Cols:", ncol(stats_data), "\n")
