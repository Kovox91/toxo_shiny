# Load libs
library(tidyverse)
theme_diss <- function() {
  theme_bw(base_family = "Liberation Sans") +
    theme(
      text = element_text(family = "Liberation Sans", colour = "black"),
      
      axis.title = element_text(size = 8),
      axis.text  = element_text(size = 6),
      
      legend.title = element_text(size = 8, face = "bold"),
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.text = element_text(size = 6),
      legend.box.spacing = unit(2, "mm"),
      legend.key.width = unit(4, "mm"),
      
      plot.title = element_blank(),   # enforce no titles
      plot.margin = margin(2, 2, 2, 2, "mm"),
      
      panel.grid = element_blank(),
      
      strip.text = element_text(size = 8, margin = margin(b = 2)),
      strip.background = element_rect(
        fill = "grey95",
        colour = NA
      )
    )
}
# load inputs
helicase <- readRDS("input/mrhel1_polyA.RDS") %>%  select(everything(), -n) %>%  mutate(experiment = "Helicase")
hpr <- readRDS("input/mrsic_polyA.RDS") %>% select(everything(), "con" = condition) %>%  mutate(experiment = "HPR")


# load references
names_raw <- read_csv("reference/naming_help.csv")
polyA_shikha_raw <- read.delim2("reference/shikha_polyA.tsv")
polyA_wang_raw <- read_csv("reference/wang_polyA_lengths.csv", skip = 2)
better_names_raw <- read.delim2("reference/cryoem_to_gtf_mapping.tsv") 

names <- better_names_raw %>% 
  select("name" = cryoem_fragment,
         "Geneid" = gtf_rRNA_id,
         "ID_long" = gtf_rRNA_id,
         "Source" = cryoem_source) %>% 
  mutate(name = base::gsub("-", "", name)) %>% 
  filter(!Geneid %in% c("No_GTF_Match", "SSU-X"))

polyA_shikha <- polyA_shikha_raw %>% 
  select("name" = Nomenclature,
         type,
         "polyA" = polya_length,
         "n_mean" = polya_pct,
         "polyA_struct" = polyA
         ) %>% 
  mutate(experiment = "Shikha et al.",
         name = base::gsub("-", "", name),
         n_mean = as.numeric(n_mean)) %>% 
  left_join(filter(names, Source == "Shikha2025") , by = "name",
            relationship = "many-to-many") %>% 
  select("Geneid" = ID_long, "con" = type, polyA, n_mean, experiment, polyA_struct, name)

polyA_wang <- polyA_wang_raw %>% 
  select("name" = Nomenclature,
         seq = `rRNA Sequence (5'-3')`) %>% 
  mutate(name = base::gsub("-", "", name),
         polya_length = nchar(sub(".*?(A+)$", "\\1", seq)),
         polya_length = ifelse(grepl("A+$", seq), polya_length, 0)) %>% 
  select(-seq) %>% 
  left_join(filter(names, Source == "Wang2024"), by = "name",
            relationship = "many-to-many") %>%
  mutate(ID_long = ifelse(is.na(ID_long), name, ID_long),
         ID_long = ifelse(ID_long == "RNA37/38_F", "RNA37-38_F", Geneid)) 


# Plot
sascha <- rbind(helicase, hpr) %>% 
  mutate(n_mean = 100* n_mean) %>% 
  rbind(polyA_shikha) %>% 
  mutate(
    con_unified = case_when(
      # controls
      con == "Untreated" ~ "Untreated",
      con == "PolyA-KO"  ~ "PolyA-KO",
      
      # parental (mR*)
      con %in% c("mRHel-FLAG -ATc", "mRSiC-FLAG -ATc") ~ "parental 0d",
      con %in% c("mRHel-FLAG 3d ATc")                  ~ "parental 3d",
      
      # promoter replaced (rmR*) — separate days
      con %in% c("rmRHel-FLAG -ATc", "rmRHel-FLAG 0d ATc",
                 "rmRSiC-FLAG -ATc")                   ~ "PR 0d",
      
      con %in% c("rmRHel-FLAG 1d ATc", "rmRSiC-FLAG 1d ATc") ~ "PR 1d",
      con %in% c("rmRHel-FLAG 2d ATc", "rmRSiC-FLAG 2d ATc") ~ "PR 2d",
      con %in% c("rmRHel-FLAG 3d ATc", "rmRSiC-FLAG 3d ATc") ~ "PR 3d",
      
      TRUE ~ NA_character_
    ),
    Geneid = ifelse(is.na(Geneid), name, Geneid)
  )

for (gene in unique(sascha$Geneid)){
  plot <- 
    ggplot(filter(sascha, Geneid == gene),
           aes(
             x = polyA,
             y = n_mean,
             color = con_unified
           )) +
    geom_line() +
    theme_bw() +
    geom_vline(data = filter(better_names_raw, gtf_rRNA_id == gene & cryoem_source == "Shikha2025"),
               aes(xintercept = cryoem_bases_trimmed-0.05),
               color = "red",
               linetype = "dashed") +
    geom_vline(data = filter(better_names_raw, gtf_rRNA_id == gene & cryoem_source == "Wang2024"),
               aes(xintercept = cryoem_bases_trimmed+0.05),
               color = "black",
               linetype = "dashed") +
    facet_wrap(~experiment,
               ncol = 1,
               scales = ifelse(gene == "RNA34_F", "free_y", "fixed")) +
    scale_x_continuous(
      breaks = function(x) floor(x[1]):ceiling(x[2])
    )  +
    theme(
      panel.grid.major = element_line(color = "grey90"),
      panel.grid.minor = element_blank()
    ) +
    labs(title = gene)
  
  ggsave(filename = paste0("out/plots/", gene, ".png"),
         plot = plot,
         height = 24,
         width = 16,
         units = "cm")
}

diss_genes <- c("RNA34_F", "RNA42", "RNA33_1", "LSUC")
plot <- 
  ggplot(filter(sascha, Geneid %in% diss_genes),
         aes(
           x = polyA,
           y = n_mean,
           color = con_unified
         )) +
  geom_line(linewidth = 0.25) +
  theme_diss() +
  geom_vline(data = filter(better_names_raw, gtf_rRNA_id %in% diss_genes & cryoem_source == "Shikha2025"),
             aes(xintercept = cryoem_bases_trimmed-0.05),
             color = "red",
             linetype = "dashed",
             linewidth = 0.22) +
  geom_vline(data = filter(better_names_raw, gtf_rRNA_id %in% diss_genes & cryoem_source == "Wang2024"),
             aes(xintercept = cryoem_bases_trimmed+0.05),
             color = "black",
             linetype = "dashed",
             linewidth = 0.22) +
  facet_grid(experiment~Geneid,
             scales = "free") +
  scale_x_continuous(
    breaks = function(x) 0:ceiling(x[2]),  # grid line every 1
    labels = function(x) ifelse(x %% 5 == 0, x, "")  # label every 2nd
  ) +
  theme(
    panel.grid.major = element_line(color = "grey90", linewidth = 0.2),
    panel.grid.minor = element_blank()
  ) +
  labs(x = "Poly(A) tail length (nt)",
       y = "Relative Read Frequency (%)")

plot
ggsave(filename = "~/diss_plots/polyA_plots.svg",
       plot = plot,
       height = 130,
       width = 155,
       units = "mm")
  