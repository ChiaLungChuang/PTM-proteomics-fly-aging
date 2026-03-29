# ============================================================
# Project: Pan-PTM Proteomics — Drosophila Aging & Longevity
# Script:  PTM_fly_aging_visualization.R
# Purpose: Visualization of JUMPptm pan-PTM proteomics data
#          from skeletal muscle of long-lived (O1, O3) and
#          parental (B3) Drosophila strains at 1-week (young)
#          and 6-week (old) ages. Generates volcano plots,
#          PTM regulation bar plots, global pattern heatmaps,
#          and UpSet plots for 8 PTM types.
# Input:   O & B strains.xlsx (JUMPptm output, Supplementary
#          Table 4 from the published dataset)
# Output:  Volcano plot PDFs per comparison, PTM bar/heatmap
#          PDFs, UpSet plot PDFs per PTM type, summary CSVs
# Author:  Chia-Lung Chuang
# Updated: 2026-03
# Reference: Poudel S*, Chuang C-L*, Shrestha HK*, Demontis F.
#   Pan-PTM profiling identifies PTMs associated with
#   exceptional longevity and skeletal muscle function in
#   Drosophila. npj Aging. 2025;11:23.
#   https://doi.org/10.1038/s41514-025-00215-2
# Dependencies: BiocManager packages — install once:
#   BiocManager::install(c("ComplexHeatmap", "enrichplot"))
# ============================================================

library(tidyverse)
library(ggrepel)
library(readxl)
library(broom)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ComplexHeatmap)
library(UpSetR)
library(enrichplot)
library(pheatmap)
library(gridExtra)
library(RColorBrewer)

# NOTE: Set your working directory to the folder containing the input file
# setwd("path/to/your/data")

# ============================================================
# 1. Data Import and Reshaping
# ============================================================

data <- read_excel("O & B strains.xlsx")
write.csv(data, "O & B strains.csv", row.names = FALSE)
data_csv <- read.csv("O & B strains.csv")

# Parse the wide-format JUMPptm output into long format.
# Columns 7 onward encode comparisons as:
# <Time1><Gene1>.<Time2><Gene2>_<Statistic>
# e.g. "OnewkO1.OnewkB3_logFC" → Time1=One week, Gene1=O1,
#                                  Time2=One week, Gene2=B3, Statistic=logFC

data_long <- data_csv %>%
  pivot_longer(cols = 7:ncol(.), names_to = "comparison", values_to = "value")

# Helper: parse a single column name into its 5 components
process_comparison <- function(name) {
  parts      <- strsplit(name, "_")[[1]]
  comparison <- parts[1]
  statistic  <- if (length(parts) > 1) paste(parts[-1], collapse = "_") else NA
  
  timepoints <- strsplit(comparison, "\\.")[[1]]
  
  first_time <- substr(timepoints[1], 1, 5)
  first_time <- ifelse(first_time == "Onewk", "One week", "Six week")
  first_gene <- substr(timepoints[1], 6, nchar(timepoints[1]))
  
  second_time <- substr(timepoints[2], 1, 5)
  second_time <- ifelse(second_time == "Onewk", "One week", "Six week")
  second_gene <- substr(timepoints[2], 6, nchar(timepoints[2]))
  
  return(c(Time1 = first_time, Time2 = second_time,
           Gene1 = first_gene, Gene2 = second_gene,
           Statistic = statistic))
}

comparison_parts <- t(sapply(data_long$comparison, process_comparison))
data_long <- data_long %>%
  bind_cols(as.data.frame(comparison_parts)) %>%
  select(-comparison) %>%
  select(Identifier, Peptides, PTM.type, Protein.ID, Gene, Modsite,
         Time1, Time2, Gene1, Gene2, Statistic, value)

write.csv(data_long, "O & B strains organized.csv", row.names = FALSE)

# ============================================================
# 2. PTM Color Palettes
# ============================================================

# Bright and bold palette used throughout
ptm_colors_bold <- c(
  "oxidation"       = "#FF0000",
  "phosphorylation" = "#0000FF",
  "methylation_mono"= "#00FF00",
  "dihydroxy"       = "#FF00FF",
  "Deamidation"     = "#FFD700",
  "Carboxylation"   = "#00FFFF",
  "acetylation"     = "#FF8C00",
  "ubiquitination"  = "#8B008B"
)

ptm_labels <- c(
  "oxidation"        = "Oxidation",
  "phosphorylation"  = "Phosphorylation",
  "methylation_mono" = "Mono-methylation",
  "dihydroxy"        = "Dihydroxylation",
  "Deamidation"      = "Deamidation",
  "Carboxylation"    = "Carboxylation",
  "acetylation"      = "Acetylation",
  "ubiquitination"   = "Ubiquitination"
)

# ============================================================
# 3. Volcano Plot Functions
# ============================================================

# Significance threshold: adj.P.Val < 0.05 & |log2FC| > 1

# With gene + modification site labels on significant points
create_individual_volcano <- function(data, time1, time2, gene1, gene2) {
  plot_data <- data %>%
    filter(Time1 == time1, Time2 == time2, Gene1 == gene1, Gene2 == gene2) %>%
    select(Identifier, PTM.type, Gene, Protein.ID, Modsite, Statistic, value) %>%
    pivot_wider(names_from = Statistic, values_from = value) %>%
    mutate(
      significant    = adj.P.Val < 0.05 & abs(logFC) > 1,
      label_distance = abs(logFC) + (-log10(adj.P.Val))
    )
  
  ggplot(plot_data, aes(x = logFC, y = -log10(adj.P.Val))) +
    geom_point(aes(color = PTM.type, alpha = significant, shape = significant), size = 2.5) +
    geom_text_repel(
      data = subset(plot_data, significant),
      aes(label = paste(Gene, Modsite, sep = "\n")),
      max.overlaps = Inf, size = 2.5, force = 10,
      box.padding = 0.2, min.segment.length = 0.4,
      segment.color = "grey50", show.legend = FALSE
    ) +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "gray50") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray50") +
    annotate("text", x = -1.1, y = 0.2, label = "FC = -2", size = 3, angle = 90) +
    annotate("text", x =  1.1, y = 0.2, label = "FC = 2",  size = 3, angle = 90) +
    annotate("text", x = -3,   y = -log10(0.05) + 0.2, label = "p = 0.05", size = 3) +
    scale_color_manual(values = ptm_colors_bold, name = "PTM Type", labels = ptm_labels) +
    scale_alpha_manual(values = c("TRUE" = 1, "FALSE" = 0.4),
                       labels = c("TRUE" = "Significant", "FALSE" = "Not Significant")) +
    scale_shape_manual(values = c("TRUE" = 16, "FALSE" = 1),
                       labels = c("TRUE" = "Significant", "FALSE" = "Not Significant")) +
    theme_minimal() +
    theme(legend.position = "right", legend.box = "vertical",
          plot.title    = element_text(size = 12, face = "bold", hjust = 0.5),
          legend.title  = element_text(size = 10, face = "bold"),
          legend.text   = element_text(size = 9),
          axis.title    = element_text(size = 10),
          axis.text     = element_text(size = 9),
          panel.grid.major = element_line(color = "grey90"),
          panel.grid.minor = element_blank()) +
    labs(title   = paste(time1, gene1, "vs", time2, gene2),
         x       = "log2 Fold Change",
         y       = "-log10(adjusted P-value)",
         caption = "Significance: adj.P.Val < 0.05 & |log2FC| > 1") +
    coord_cartesian(clip = "off", expand = TRUE)
}

# Without labels (clean version)
create_individual_volcano_nolabel <- function(data, time1, time2, gene1, gene2) {
  plot_data <- data %>%
    filter(Time1 == time1, Time2 == time2, Gene1 == gene1, Gene2 == gene2) %>%
    select(Identifier, PTM.type, Gene, Protein.ID, Modsite, Statistic, value) %>%
    pivot_wider(names_from = Statistic, values_from = value) %>%
    mutate(significant = adj.P.Val < 0.05 & abs(logFC) > 1)
  
  ggplot(plot_data, aes(x = logFC, y = -log10(adj.P.Val))) +
    geom_point(aes(color = PTM.type, alpha = significant, shape = significant), size = 2.5) +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "gray50") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray50") +
    annotate("text", x = -1.1, y = 0.2, label = "FC = -2", size = 3, angle = 90) +
    annotate("text", x =  1.1, y = 0.2, label = "FC = 2",  size = 3, angle = 90) +
    annotate("text", x = -3,   y = -log10(0.05) + 0.2, label = "p = 0.05", size = 3) +
    scale_color_manual(values = ptm_colors_bold, name = "PTM Type", labels = ptm_labels) +
    scale_alpha_manual(values = c("TRUE" = 1, "FALSE" = 0.4),
                       labels = c("TRUE" = "Significant", "FALSE" = "Not Significant")) +
    scale_shape_manual(values = c("TRUE" = 16, "FALSE" = 1),
                       labels = c("TRUE" = "Significant", "FALSE" = "Not Significant")) +
    theme_minimal() +
    theme(legend.position = "right", legend.box = "vertical",
          plot.title   = element_text(size = 12, face = "bold", hjust = 0.5),
          legend.title = element_text(size = 10, face = "bold"),
          legend.text  = element_text(size = 9),
          axis.title   = element_text(size = 10),
          axis.text    = element_text(size = 9),
          panel.grid.major = element_line(color = "grey90"),
          panel.grid.minor = element_blank()) +
    labs(title   = paste(time1, gene1, "vs", time2, gene2),
         x       = "log2 Fold Change",
         y       = "-log10(adjusted P-value)",
         caption = "Significance: adj.P.Val < 0.05 & |log2FC| > 1") +
    coord_cartesian(clip = "off", expand = TRUE)
}

# ============================================================
# 4. Generate All Volcano Plots
# ============================================================

save_volcano <- function(plot_obj, filename) {
  print(plot_obj)
  ggsave(filename, plot = plot_obj, width = 5800, height = 5800, units = "px", dpi = 600)
}

# --- One week: O1 vs B3 ---
onewk_O1_B3        <- create_individual_volcano(        data_long, "One week", "One week", "O1", "B3")
onewk_O1_B3_nolabel <- create_individual_volcano_nolabel(data_long, "One week", "One week", "O1", "B3")
save_volcano(onewk_O1_B3,         "onewk_O1_B3.pdf")
save_volcano(onewk_O1_B3_nolabel, "onewk_O1_B3_nolabel.pdf")

# --- One week: O3 vs B3 ---
onewk_O3_B3        <- create_individual_volcano(        data_long, "One week", "One week", "O3", "B3")
onewk_O3_B3_nolabel <- create_individual_volcano_nolabel(data_long, "One week", "One week", "O3", "B3")
save_volcano(onewk_O3_B3,         "onewk_O3_B3.pdf")
save_volcano(onewk_O3_B3_nolabel, "onewk_O3_B3_nolabel.pdf")

# --- Six week: O1 vs B3 ---
sixwk_O1_B3        <- create_individual_volcano(        data_long, "Six week", "Six week", "O1", "B3")
sixwk_O1_B3_nolabel <- create_individual_volcano_nolabel(data_long, "Six week", "Six week", "O1", "B3")
save_volcano(sixwk_O1_B3,         "sixwk_O1_B3.pdf")
save_volcano(sixwk_O1_B3_nolabel, "sixwk_O1_B3_nolabel.pdf")

# --- Six week: O3 vs B3 ---
sixwk_O3_B3        <- create_individual_volcano(        data_long, "Six week", "Six week", "O3", "B3")
sixwk_O3_B3_nolabel <- create_individual_volcano_nolabel(data_long, "Six week", "Six week", "O3", "B3")
save_volcano(sixwk_O3_B3,         "sixwk_O3_B3.pdf")
save_volcano(sixwk_O3_B3_nolabel, "sixwk_O3_B3_nolabel.pdf")

# --- Aging within strain: Six week vs One week ---
B3_sixwk_onewk        <- create_individual_volcano(        data_long, "Six week", "One week", "B3", "B3")
B3_sixwk_onewk_nolabel <- create_individual_volcano_nolabel(data_long, "Six week", "One week", "B3", "B3")
save_volcano(B3_sixwk_onewk,         "B3_sixwk_onewk.pdf")
save_volcano(B3_sixwk_onewk_nolabel, "B3_sixwk_onewk_nolabel.pdf")

O1_sixwk_onewk        <- create_individual_volcano(        data_long, "Six week", "One week", "O1", "O1")
O1_sixwk_onewk_nolabel <- create_individual_volcano_nolabel(data_long, "Six week", "One week", "O1", "O1")
save_volcano(O1_sixwk_onewk,         "O1_sixwk_onewk.pdf")
save_volcano(O1_sixwk_onewk_nolabel, "O1_sixwk_onewk_nolabel.pdf")

O3_sixwk_onewk        <- create_individual_volcano(        data_long, "Six week", "One week", "O3", "O3")
O3_sixwk_onewk_nolabel <- create_individual_volcano_nolabel(data_long, "Six week", "One week", "O3", "O3")
save_volcano(O3_sixwk_onewk,         "O3_sixwk_onewk.pdf")
save_volcano(O3_sixwk_onewk_nolabel, "O3_sixwk_onewk_nolabel.pdf")

# ============================================================
# 5. PTM Regulation Bar Plots (per time point)
# ============================================================

# Bar plot showing count of up- vs. down-regulated PTM sites
# (|logFC| > 1) per PTM type, faceted by genotype comparison

analyze_ptm_regulation <- function(data, time_point) {
  data %>%
    filter(Time1 == time_point, Time2 == time_point, Statistic == "logFC") %>%
    group_by(PTM.type, Gene1, Gene2) %>%
    summarize(
      significant_up   = sum(value >  1),
      significant_down = sum(value < -1),
      .groups = "drop"
    ) %>%
    ggplot(aes(x = PTM.type)) +
    geom_bar(aes(y =  significant_up,   fill = "Up-regulated"),   stat = "identity") +
    geom_bar(aes(y = -significant_down, fill = "Down-regulated"), stat = "identity") +
    facet_grid(Gene1 ~ Gene2) +
    coord_flip() +
    scale_fill_manual(values = c("Up-regulated" = "#d95f02", "Down-regulated" = "#1b9e77")) +
    labs(title = paste("PTM Regulation Patterns at", time_point),
         y = "Number of Significantly Changed PTMs (|log2FC| > 1)",
         x = "PTM Type", fill = "Regulation") +
    theme_minimal() +
    theme(axis.text     = element_text(size = 10),
          axis.title    = element_text(size = 12),
          plot.title    = element_text(size = 14, face = "bold"),
          legend.title  = element_text(size = 10, face = "bold"),
          legend.text   = element_text(size = 9),
          strip.text    = element_text(size = 10, face = "bold"))
}

onewk_ptm_reg <- analyze_ptm_regulation(data_long, "One week")
sixwk_ptm_reg <- analyze_ptm_regulation(data_long, "Six week")

print(onewk_ptm_reg)
ggsave("onewk_ptm_reg.pdf",  plot = onewk_ptm_reg, width = 7200, height = 10800, units = "px", dpi = 600)
print(sixwk_ptm_reg)
ggsave("sixwk_ptm_reg.pdf",  plot = sixwk_ptm_reg, width = 7200, height = 10800, units = "px", dpi = 600)

# ============================================================
# 6. Aging-Specific PTM Change Patterns
# ============================================================

# Shows how many PTM sites changed with aging and whether
# the change was B3-specific, O1/O3-specific, or shared

analyze_aging_ptm_changes <- function(data) {
  aging_changes <- data %>%
    filter(Time1 != Time2, Statistic == "logFC") %>%
    group_by(PTM.type, Identifier, Gene1) %>%
    summarize(is_significant = abs(value) > 1, .groups = "drop") %>%
    pivot_wider(names_from = Gene1, values_from = is_significant, names_prefix = "changed_in_")
  
  plot_data <- aging_changes %>%
    group_by(PTM.type) %>%
    summarize(
      B3_specific  = sum( changed_in_B3 & !changed_in_O1 & !changed_in_O3),
      O1_specific  = sum(!changed_in_B3 &  changed_in_O1 & !changed_in_O3),
      O3_specific  = sum(!changed_in_B3 & !changed_in_O1 &  changed_in_O3),
      O1_O3_common = sum(!changed_in_B3 &  changed_in_O1 &  changed_in_O3),
      shared_all   = sum( changed_in_B3 &  changed_in_O1 &  changed_in_O3),
      .groups = "drop"
    ) %>%
    pivot_longer(cols = c(B3_specific, O1_specific, O3_specific, O1_O3_common, shared_all),
                 names_to = "pattern", values_to = "count") %>%
    mutate(pattern = factor(pattern,
                            levels = c("B3_specific","O1_specific","O3_specific","O1_O3_common","shared_all"),
                            labels = c("B3 Specific","O1 Specific","O3 Specific","O1-O3 Common","Shared All")),
           PTM.type = recode(PTM.type,
                             oxidation = "Oxidation", phosphorylation = "Phosphorylation",
                             methylation_mono = "Methylation_mono", dihydroxy = "Dihydroxy",
                             Deamidation = "Deamidation", Carboxylation = "Carboxylation",
                             acetylation = "Acetylation", ubiquitination = "Ubiquitination"))
  
  ggplot(plot_data, aes(x = PTM.type, y = count, fill = pattern)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.9),
             color = "black", linewidth = 0.2) +
    coord_flip() +
    labs(title    = "Aging-Related PTM Changes by Genotype",
         subtitle = "Distribution of aging-associated PTM changes across B3, O1, O3",
         x = "PTM Type", y = "Number of Modified Sites", fill = "Change Pattern") +
    theme_minimal() +
    theme(axis.text.y   = element_text(size = 12, face = "bold"),
          plot.title    = element_text(size = 14, face = "bold"),
          plot.subtitle = element_text(size = 12),
          legend.title  = element_text(face = "bold")) +
    scale_fill_brewer(palette = "Set2")
}

aging_changes <- analyze_aging_ptm_changes(data_long)
print(aging_changes)
ggsave("aging_changes.pdf", plot = aging_changes, width = 7200, height = 10800, units = "px", dpi = 600)

# ============================================================
# 7. Global PTM Pattern Heatmap (% sites changed, by adj.P.Val)
# ============================================================

# Heatmap showing percentage of significantly changed PTM sites
# (adj.P.Val < 0.05) across all 7 genotype/age comparisons

analyze_global_ptm_heatmap <- function(data_long) {
  # Standardized comparison labels and ordering
  comp_levels <- c("One week O1 vs B3", "One week O3 vs B3",
                   "Six week O1 vs B3", "Six week O3 vs B3",
                   "One week vs Six week B3", "One week vs Six week O1",
                   "One week vs Six week O3")
  
  processed <- data_long %>%
    filter(Statistic == "adj.P.Val") %>%
    group_by(PTM.type, Time1, Time2, Gene1, Gene2) %>%
    summarize(
      n_significant  = n_distinct(Identifier[value < 0.05]),
      total_proteins = n_distinct(Identifier),
      percent_changed = n_significant / total_proteins * 100,
      .groups = "drop"
    ) %>%
    mutate(
      comparison = case_when(
        Time1 == Time2 ~ paste(Time1, paste(Gene1, "vs", Gene2)),
        TRUE           ~ paste(Time2, "vs", Time1, Gene1)
      ),
      comparison = factor(comparison, levels = comp_levels),
      PTM.type   = recode(PTM.type,
                          dihydroxy = "Dihydroxy", oxidation = "Oxidation",
                          phosphorylation = "Phosphorylation", ubiquitination = "Ubiquitination",
                          acetylation = "Acetylation", methylation_mono = "Methylation_mono")
    )
  
  ggplot(processed, aes(x = comparison, y = PTM.type)) +
    geom_tile(aes(fill = percent_changed), color = "white", width = 0.9) +
    geom_text(aes(label = sprintf("%.1f%%\n(n=%d)", percent_changed, n_significant)),
              size = 3, color = "black") +
    scale_fill_gradient(low = "#4575B4", high = "#D73027") +
    theme_minimal() +
    theme(axis.text.x  = element_text(angle = 45, hjust = 1, margin = margin(t = 10)),
          axis.text    = element_text(size = 10),
          panel.grid   = element_blank(),
          panel.border = element_rect(fill = NA, color = "grey80"),
          legend.title = element_text(size = 10, face = "bold"),
          plot.margin  = margin(20, 20, 20, 20)) +
    labs(title    = "Global PTM Regulation Patterns",
         subtitle = "PTM changes across all genotype/age comparisons",
         caption  = "Significant changes: adj.P.Val < 0.05\nn = number of significantly changed sites",
         x = "Comparison", y = "PTM Type", fill = "% Sites Changed")
}

global_patterns_heatmap_sig <- analyze_global_ptm_heatmap(data_long)
print(global_patterns_heatmap_sig)
ggsave("global_patterns_heatmap_sig.pdf", plot = global_patterns_heatmap_sig,
       width = 5400, height = 5400, units = "px", dpi = 600)

# ============================================================
# 8. Global PTM Pattern Bar Plot (% sites, colored by mean FC)
# ============================================================

analyze_global_ptm_bars <- function(data_long) {
  comp_levels <- c("One week O1 vs B3", "One week O3 vs B3",
                   "Six week O1 vs B3", "Six week O3 vs B3",
                   "One week vs Six week B3", "One week vs Six week O1",
                   "One week vs Six week O3")
  
  processed <- data_long %>%
    group_by(PTM.type, Time1, Time2, Gene1, Gene2, Identifier) %>%
    summarize(
      sig     = any(Statistic == "adj.P.Val" & value < 0.05),
      mean_fc = mean(value[Statistic == "logFC"], na.rm = TRUE),
      .groups = "drop"
    ) %>%
    group_by(PTM.type, Time1, Time2, Gene1, Gene2) %>%
    summarize(
      n_significant    = sum(sig),
      total_proteins   = n(),
      percent_changed  = n_significant / total_proteins * 100,
      mean_fold_change = mean(mean_fc[sig], na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      comparison = case_when(
        Time1 == Time2 ~ paste(Time1, paste(Gene1, "vs", Gene2)),
        TRUE           ~ paste(Time2, "vs", Time1, Gene1)
      ),
      comparison = factor(comparison, levels = comp_levels),
      PTM.type   = recode(PTM.type,
                          dihydroxy = "Dihydroxy", oxidation = "Oxidation",
                          phosphorylation = "Phosphorylation", ubiquitination = "Ubiquitination",
                          acetylation = "Acetylation", methylation_mono = "Methylation_mono")
    )
  
  fc_limit <- max(abs(processed$mean_fold_change), na.rm = TRUE)
  
  ggplot(processed, aes(x = comparison, y = percent_changed)) +
    geom_bar(aes(fill = mean_fold_change), stat = "identity", width = 0.7,
             color = "black", linewidth = 0.2) +
    geom_text(aes(label = sprintf("%.1f%%", percent_changed)),
              vjust = -0.5, size = 3) +
    geom_text(aes(label = sprintf("(n=%d)", n_significant)),
              vjust = -2, size = 3) +
    facet_wrap(~PTM.type, scales = "free_y", ncol = 2) +
    scale_fill_gradientn(colors = c("#4575B4", "#D73027"),
                         limits = c(-fc_limit, fc_limit),
                         oob = scales::squish) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.4)), limits = c(0, NA)) +
    theme_minimal() +
    theme(axis.text.x   = element_text(angle = 45, hjust = 1, vjust = 1),
          axis.text     = element_text(size = 8),
          strip.text    = element_text(size = 12, face = "bold"),
          strip.background = element_rect(fill = "grey95"),
          panel.border  = element_rect(fill = NA, color = "grey80"),
          plot.margin   = margin(30, 30, 30, 30)) +
    labs(title    = "Global PTM Regulation Patterns",
         subtitle = "Comparing PTM changes between and within genotypes",
         caption  = "Significant: adj.P.Val < 0.05\nn = number of significantly changed sites",
         x = "Comparison", y = "Percentage of Sites Changed (%)",
         fill = "Mean Fold Change")
}

global_patterns_bars <- analyze_global_ptm_bars(data_long)
print(global_patterns_bars)
ggsave("global_patterns_bars_p_sig.pdf", plot = global_patterns_bars,
       width = 5400, height = 7200, units = "px", dpi = 600)

# ============================================================
# 9. Per-PTM-Type Heatmaps (log2FC across genes)
# ============================================================

create_ptm_heatmap <- function(data, selected_ptm_type) {
  data %>%
    filter(PTM.type == selected_ptm_type, Statistic == "logFC") %>%
    ggplot(aes(x = paste(Time1, Gene1, "vs", Time2, Gene2), y = Gene)) +
    geom_tile(aes(fill = value)) +
    scale_fill_gradient2(low = "#4575B4", mid = "white", high = "#D73027") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          axis.text.y = element_text(size = 8)) +
    labs(title = paste("Heatmap of", selected_ptm_type, "Sites"),
         x = "Comparison", y = "Gene", fill = "log2FC")
}

for (ptm in c("oxidation","phosphorylation","methylation_mono","dihydroxy",
              "Deamidation","Carboxylation","acetylation","ubiquitination")) {
  p <- create_ptm_heatmap(data_long, ptm)
  print(p)
}

# ============================================================
# 10. Global UpSet Plot (all PTM types, adj.P.Val < 0.05)
# ============================================================

create_global_upset_plot <- function(data) {
  upset_data <- data %>%
    filter(Statistic == "adj.P.Val", value < 0.05) %>%
    mutate(comparison = case_when(
      Time1 == Time2                              ~ paste(Time1, Gene1, "vs", Gene2),
      Time1 == "Six week" & Time2 == "One week"  ~ paste("Six week vs One week", Gene2)
    )) %>%
    select(Identifier, comparison) %>%
    distinct() %>%
    mutate(present = 1) %>%
    pivot_wider(id_cols = Identifier, names_from = comparison,
                values_from = present, values_fill = 0)
  
  matrix_data <- as.data.frame(upset_data[, -1])
  rownames(matrix_data) <- upset_data$Identifier
  
  upset(matrix_data,
        nsets = ncol(matrix_data), nintersects = 60, order.by = "freq",
        main.bar.color = "darkred", sets.bar.color = "darkblue",
        text.scale = 1.2, point.size = 2.5, line.size = 1,
        mb.ratio = c(0.6, 0.4), show.numbers = "yes",
        mainbar.y.label = "Number of Protein-PTM combinations",
        sets.x.label    = "Proteins per Comparison")
}

pdf("global_upset_plot_sig.pdf", width = 12, height = 8, onefile = FALSE)
create_global_upset_plot(data_long)
dev.off()

# ============================================================
# 11. Per-PTM-Type UpSet Plots (ComplexHeatmap)
# ============================================================

# Filter threshold: |logFC| > 1
# Generates one PDF per PTM type showing intersection of
# comparisons, with bar annotations and summary CSV

create_ptm_distribution_upset <- function(data, filter_type = "logFC") {
  if (filter_type == "logFC") {
    upset_data <- data %>%
      filter(Statistic == "logFC", abs(value) > 1) %>%
      mutate(comparison = case_when(
        Time1 == Time2 ~ paste(Time1, Gene1, "vs", Gene2),
        TRUE           ~ paste(Time1, "vs", Time2, Gene1)
      ))
    threshold_label <- "proteins with |logFC| > 1"
    file_prefix_fn  <- function(ptm) paste0("upset_plot_", tolower(gsub(" ", "_", ptm)))
  } else {
    upset_data <- data %>%
      filter(Statistic == "adj.P.Val", value < 0.05) %>%
      mutate(comparison = case_when(
        Time1 == Time2 ~ paste(Time1, Gene1, "vs", Gene2),
        TRUE           ~ paste(Time1, "vs", Time2, Gene1)
      ))
    threshold_label <- "proteins with adj.P.Val < 0.05"
    file_prefix_fn  <- function(ptm) paste0("upset_plot_sig_", tolower(gsub(" ", "_", ptm)))
  }
  
  ptm_types <- if (filter_type == "logFC") {
    setdiff(unique(upset_data$PTM.type), "ubiquitination")
  } else {
    unique(upset_data$PTM.type)
  }
  
  for (ptm in ptm_types) {
    presence_df <- upset_data %>%
      filter(PTM.type == ptm) %>%
      select(Identifier, comparison) %>%
      distinct() %>%
      mutate(present = 1) %>%
      spread(comparison, present, fill = 0)
    
    if (nrow(presence_df) < 2) next
    
    matrix_data    <- as.matrix(presence_df[, -1])
    rownames(matrix_data) <- presence_df$Identifier
    
    gene_modsite_map <- upset_data %>%
      filter(PTM.type == ptm) %>%
      select(Identifier, Gene, Modsite) %>%
      distinct() %>%
      mutate(Gene_Modsite = paste(Gene, Modsite, sep = "_"))
    
    m          <- make_comb_mat(matrix_data, mode = "distinct")
    set_sz     <- set_size(m)
    comb_sz    <- comb_size(m)
    plot_order <- order(comb_degree(m), -comb_sz)
    
    # Identify which identifiers fall in each intersection pattern
    get_gene_modsites <- function(pattern, mat, gm_map) {
      pb   <- as.logical(as.numeric(strsplit(pattern, "")[[1]]))
      rows <- which(rowSums(mat[, pb,  drop = FALSE]) == sum(pb) &
                      rowSums(mat[, !pb, drop = FALSE]) == 0)
      gm_map %>% filter(Identifier %in% rownames(mat)[rows]) %>% pull(Gene_Modsite)
    }
    
    ordered_patterns    <- names(comb_sz)[plot_order]
    ordered_counts      <- comb_sz[plot_order]
    ordered_gene_modsites <- lapply(ordered_patterns, get_gene_modsites,
                                    mat = matrix_data, gm_map = gene_modsite_map)
    
    top_anno <- HeatmapAnnotation(
      "Protein-PTM" = anno_barplot(
        ordered_counts, ylim = c(0, max(ordered_counts) * 1.2),
        border = FALSE, gp = gpar(fill = ptm_colors_bold[ptm]),
        height = unit(5, "cm"), axis = TRUE,
        labels = ordered_counts, labels_rot = 0,
        labels_offset = unit(0.5, "cm")
      ),
      annotation_name_side = "left", annotation_name_rot = 90,
      height = unit(6, "cm")
    )
    
    left_anno <- rowAnnotation(
      set_name = anno_text(colnames(matrix_data), location = 0.1, just = "left",
                           width = unit(4, "cm"), gp = gpar(fontsize = 8))
    )
    
    right_anno <- rowAnnotation(
      "Protein-PTM" = anno_barplot(set_sz, border = FALSE,
                                   gp = gpar(fill = ptm_colors_bold[ptm]),
                                   width = unit(3, "cm"), axis = TRUE)
    )
    
    pdf(paste0(file_prefix_fn(ptm), ".pdf"), width = 10, height = 8)
    ht <- UpSet(m,
                set_order    = order(set_sz),
                comb_order   = plot_order,
                top_annotation  = top_anno,
                left_annotation = left_anno,
                right_annotation = right_anno,
                show_row_names = FALSE,
                column_title = paste0(
                  toupper(substr(ptm, 1, 1)), substr(ptm, 2, nchar(ptm)),
                  " Modifications\n(", nrow(presence_df), " ", threshold_label, ")"
                ),
                height = unit(5, "cm"), width = unit(15, "cm"),
                pt_size = unit(5, "mm")
    )
    draw(ht)
    dev.off()
    
    # Save intersection summary CSV
    summary_df <- data.frame(
      Column      = seq_along(ordered_gene_modsites),
      Count       = ordered_counts,
      Gene_Modsites = sapply(ordered_gene_modsites, paste, collapse = "; ")
    )
    write.csv(summary_df, paste0(file_prefix_fn(ptm), "_summary.csv"), row.names = FALSE)
  }
}

# |logFC| > 1 threshold
create_ptm_distribution_upset(data_long, filter_type = "logFC")

# adj.P.Val < 0.05 threshold
create_ptm_distribution_upset(data_long, filter_type = "adj.P.Val")

# ============================================================
# 12. Summary Statistics and CSV Exports
# ============================================================

# Overall comparison summary
comparison_summary <- data_long %>%
  group_by(Time1, Time2, Gene1, Gene2) %>%
  summarize(
    n_significant    = sum(Statistic == "adj.P.Val" & value < 0.05),
    max_fold_change  = max(abs(value[Statistic == "logFC"])),
    .groups = "drop"
  )
print(comparison_summary)

# Visualization summary CSVs (for heatmap and bar plot)
basic_summary <- data_long %>%
  filter(Statistic == "adj.P.Val") %>%
  group_by(PTM.type, Time1, Time2, Gene1, Gene2) %>%
  summarize(
    n_significant   = n_distinct(Identifier[value < 0.05]),
    total_proteins  = n_distinct(Identifier),
    percent_changed = n_significant / total_proteins * 100,
    sig_proteins    = paste(sort(unique(Gene[value < 0.05])), collapse = "; "),
    .groups = "drop"
  )
write.csv(basic_summary, "heatmap_and_basic_barplot_summary.csv", row.names = FALSE)

# UpSet summary across all PTM types
ptm_upset_summary <- data_long %>%
  filter(Statistic == "adj.P.Val", value < 0.05) %>%
  mutate(comparison = case_when(
    Time1 == Time2 ~ paste(Time1, Gene1, "vs", Gene2),
    TRUE           ~ paste(Time1, "vs", Time2, Gene1)
  )) %>%
  group_by(PTM.type) %>%
  summarize(
    total_proteins = n_distinct(Identifier),
    unique_genes   = n_distinct(Gene),
    comparisons    = paste(sort(unique(comparison)), collapse = "; "),
    protein_sites  = paste(sort(unique(paste(Gene, Modsite, sep = "_"))), collapse = "; "),
    .groups = "drop"
  )
write.csv(ptm_upset_summary, "ptm_upset_summary.csv", row.names = FALSE)

# Global UpSet intersections
global_upset_intersection <- data_long %>%
  filter(Statistic == "adj.P.Val", value < 0.05) %>%
  mutate(comparison = case_when(
    Time1 == Time2                             ~ paste(Time1, Gene1, "vs", Gene2),
    Time1 == "Six week" & Time2 == "One week" ~ paste("Six week vs One week", Gene2)
  )) %>%
  group_by(Identifier, Gene) %>%
  summarize(
    comparisons      = paste(sort(unique(comparison)), collapse = "; "),
    num_comparisons  = n_distinct(comparison),
    ptm_types        = paste(sort(unique(PTM.type)), collapse = "; "),
    .groups = "drop"
  ) %>%
  arrange(desc(num_comparisons))

write.csv(global_upset_intersection, "global_upset_intersection_data.csv", row.names = FALSE)