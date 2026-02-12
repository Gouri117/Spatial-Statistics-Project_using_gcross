
# ------------------------- Load Required Libraries -------------------------
library(dplyr)
library(readr)
library(ggplot2)
library(rstatix)
library(ggpubr)

# ------------------------- Load Metadata -------------------------
meta_file <- "/nfs/turbo/umms-ukarvind/shared_data/gcross_data/example/results/merged_metadata_clusters_k20.csv"
meta_df <- read_csv(meta_file)

# ------------------------- Count CN Frequencies Per Sample -------------------------
cn_counts <- meta_df %>%
  group_by(sample, neighborhood_cluster) %>%
  summarise(cn_count = n(), .groups = "drop")

sample_totals <- meta_df %>%
  group_by(sample) %>%
  summarise(total_cells = n(), .groups = "drop")

cn_freq <- cn_counts %>%
  left_join(sample_totals, by = "sample") %>%
  mutate(frequency = cn_count / total_cells)

# ------------------------- Attach Outcome Labels -------------------------
cn_freq <- cn_freq %>%
  left_join(meta_df %>% select(sample, primary_outcome) %>% distinct(), by = "sample") %>%
  mutate(
    outcome_label = case_when(
      primary_outcome == 0 ~ "Dead/recurred in 60 months",
      primary_outcome == 1 ~ "Alive after 60 months"
    ),
    outcome_label = factor(outcome_label, levels = c("Dead/recurred in 60 months", "Alive after 60 months")),
    neighborhood_cluster = as.character(neighborhood_cluster)
  )

# ------------------------- Wilcoxon Test with rstatix -------------------------
stat.test <- cn_freq %>%
  group_by(neighborhood_cluster) %>%
  wilcox_test(frequency ~ outcome_label, detailed = TRUE) %>%
  adjust_pvalue(method = "fdr") %>%
  add_significance("p.adj") %>%
  add_xy_position(x = "outcome_label")

# Save full stats table
stat.test_sub <- stat.test[, !names(stat.test) %in% c("y.position", "x", "groups", "xmin", "xmax")]
write.csv(stat.test_sub, "/nfs/turbo/umms-ukarvind/shared_data/gcross_data/example/results/cn_wilcox_stats_rstatix.csv", row.names = FALSE)

# ------------------------- Filter Significant CNs -------------------------
signif_cns <- stat.test %>% filter(p.adj.signif != "ns") %>% pull(neighborhood_cluster)
signif_df <- cn_freq %>% filter(neighborhood_cluster %in% signif_cns)
stat.test.sig <- stat.test %>% filter(neighborhood_cluster %in% signif_cns)

# ------------------------- Create Output Directory -------------------------
output_dir <- "/nfs/turbo/umms-ukarvind/shared_data/gcross_data/example/results/individual_cn_boxplots_rstatix/"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# ------------------------- Generate & Save Individual Plots -------------------------
for (cn in unique(signif_df$neighborhood_cluster)) {
  df_cn <- signif_df %>% filter(neighborhood_cluster == cn)
  stat_cn <- stat.test.sig %>% filter(neighborhood_cluster == cn)
  
  p <- ggboxplot(df_cn, x = "outcome_label", y = "frequency", fill = "outcome_label", palette = c("red", "blue")) +
    geom_jitter(width = 0.15, size = 1.5, alpha = 0.5) +
    stat_pvalue_manual(stat_cn, label = "p.adj.signif", size = 4) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
    labs(
      title = paste("CN", cn, "- Frequency by Outcome"),
      x = "Clinical Outcome",
      y = "CN Frequency",
      fill = "Outcome"
    ) +
    theme_bw(base_size = 13) +
    theme(axis.text.x = element_text(angle = 30, hjust = 1))
  
  # Display
  print(p)
  
  # Save
  ggsave(filename = paste0(output_dir, "CN_", cn, "_boxplot.png"), plot = p, width = 8, height = 5, dpi = 300)
}

cat(length(signif_cns), "significant CNs saved in:", output_dir, "\n")
