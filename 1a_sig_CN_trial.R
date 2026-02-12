library(tidyverse)

# Load your CSV file
gcross_df <- read_csv("/nfs/turbo/umms-ukarvind/shared_data/mfIHC_collab_data/stanford_crc/results/g_cross_cleaned_source_target.csv")

# Inspect the outcome column (e.g., "outcome_label") and CN column (e.g., "neighborhood_cluster")
glimpse(gcross_df)

gcross_df <- gcross_df %>%
  rename(sample = image)

charville_df <- read_csv("/nfs/turbo/umms-ukarvind/shared_data/mfIHC_collab_data/stanford_crc/charville_labels_samples_with_outcome.csv")

# ------------------ Filter and Merge ------------------ #
# Keep only relevant columns from charville
charville_subset <- charville_df %>% select(sample, primary_outcome)

# Merge outcome into gcross
gcross_df <- gcross_df %>%
  filter(sample %in% charville_subset$sample) %>%
  left_join(charville_subset, by = "sample")

# ------------------ Check result ------------------ #
glimpse(gcross_df)

# Count number of cells per CN per sample
cn_counts <- gcross_df %>%
  group_by(sample, neighborhood_cluster) %>%
  summarise(cn_count = n(), .groups = "drop")

# Total cell count per sample
sample_totals <- gcross_df %>%
  group_by(sample) %>%
  summarise(total_cells = n(), .groups = "drop")

# Merge and calculate frequency
cn_freq <- left_join(cn_counts, sample_totals, by = "sample") %>%
  mutate(cn_frequency = cn_count / total_cells)
library(rstatix)
library(dplyr)

# Ensure primary_outcome is in cn_freq
# Perform Wilcoxon test across CNs using rstatix
stat_results <- cn_freq %>%
  group_by(neighborhood_cluster) %>%
  wilcox_test(cn_frequency ~ primary_outcome, detailed = TRUE) %>%
  adjust_pvalue(method = "fdr") %>%
  add_significance("p.adj") %>%
  rename(CN = neighborhood_cluster) %>%
  arrange(p.adj)

print(stat.test)

sig_results <- stat.test %>%
  filter(p < 0.05)

library(ggpubr)

# Plot each significant CN
for (cn in sig_results$CN) {
  plot_df <- cn_freq %>% filter(neighborhood_cluster == cn)
  
  p <- ggbarplot(
    plot_df, 
    x = "primary_outcome", 
    y = "cn_frequency",
    add = "mean_se", 
    fill = "primary_outcome",
    palette = "jco",
    title = paste("CN", cn, "- Frequency by Primary Outcome"),
    ylab = "Frequency",
    xlab = "Primary Outcome"
  ) +
    stat_pvalue_manual(
      stat_results %>% filter(CN == cn),
      label = "p.adj.signif",
      y.position = max(plot_df$cn_frequency) * 1.05
    )
  
  print(p)
}
