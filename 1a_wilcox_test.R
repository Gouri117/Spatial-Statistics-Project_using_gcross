# ------------------------- Load Required Libraries -------------------------
library(dplyr)
library(readr)
library(ggplot2)

# ------------------------- Set Output Directory -------------------------
output_dir <- "/nfs/turbo/umms-ukarvind/shared_data/gcross_data/example/results/cn_boxplots/"
dir.create(output_dir, showWarnings = FALSE)

# ------------------------- Load Metadata -------------------------
meta_file <- "/nfs/turbo/umms-ukarvind/shared_data/gcross_data/example/results/merged_metadata_clusters_k20.csv"
meta_df <- read_csv(meta_file)

# ------------------------- Get All Unique CNs and Samples -------------------------
all_CNs <- unique(meta_df$neighborhood_cluster)
all_samples <- unique(meta_df$sample)

# ------------------------- Create All Sample × CN Combinations -------------------------
sample_cn_grid <- expand.grid(
  sample = all_samples,
  neighborhood_cluster = all_CNs,
  stringsAsFactors = FALSE
)

# ------------------------- Count CN Frequencies Per Sample -------------------------
cn_counts <- meta_df %>%
  group_by(sample, neighborhood_cluster) %>%
  summarise(cn_count = n(), .groups = "drop")

sample_totals <- meta_df %>%
  group_by(sample) %>%
  summarise(total_cells = n(), .groups = "drop")

# ------------------------- Compute CN frequency for all sample × CNs -------------------------
cn_freq <- sample_cn_grid %>%
  left_join(cn_counts, by = c("sample", "neighborhood_cluster")) %>%
  left_join(sample_totals, by = "sample") %>%
  mutate(
    cn_count = ifelse(is.na(cn_count), 0, cn_count),
    frequency = cn_count / total_cells
  ) %>%
  left_join(meta_df %>% select(sample, primary_outcome) %>% distinct(), by = "sample")

# ------------------------- Run Wilcoxon Rank-Sum Test Per CN -------------------------
wilcox_results <- cn_freq %>%
  group_by(neighborhood_cluster) %>%
  summarise(
    p_value = wilcox.test(frequency ~ primary_outcome)$p.value,
    .groups = "drop"
  ) %>%
  mutate(
    p_adj = p.adjust(p_value, method = "BH"),
    significant = p_adj < 0.05
  )

# Save Wilcoxon test results
write_csv(wilcox_results, "/nfs/turbo/umms-ukarvind/shared_data/gcross_data/example/results/cn_wilcox_fdr.csv")

# ------------------------- Filter Significant CNs -------------------------
sig_CNs <- wilcox_results %>% filter(significant == TRUE)
sig_CN_list <- unique(sig_CNs$neighborhood_cluster)

# ------------------------- Prepare Data for Plotting -------------------------
cn_freq_sig <- cn_freq %>%
  filter(neighborhood_cluster %in% sig_CN_list) %>%
  mutate(
    outcome_label = case_when(
      primary_outcome == 0 ~ "Dead/recurred in 60 months",
      primary_outcome == 1 ~ "Alive after 60 months",
      TRUE ~ NA_character_
    ),
    outcome_label = factor(outcome_label, levels = c("Dead/recurred in 60 months", "Alive after 60 months")),
    neighborhood_cluster = as.character(neighborhood_cluster)
  )

# ------------------------- Compute y position and FDR p-value label per CN -------------------------
p_labels <- cn_freq_sig %>%
  left_join(
    sig_CNs %>% mutate(neighborhood_cluster = as.character(neighborhood_cluster)) %>%
      select(neighborhood_cluster, p_adj),
    by = "neighborhood_cluster"
  ) %>%
  group_by(neighborhood_cluster) %>%
  summarise(
    y_pos = max(frequency, na.rm = TRUE) * 1.05,
    p_label = paste0("FDR p = ", signif(unique(p_adj), 2)),
    .groups = "drop"
  )

# ------------------------- Plot Faceted Summary -------------------------
ggplot(cn_freq_sig, aes(x = outcome_label, y = frequency, fill = outcome_label)) +
  geom_boxplot(width = 0.6, alpha = 0.9, outlier.shape = NA) +
  geom_jitter(width = 0.3, size = 1, alpha = 0.7, shape = 21, stroke = 0.2, color = "black") +
  facet_wrap(~ neighborhood_cluster, scales = "free_y", ncol = 3) +
  geom_text(
    data = p_labels,
    aes(x = 1.5, y = y_pos, label = p_label),
    inherit.aes = FALSE,
    size = 3.5
  ) +
  scale_fill_manual(
    values = c(
      "Dead/recurred in 60 months" = "red",
      "Alive after 60 months" = "blue"
    )
  ) +
  labs(
    title = "Significant CN Frequencies by 60-Month Survival",
    x = "Clinical Outcome",
    y = "CN Frequency",
    color = "Outcome Group"
  ) +
  theme_bw(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1),
    strip.text = element_text(size = 10),
    legend.position = "right"
  )

# ------------------------- Save Individual Boxplots -------------------------
for (cn in sig_CN_list) {
  cn_data <- cn_freq_sig %>% filter(neighborhood_cluster == as.character(cn))
  cn_label <- p_labels %>% filter(neighborhood_cluster == as.character(cn))
  
  p <- ggplot(cn_data, aes(x = outcome_label, y = frequency, fill = outcome_label)) +
    geom_boxplot(width = 0.6, alpha = 0.9, outlier.shape = NA) +
    geom_jitter(width = 0.3, size = 1, alpha = 0.7, shape = 21, stroke = 0.2, color = "black") +
    geom_text(
      data = cn_label,
      aes(x = 1.5, y = y_pos, label = p_label),
      inherit.aes = FALSE,
      size = 3.5
    ) +
    scale_fill_manual(
      values = c(
        "Dead/recurred in 60 months" = "red",
        "Alive after 60 months" = "blue"
      )
    ) +
    labs(
      title = paste("CN", cn, "- Frequency by 60-Month Survival"),
      x = "Clinical Outcome",
      y = "CN Frequency",
      fill = "Outcome Group"
    ) +
    theme_bw(base_size = 12) +
    theme(
      axis.text.x = element_text(angle = 30, hjust = 1),
      strip.text = element_text(size = 10),
      legend.position = "right"
    )
  
  ggsave(
    filename = paste0(output_dir, "CN_", cn, "_boxplot.png"),
    plot = p,
    width = 10,
    height = 5,
    dpi = 300
  )
}

# ------------------------- Count and Print Number of Samples per CN -------------------------
sample_counts_per_cn <- cn_freq_sig %>%
  count(neighborhood_cluster, name = "num_samples")

print(sample_counts_per_cn)

