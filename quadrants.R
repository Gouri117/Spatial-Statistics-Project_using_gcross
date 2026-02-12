###################################################################################################
# Divides images into quadrants and tries to match CNs with patient outcomes

# -------------------- SETUP --------------------
library(tidyverse)
library(ggthemes)

# Step 1: Load data
file_path <- "/nfs/turbo/umms-ukarvind/shared_data/mfIHC_collab_data/stanford_crc/results/g_cross_clustered_multiple_K.csv"
df <- read_csv(file_path)

# Step 2: Assign quadrants within each image
df <- df %>%
  group_by(image) %>%
  mutate(
    quadrant = case_when(
      X <= median(X, na.rm = TRUE) & Y <= median(Y, na.rm = TRUE) ~ "Q1",
      X >  median(X, na.rm = TRUE) & Y <= median(Y, na.rm = TRUE) ~ "Q2",
      X <= median(X, na.rm = TRUE) & Y >  median(Y, na.rm = TRUE) ~ "Q3",
      X >  median(X, na.rm = TRUE) & Y >  median(Y, na.rm = TRUE) ~ "Q4"
    )
  ) %>% ungroup()

# Step 3: Create quadrant-level CN composition table
cn_counts <- df %>%
  group_by(image, quadrant, cluster_k25) %>%
  summarise(count = n(), .groups = "drop")

# Step 4: Calculate proportions per quadrant
cn_props <- cn_counts %>%
  group_by(image, quadrant) %>%
  mutate(prop = count / sum(count)) %>%
  ungroup()

# Step 5: Merge with primary outcome
# Get unique outcome per image (assuming it's the same for all cells in one image)
outcome_df <- df %>%
  distinct(image, primary_outcome)

cn_props <- cn_props %>%
  left_join(outcome_df, by = "image")

# Step 6: Reshape data to wide format for plotting/statistics
cn_wide <- cn_props %>%
  pivot_wider(
    id_cols = c(image, quadrant, primary_outcome),
    names_from = cluster_k25,
    values_from = prop,
    values_fill = 0  # if a CN is absent, set its prop to 0
  )

# Step 7: Run Wilcoxon test for each CN (column)
wilcox_results <- cn_wide %>%
  pivot_longer(cols = -c(image, quadrant, primary_outcome), names_to = "CN", values_to = "prop") %>%
  group_by(CN) %>%
  summarise(
    p_value = wilcox.test(prop ~ primary_outcome)$p.value,
    .groups = "drop"
  ) %>%
  arrange(p_value)

# Optional: Add FDR correction
wilcox_results <- wilcox_results %>%
  mutate(p_adj = p.adjust(p_value, method = "fdr"))

# Step 8: View significant CNs (FDR < 0.05)
significant_cns_table <- wilcox_results %>%
  filter(p_adj < 0.05)

print(significant_cns_table)

# Format
df <- df %>%
  mutate(
    cluster_k25 = as.factor(cluster_k25),
    primary_outcome = as.factor(primary_outcome)
  )

# Significant CNs
significant_cns <- c("6", "5", "23", "16", "8")


# Frequency per image per CN
cn_image_freqs <- df %>%
  filter(cluster_k25 %in% significant_cns) %>%
  group_by(image, primary_outcome, cluster_k25) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(image) %>%
  mutate(freq = count / sum(count))

# Plot all in one faceted figure
ggplot(cn_image_freqs, aes(x = primary_outcome, y = freq, fill = primary_outcome)) +
  geom_boxplot(alpha = 0.8, outlier.shape = NA) +
  geom_jitter(width = 0.1, size = 1.5, alpha = 0.6, color = "black") +
  facet_wrap(~ cluster_k25, nrow = 1, scales = "free_y") +
  labs(
    title = "Relative Frequency of Significant CNs by Outcome",
    x = "Patient Outcome",
    y = "Relative Frequency"
  ) +
  scale_fill_brewer(palette = "Set2") +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "none",
    strip.text = element_text(size = 12, face = "bold")
  )

