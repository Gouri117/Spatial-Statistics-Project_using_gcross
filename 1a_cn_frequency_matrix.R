library(dplyr)
library(tidyr)
library(readr)

# -------- Load Data --------
df <- read_csv("/nfs/turbo/umms-ukarvind/shared_data/mfIHC_collab_data/stanford_crc/results/g_cross_clustered_k20.csv")

# -------- Count Cells per CN per Sample --------
cn_counts <- df %>%
  group_by(image, cluster_k20) %>%
  summarise(count = n(), .groups = "drop")

# -------- Total Cells per Sample --------
sample_totals <- df %>%
  group_by(image) %>%
  summarise(total_cells = n(), .groups = "drop")

# -------- Merge and Calculate Proportions --------
cn_props <- cn_counts %>%
  left_join(sample_totals, by = "image") %>%
  mutate(proportion = count / total_cells)

# -------- Convert to CN × Sample Matrix --------
cn_matrix <- cn_props %>%
  select(image, cluster_k20, proportion) %>%
  pivot_wider(names_from = cluster_k20,
              values_from = proportion,
              values_fill = list(proportion = 0))

# -------- Optional: Save to CSV --------
write_csv(cn_matrix, "/nfs/turbo/umms-ukarvind/shared_data/mfIHC_collab_data/stanford_crc/results/cn_sample_matrix_k20.csv")
outcome_df <- read_csv("/nfs/turbo/umms-ukarvind/shared_data/mfIHC_collab_data/stanford_crc/charville_labels_samples_with_outcome.csv")

# ---------- Match and Join ----------
cn_matrix_annotated <- cn_matrix %>%
  left_join(outcome_df %>% select(sample, primary_outcome),
            by = c("image" = "sample"))

# ---------- Save to CSV ----------
write_csv(cn_matrix_annotated,
          "/nfs/turbo/umms-ukarvind/shared_data/mfIHC_collab_data/stanford_crc/results/cn_sample_matrix_k20_with_outcome.csv")
