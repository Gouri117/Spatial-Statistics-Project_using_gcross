# ------------------------- Load Libraries -------------------------
library(dplyr)
library(tidyr)
library(readr)
library(stringr)

# ------------------------- Load Cleaned File -------------------------
input_file <- "/nfs/turbo/umms-ukarvind/shared_data/mfIHC_collab_data/stanford_crc/results/g_cross_cleaned_source_target.csv"
cluster_df <- read_csv(input_file)

# ------------------------- Identify Interaction Columns -------------------------
meta_cols <- c("cell_id", "image", "window_center_id", "neighborhood_cluster")
interaction_cols <- setdiff(colnames(cluster_df), meta_cols)

# ------------------------- Compute Mean AUC per Interaction per Cluster -------------------------
pairwise_avg <- cluster_df %>%
  group_by(neighborhood_cluster) %>%
  summarise(across(all_of(interaction_cols), mean, na.rm = TRUE), .groups = "drop")

# ------------------------- Compute Global Mean per Interaction -------------------------
global_avg <- colMeans(cluster_df[, interaction_cols], na.rm = TRUE)

# ------------------------- Pivot for log2 Fold-Change -------------------------
pairwise_long <- pairwise_avg %>%
  pivot_longer(cols = all_of(interaction_cols),
               names_to = "interaction_pair",
               values_to = "cluster_mean")

global_df <- data.frame(interaction_pair = names(global_avg), global_mean = as.numeric(global_avg))

pairwise_enrichment <- pairwise_long %>%
  left_join(global_df, by = "interaction_pair") %>%
  mutate(log2_fc = log2((cluster_mean + 1e-9) / (global_mean + 1e-9))) %>%
  arrange(desc(log2_fc))

# ------------------------- Save Output -------------------------
out_dir <- "/nfs/turbo/umms-ukarvind/shared_data/mfIHC_collab_data/stanford_crc/results"
write_csv(pairwise_enrichment, file.path(out_dir, "pairwise_interaction_enrichment.csv"))

# ------------------------- Top 3 Per Cluster -------------------------
top3_pairs_per_cluster <- pairwise_enrichment %>%
  group_by(neighborhood_cluster) %>%
  arrange(desc(log2_fc)) %>%
  slice_head(n = 3) %>%
  ungroup()

write_csv(top3_pairs_per_cluster, file.path(out_dir, "top3_interactions_per_cluster.csv"))

message(" Pairwise enrichment complete. Files saved to: ", out_dir)
