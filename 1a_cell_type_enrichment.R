# ------------------------- Load Required Libraries -------------------------
library(dplyr)
library(readr)
library(stringr)
library(tidyr)

# ------------------------- Load Clustered G-Cross Data -------------------------
data_fd <- '/nfs/turbo/umms-ukarvind/shared_data/gcross_data/'
proj_fd <- 'example'
cluster_file <- paste0(data_fd, proj_fd, '/results/g_cross_auc_clustered_k20.csv')

cluster_df <- read_csv(cluster_file)

# ------------------------- Get Outgoing Cell Types -------------------------
interaction_cols <- colnames(cluster_df)[grepl("→", colnames(cluster_df))]
source_cell_types <- unique(str_extract(interaction_cols, "^[^→]+"))

# ------------------------- Compute Outgoing Composition per Row -------------------------
source_activity_df <- sapply(source_cell_types, function(ct) {
  cols <- grep(paste0("^", ct, "→"), colnames(cluster_df), value = TRUE)
  rowSums(cluster_df[, cols, drop = FALSE], na.rm = TRUE)
}) %>% as.data.frame()

# Attach cluster label
source_activity_df$neighborhood_cluster <- cluster_df$neighborhood_cluster

# ------------------------- Compute Average Composition per Cluster -------------------------
cluster_avg <- source_activity_df %>%
  group_by(neighborhood_cluster) %>%
  summarise(across(all_of(source_cell_types), mean, na.rm = TRUE), .groups = "drop")

# ------------------------- Compute Global Composition -------------------------
global_avg <- colMeans(source_activity_df[, source_cell_types], na.rm = TRUE)

# ------------------------- Compute log2 Fold Change Enrichment -------------------------
cluster_long <- cluster_avg %>%
  pivot_longer(cols = all_of(source_cell_types), names_to = "cell_type", values_to = "cluster_avg")

global_df <- data.frame(cell_type = names(global_avg), global_avg = as.numeric(global_avg))

enrichment_df <- cluster_long %>%
  left_join(global_df, by = "cell_type") %>%
  mutate(log2_fc = log2((cluster_avg + 1e-9) / (global_avg + 1e-9))) %>%
  arrange(desc(log2_fc))

# ------------------------- Save Outputs -------------------------
write_csv(enrichment_df, paste0(data_fd, proj_fd, "/results/cell_type_enrichment_full.csv"))


# ------------------------- Subset: Top 3 Cell Types per Cluster (No Threshold) -------------------------
top3_enrichment_df <- enrichment_df %>%
  group_by(neighborhood_cluster) %>%
  arrange(desc(log2_fc)) %>%
  slice_head(n = 3) %>%
  ungroup()

# ------------------------- Save Top 3 Enriched Cell Types per Cluster -------------------------
write_csv(top3_enrichment_df, paste0(data_fd, proj_fd, "/results/top3_enrichment_per_cluster.csv"))




