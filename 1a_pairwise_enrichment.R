# ------------------------- Load Libraries -------------------------
library(dplyr)
library(tidyr)
library(readr)
library(stringr)

# ------------------------- Load Files -------------------------
input_file <- "/nfs/turbo/umms-ukarvind/shared_data/mfIHC_collab_data/stanford_crc/results/g_cross_clustered_multiple_K.csv"
top3_file <- "/nfs/turbo/umms-ukarvind/shared_data/mfIHC_collab_data/stanford_crc/results/top3_enrichment_per_cluster.csv"
cluster_df <- read_csv(input_file)
top3_df <- read_csv(top3_file)

# ------------------------- Identify Interaction Columns -------------------------
meta_cols <- c("cell_id", "image", "window_center_id", "cluster_k25", "CELL_TYPE", "X", "Y", "primary_outcome")
exclude_cols <- c(meta_cols, "cluster_k20", "cluster_k30", "cluster_k40", "cluster_k45")
interaction_cols <- setdiff(colnames(cluster_df), exclude_cols)

# ------------------------- Compute Mean AUC per Interaction per Cluster -------------------------
pairwise_avg <- cluster_df %>%
  group_by(cluster_k25) %>%
  summarise(across(all_of(interaction_cols), mean, na.rm = TRUE), .groups = "drop")

# ------------------------- Compute Global Mean -------------------------
global_avg <- colMeans(cluster_df[, interaction_cols], na.rm = TRUE)
global_df <- data.frame(interaction = names(global_avg), global_mean = as.numeric(global_avg))

# ------------------------- Pivot and Compute log2FC -------------------------
pairwise_long <- pairwise_avg %>%
  pivot_longer(cols = all_of(interaction_cols),
               names_to = "interaction",
               values_to = "cluster_mean") %>%
  left_join(global_df, by = "interaction") %>%
  mutate(log2_fc = log2((cluster_mean + 1e-9) / (global_mean + 1e-9)))

# ------------------------- Build 3x3 Interaction Grid from Top 3 Cell Types -------------------------
top3_cells_per_cluster <- top3_df %>%
  group_by(cluster_k25) %>%
  summarise(cell_types = list(unique(cell_type)), .groups = "drop")

top3_interactions <- top3_cells_per_cluster %>%
  rowwise() %>%
  mutate(pairs = list(expand.grid(source = cell_types, target = cell_types) %>%
                        mutate(interaction = paste0(source, "→", target)))) %>%
  unnest(pairs)


# ------------------------- Filter Final 9 Interactions per Cluster -------------------------
filtered_interactions <- pairwise_long %>%
  inner_join(top3_interactions %>% select(cluster_k25, interaction), by = c("cluster_k25", "interaction"))

# ------------------------- Save Output -------------------------
out_file <- "/nfs/turbo/umms-ukarvind/shared_data/mfIHC_collab_data/stanford_crc/results/top9_pairwise_interactions_per_cluster.csv"
write_csv(filtered_interactions, out_file)

message(" Saved top 9 interactions per cluster to: ", out_file)

message("Saved top 9 interactions for selected clusters to: ", selected_out_file)

