# ------------------------- Load Libraries -------------------------
library(dplyr)
library(tidyr)
library(readr)
library(stringr)

# ------------------------- Load Cleaned File -------------------------
input_file <- "/nfs/turbo/umms-ukarvind/shared_data/mfIHC_collab_data/stanford_crc/results/g_cross_clustered_multiple_K.csv"
cluster_df <- read_csv(input_file)

# ------------------------- Identify Metadata and Interaction Columns -------------------------
meta_cols <- c("cell_id", "image", "window_center_id", "cluster_k25", "CELL_TYPE", "X", "Y", "primary_outcome")
exclude_cols <- c(meta_cols, "cluster_k20", "cluster_k30", "cluster_k40", "cluster_k45")

interaction_cols_clean <- setdiff(colnames(cluster_df), exclude_cols)

# ------------------------- Extract Source Cell Types -------------------------
source_cell_types <- unique(str_extract(interaction_cols_clean, "^[^→]+"))
source_cell_types <- source_cell_types[!is.na(source_cell_types) & source_cell_types != "Other"]

# ------------------------- Compute Outgoing Activity per Row -------------------------
source_activity_df <- sapply(source_cell_types, function(ct) {
  cols <- grep(paste0("^", ct, "→"), interaction_cols_clean, value = TRUE)
  if (length(cols) == 0) return(rep(0, nrow(cluster_df)))  # fallback if no columns match
  rowSums(cluster_df[, cols, drop = FALSE], na.rm = TRUE)
}) %>% as.data.frame()

# Attach cluster labels
source_activity_df$cluster_k25 <- cluster_df$cluster_k25

# ------------------------- Compute Average Activity per Cluster -------------------------
cluster_avg <- source_activity_df %>%
  group_by(cluster_k25) %>%
  summarise(across(all_of(source_cell_types), ~mean(.x, na.rm = TRUE)), .groups = "drop")

# ------------------------- Compute Global Average -------------------------
global_avg <- colMeans(source_activity_df[, source_cell_types], na.rm = TRUE)

# ------------------------- Compute log2 Fold-Change -------------------------
cluster_long <- cluster_avg %>%
  pivot_longer(cols = all_of(source_cell_types), names_to = "cell_type", values_to = "cluster_avg")

global_df <- data.frame(cell_type = names(global_avg), global_avg = as.numeric(global_avg))

enrichment_df <- cluster_long %>%
  left_join(global_df, by = "cell_type") %>%
  mutate(log2_fc = log2((cluster_avg + 1e-9) / (global_avg + 1e-9))) %>%
  arrange(desc(log2_fc))

# ------------------------- Save Enrichment Table -------------------------
out_dir <- "/nfs/turbo/umms-ukarvind/shared_data/mfIHC_collab_data/stanford_crc/results"
write_csv(enrichment_df, file.path(out_dir, "cell_type_enrichment_full.csv"))

# ------------------------- Extract Top 3 Enriched Cell Types per Cluster -------------------------
top3_enrichment_df <- enrichment_df %>%
  group_by(cluster_k25) %>%
  arrange(desc(log2_fc)) %>%
  slice_head(n = 3) %>%
  ungroup()

write_csv(top3_enrichment_df, file.path(out_dir, "top3_enrichment_per_cluster.csv"))

message(" Enrichment complete. Files saved to: ", out_dir)
