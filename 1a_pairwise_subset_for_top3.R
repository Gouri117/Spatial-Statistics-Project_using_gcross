# ------------------------- Load Required Libraries -------------------------
library(dplyr)
library(readr)
library(purrr)
library(tidyr)
library(stringr)

# ------------------------- Paths -------------------------
data_fd <- '/nfs/turbo/umms-ukarvind/shared_data/gcross_data/'
proj_fd <- 'example'

cell_enrichment_file <- paste0(data_fd, proj_fd, "/results/cell_type_enrichment_k20.csv")
pairwise_file <- paste0(data_fd, proj_fd, "/results/pairwise_enrichment_full.csv")
output_dir <- paste0(data_fd, proj_fd, "/results/")
output_path <- paste0(output_dir, "pairwise_top3_subset.csv")

# ------------------------- Load the Data -------------------------
cell_type_enrichment <- read_csv(cell_enrichment_file)
pairwise_df <- read_csv(pairwise_file)

# ------------------------- Step 1: Subset Top 3 Cell Types per Cluster -------------------------
top3_df <- cell_type_enrichment %>%
  group_by(neighborhood_cluster) %>%
  arrange(desc(log2_fc)) %>%
  slice_head(n = 3) %>%
  ungroup()

# ------------------------- Step 2: Generate All Pairwise and Self-Interactions -------------------------
top3_pairs <- top3_df %>%
  group_by(neighborhood_cluster) %>%
  summarise(
    top3 = list(CELL_TYPE),
    .groups = "drop"
  ) %>%
  mutate(
    interactions = map(top3, ~ {
      ct <- .x
      # Generate all combinations including self-pairs
      expand.grid(from = ct, to = ct, stringsAsFactors = FALSE) %>%
        mutate(pair = paste0(from, "→", to)) %>%
        pull(pair)
    })
  ) %>%
  select(neighborhood_cluster, interactions) %>%
  unnest(interactions)

# ------------------------- Step 3: Subset Pairwise File -------------------------
subset_pairwise <- pairwise_df %>%
  inner_join(top3_pairs, by = c("neighborhood_cluster", "pair" = "interactions"))

# ------------------------- Step 4: Ensure Output Directory Exists -------------------------
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# ------------------------- Step 5: Save Output -------------------------
write_csv(subset_pairwise, output_path)

# ------------------------- Confirm -------------------------
message("Saved subset with self-pairs included: ", output_path)



