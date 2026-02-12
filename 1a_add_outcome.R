# ------------------------- Load Required Libraries -------------------------
library(readr)
library(dplyr)

# ------------------------- Load Clustered G-Cross Data -------------------------
cluster_file <- "/nfs/turbo/umms-ukarvind/shared_data/mfIHC_collab_data/stanford_crc/results/g_cross_cleaned_source_target.csv"
cluster_df <- read_csv(cluster_file)

# ------------------------- Load Outcome Metadata -------------------------
outcome_file <- "/nfs/turbo/umms-ukarvind/shared_data/mfIHC_collab_data/stanford_crc/subset_data_samples_with_outcome.csv"
outcome_df <- read_csv(outcome_file)

# ------------------------- Normalize column names -------------------------
outcome_df <- outcome_df %>% rename(cell_id = CELL_ID)

# ------------------------- Merge on cell_id -------------------------
merged_df <- left_join(cluster_df, outcome_df, by = "cell_id")

# ------------------------- Save Merged Output -------------------------
output_file <- "/nfs/turbo/umms-ukarvind/shared_data/mfIHC_collab_data/stanford_crc/results/g_cross_with_outcome.csv"
write_csv(merged_df, output_file)

message(" Merged cluster_df with outcome data and saved to: ", output_file)

