library(dplyr)
library(readr)

# -------------------- Load Files --------------------
gcross1 <- read_csv("/nfs/turbo/umms-ukarvind/shared_data/mfIHC_collab_data/stanford_crc/results/g_cross_full_matrix_by_window_split1.csv")
gcross2 <- read_csv("/nfs/turbo/umms-ukarvind/shared_data/mfIHC_collab_data/stanford_crc/results/g_cross_full_matrix_by_window_split2.csv")

subset1 <- read_csv("/nfs/turbo/umms-ukarvind/shared_data/mfIHC_collab_data/stanford_crc/subset_data_samples_with_outcome_split1.csv")
subset2 <- read_csv("/nfs/turbo/umms-ukarvind/shared_data/mfIHC_collab_data/stanford_crc/subset_data_samples_with_outcome_split2.csv")

# -------------------- Recreate window_center_id in subset files --------------------
subset1_coords <- subset1 %>%
  group_by(sample) %>%
  mutate(window_center_id = row_number()) %>%
  ungroup() %>%
  select(image = sample, window_center_id, X, Y, primary_outcome, CELL_TYPE)

subset2_coords <- subset2 %>%
  group_by(sample) %>%
  mutate(window_center_id = row_number()) %>%
  ungroup() %>%
  select(image = sample, window_center_id, X, Y, primary_outcome, CELL_TYPE)

# -------------------- Join onto G-cross matrices --------------------
gcross1_with_xy <- gcross1 %>%
  left_join(subset1_coords, by = c("image", "window_center_id"))

gcross2_with_xy <- gcross2 %>%
  left_join(subset2_coords, by = c("image", "window_center_id"))

# -------------------- Merge --------------------
gcross_merged <- bind_rows(gcross1_with_xy, gcross2_with_xy)
write_csv(gcross_merged, "/nfs/turbo/umms-ukarvind/shared_data/mfIHC_collab_data/stanford_crc/results/g_cross_merged.csv")

# Count number of 0s and 1s in primary_outcome
gcross_merged %>%
  select(image, primary_outcome) %>%
  distinct() %>%                     # keep unique image-primary_outcome pairs
  count(primary_outcome)

