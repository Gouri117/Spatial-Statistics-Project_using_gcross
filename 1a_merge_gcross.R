library(dplyr)
library(readr)

# -------------------- Load G-cross CSVs --------------------
file1 <- "/nfs/turbo/umms-ukarvind/shared_data/mfIHC_collab_data/stanford_crc/results/g_cross_full_matrix_by_window_split1.csv"
file2 <- "/nfs/turbo/umms-ukarvind/shared_data/mfIHC_collab_data/stanford_crc/results/g_cross_full_matrix_by_window_split2.csv"

df1 <- read_csv(file1)
df2 <- read_csv(file2)

# -------------------- Merge Datasets --------------------
merged_df <- bind_rows(df1, df2)

# -------------------- Count Unique Samples --------------------
# Assumes 'image' column represents sample
n_unique_samples <- merged_df %>% pull(image) %>% unique() %>% length()

cat("Number of unique samples:", n_unique_samples, "\n")
