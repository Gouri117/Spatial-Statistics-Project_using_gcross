# ------------------------- Load Required Libraries -------------------------
library(readr)
library(dplyr)
library(stringr)

# ------------------------- Step 1: Load CSV and Drop 'Other' Columns -------------------------

# Define metadata columns
meta_cols <- c("cell_id", "image", "window_center_id", "neighborhood_cluster")

# Load the full clustered file
cluster_file <- "/nfs/turbo/umms-ukarvind/shared_data/mfIHC_collab_data/stanford_crc/results/g_cross_auc_clustered_k20_all.csv"
cluster_df <- read_csv(cluster_file)

# Identify interaction columns (not metadata)
interaction_cols <- setdiff(colnames(cluster_df), meta_cols)

# Drop any columns that contain "Other"
interaction_cols_no_other <- interaction_cols[!grepl("Other", interaction_cols)]
cluster_df <- cluster_df[, c(meta_cols, interaction_cols_no_other)]

# ------------------------- Step 2: Normalize Dot Notation -------------------------

# Function to collapse multiple dots, remove trailing/leading ones
normalize_dots <- function(name) {
  name <- gsub("\\.+", ".", name)      # collapse multiple dots
  name <- gsub("^\\.|\\.$", "", name)  # remove leading/trailing dots
  name <- trimws(name)
  return(name)
}

# Normalize interaction column names
interaction_cols_original <- setdiff(colnames(cluster_df), meta_cols)
interaction_cols_normalized <- sapply(interaction_cols_original, normalize_dots, USE.NAMES = FALSE)

# Apply normalized names back to cluster_df
colnames(cluster_df)[match(interaction_cols_original, colnames(cluster_df))] <- interaction_cols_normalized

# Update interaction column names after normalization
interaction_cols_original <- interaction_cols_normalized

# ------------------------- Step 3: Parse and Rename as 'Source → Target' -------------------------

# Define final valid cell types including subtypes
valid_cell_types <- c(
  "Tumor 1", "Tumor 2 Ki67 Proliferating", "Tumor 3", "Tumor 4", "Tumor 5",
  "Tumor 6 DC", "Tumor 7", "CD4 T cell", "CD8 T cell", "B cell",
  "Blood vessel", "Granulocyte", "Macrophage", "Stroma"
)

# Function to parse and rename interaction columns to "Source → Target"
parse_interaction_name <- function(colname, valid_types) {
  terms <- unlist(strsplit(colname, "\\."))
  n <- length(terms)
  
  for (i in 1:(n - 1)) {
    source <- paste(terms[1:i], collapse = " ")
    target <- paste(terms[(i + 1):n], collapse = " ")
    source <- trimws(source)
    target <- trimws(target)
    
    if (source %in% valid_types && target %in% valid_types) {
      return(paste0(source, "→", target))
    }
  }
  return(NA)
}

# Apply parsing
renamed_cols <- sapply(interaction_cols_original, parse_interaction_name, valid_types = valid_cell_types)

# Keep only successfully parsed interactions
valid_idx <- which(!is.na(renamed_cols))
interaction_cols_final <- interaction_cols_original[valid_idx]
renamed_cols_final <- renamed_cols[valid_idx]

# Rename columns in the data frame
colnames(cluster_df)[match(interaction_cols_final, colnames(cluster_df))] <- renamed_cols_final
interaction_cols_clean <- renamed_cols_final  # final list of renamed valid interactions

# ------------------------- Confirm -------------------------
message(" Final number of valid 'Source → Target' interactions: ", length(interaction_cols_clean))
print(head(interaction_cols_clean, 10))

# ------------------------- Save Final Parsed Data -------------------------

# Define output path
output_file <- "/nfs/turbo/umms-ukarvind/shared_data/mfIHC_collab_data/stanford_crc/results/g_cross_cleaned_source_target.csv"

# Save to CSV
write_csv(cluster_df, output_file)

message(" cluster_df saved to: ", output_file)
