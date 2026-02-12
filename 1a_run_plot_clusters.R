####################################################################################
## Overview: Plot example point patterns by CN clusters
## Updated: May 29, 2025
####################################################################################

# -------------------- Load libraries --------------------
library(ggplot2)
library(dplyr)
library(tidyverse)
library(yaml)

# -------------------- Path setup --------------------
code_fd <- '/nfs/turbo/umms-ukarvind/gouri/gcross_pkg/'
proj_fd <- 'example'
data_fd <- '/nfs/turbo/umms-ukarvind/shared_data/gcross_data/'
csv_name <- 'subset_data.csv'

spatial_file <- paste0(data_fd, proj_fd, '/processed/', csv_name)
cluster_file <- paste0(data_fd, proj_fd, '/results/g_cross_auc_clustered_k30.csv')

# -------------------- Load parameter YAML --------------------
params <- yaml::yaml.load_file(paste0(code_fd, 'param_files/', proj_fd, '.yaml'))
for (param_name in names(params)) {
  assign(param_name, params[[param_name]])
}

# -------------------- Load spatial data and add cell_id --------------------
tbl <- read.csv(spatial_file, header = TRUE, check.names = TRUE)
tbl$cell_id <- as.character(seq_len(nrow(tbl)))

# -------------------- Load cluster assignments --------------------
clusters <- read.csv(cluster_file, header = TRUE)
clusters$cell_id <- as.character(clusters$cell_id)

# -------------------- Merge spatial data with clusters --------------------
merged_tbl <- left_join(tbl, clusters[, c("cell_id", "neighborhood_cluster")], by = "cell_id")

# -------------------- Plotting --------------------
sample_names <- unique(merged_tbl[[im_col]])
num_samples <- 1  # How many images to plot


for (sample in sample_names[1:num_samples]) {
  this_tbl <- merged_tbl[merged_tbl[[im_col]] == sample, ]
  this_tbl[[ct_col]][this_tbl[[ct_col]] == ""] <- "Other"
  
  p <- ggplot(this_tbl, aes_string(x = x_col, y = y_col, color = "factor(neighborhood_cluster)")) +
    geom_point(size = 1.5, alpha = 0.8) +
    theme_minimal() +
    labs(
      title = paste0("Sample: ", sample, 
                     "\nPatient: ", this_tbl[[pat_col]][1],
                     ", Pathology: ", this_tbl[[resp_col]][1]),
      color = "CN Cluster", x = "X", y = "Y"
    ) +
    theme(
      strip.text = element_text(face = "bold", size = 12),
      legend.position = "right"
    )
  
  print(p)
}

