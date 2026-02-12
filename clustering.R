library(dplyr)
library(ggplot2)
library(cluster)    # for silhouette
library(factoextra) # for clustering visuals
library(readr)
library(umap)
library(klaR)       # optional for MiniBatch k-means

auc_matrix <- read_csv("/nfs/turbo/umms-ukarvind/shared_data/gcross_data/example/results/g_cross_full_matrix_by_window.csv")

# Separate metadata and features
cell_info <- auc_matrix %>% select(cell_id, image, window_center_id)
auc_data <- auc_matrix %>% select(-cell_id, -image, -window_center_id)

# Scale feature vectors
scaled_auc <- scale(auc_data)
