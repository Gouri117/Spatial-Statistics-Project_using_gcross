library(dplyr)
library(readr)
library(ggplot2)
library(ggforce)

# -------------------- Load Cluster and Cell Type Data --------------------
df_cluster <- read_csv("/nfs/turbo/umms-ukarvind/shared_data/mfIHC_collab_data/stanford_crc/results/g_cross_clustered_k30.csv")
df_celltype <- read_csv("/nfs/turbo/umms-ukarvind/shared_data/mfIHC_collab_data/stanford_crc/subset_data_samples_with_outcome_tumor_category.csv")

# -------------------- Create Output Folders --------------------
output_dir_clusters <- "/nfs/turbo/umms-ukarvind/shared_data/mfIHC_collab_data/stanford_crc/results/plots_by_image/cluster"
output_dir_celltype <- "/nfs/turbo/umms-ukarvind/shared_data/mfIHC_collab_data/stanford_crc/results/plots_by_image/celltype"
dir.create(output_dir_clusters, showWarnings = FALSE, recursive = TRUE)
dir.create(output_dir_celltype, showWarnings = FALSE, recursive = TRUE)

# -------------------- Plot by Cluster --------------------
for (img in unique(df_cluster$image)) {
  df_img <- df_cluster %>% filter(image == img)
  
  p <- ggplot(df_img, aes(x = x, y = y, color = factor(cluster))) +
    geom_point(size = 0.5, alpha = 0.8) +
    coord_fixed() +
    theme_void() +
    scale_color_viridis_d(name = "Cluster") +
    ggtitle(paste("Image:", img, "- Clusters")) +
    theme(legend.position = "right", plot.title = element_text(hjust = 0.5))
  
  ggsave(file.path(output_dir_clusters, paste0(img, "_clusters.png")),
         plot = p, width = 6, height = 5, dpi = 300)
}

# -------------------- Plot by Cell Type --------------------
for (img in unique(df_celltype$sample)) {
  df_img <- df_celltype %>% filter(sample == img)
  
  p <- ggplot(df_img, aes(x = X, y = Y, color = CELL_TYPE)) +
    geom_point(size = 0.5, alpha = 0.8) +
    coord_fixed() +
    theme_void() +
    scale_color_manual(values = scales::hue_pal()(length(unique(df_celltype$CELL_TYPE)))) +
    ggtitle(paste("Image:", img, "- Cell Types")) +
    theme(legend.position = "right", plot.title = element_text(hjust = 0.5))
  
  ggsave(file.path(output_dir_celltype, paste0(img, "_celltypes.png")),
         plot = p, width = 6, height = 5, dpi = 300)
}


