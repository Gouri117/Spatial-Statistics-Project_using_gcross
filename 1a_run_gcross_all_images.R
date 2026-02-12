# -------------------- Load libraries --------------------
library(pracma)
library(spatstat)
library(ggplot2)
library(dplyr)
library(tidyverse)
library(yaml)
library(RANN)
library(ggforce)
library(tidyr)
library(progress)

# -------------------- Paths and Config --------------------
code_fd <- '/nfs/turbo/umms-ukarvind/gouri/gcross_pkg/'

# Define project-specific path components
proj_fd <- 'stanford_crc'
data_file <- '/nfs/turbo/umms-ukarvind/shared_data/mfIHC_collab_data/stanford_crc/subset_data_samples_with_outcome.csv'

# YAML param file is reused from 'example' project
param_file <- file.path(code_fd, 'param_files', 'example.yaml')

# -------------------- Load Functions and Parameters --------------------
source(file.path(code_fd, 'main_code', 'gcross_functions.R'))

params <- yaml::yaml.load_file(param_file)
for (param_name in names(params)) {
  assign(param_name, params[[param_name]])
}

# -------------------- G-cross Parameters --------------------
R <- 0:60
radius_px <- 60
min_neighbors <- 2
local_padding <- 60

# -------------------- Load Data --------------------
tbl <- read.table(data_file, header = TRUE, sep = ",", comment.char = "%", check.names = TRUE)
tbl$cell_id <- paste0(seq_len(nrow(tbl)))

# Split image data
img_list <- split(tbl, f = tbl[[im_col]])
all_types <- unique(tbl[[ct_col]])

# -------------------- Output Directory --------------------
out_dir <- file.path(dirname(data_file), 'results')
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# -------------------- Initialize Storage --------------------
g_auc_full_list <- list()
composition_long <- list()
row_counter <- 1

# -------------------- Loop Through All Images --------------------
message("\n--- Running full G-cross AUC extraction across all images ---\n")
total_images <- length(img_list)
image_counter <- 1

for (image_name in names(img_list)) {
  message(paste0("Processing image ", image_counter, "/", total_images, ": ", image_name))
  
  curr_img <- img_list[[image_name]]
  x_vals <- curr_img[[x_col]]
  y_vals <- curr_img[[y_col]]
  labels <- as.character(curr_img[[ct_col]])
  cell_ids <- curr_img$cell_id
  
  radius_search <- nn2(
    data = cbind(x_vals, y_vals),
    query = cbind(x_vals, y_vals),
    k = nrow(curr_img),
    searchtype = "radius",
    radius = radius_px
  )
  
  pb_full <- progress_bar$new(format = paste0(image_name, " [:bar] :current/:total (:percent) ETA: :eta"),
                              total = nrow(curr_img), clear = FALSE, width = 80)
  
  for (i in seq_len(nrow(curr_img))) {
    pb_full$tick()
    
    if (i %% 100 == 0) {
      message(paste0("  - Processed ", i, "/", nrow(curr_img), " windows for image ", image_name))
    }
    
    neighbors_idx <- radius_search$nn.idx[i, ]
    neighbors_idx <- neighbors_idx[neighbors_idx != 0]
    if (length(neighbors_idx) < min_neighbors) next
    
    window_cells <- curr_img[neighbors_idx, ]
    total_cells <- nrow(window_cells)
    
    comp_table <- table(window_cells[[ct_col]])
    comp_df <- data.frame(
      image = image_name,
      window_center_id = i,
      cell_id = curr_img$cell_id[i],
      cell_type = names(comp_table),
      count = as.integer(comp_table),
      n_neighbors = total_cells,
      composition = round(as.integer(comp_table) / total_cells, 5)
    )
    composition_long[[length(composition_long) + 1]] <- comp_df
    
    x_range <- range(window_cells[[x_col]])
    y_range <- range(window_cells[[y_col]])
    padded_window <- owin(
      xrange = c(x_range[1] - local_padding, x_range[2] + local_padding),
      yrange = c(y_range[1] - local_padding, y_range[2] + local_padding)
    )
    
    window_ppp <- ppp(
      x = window_cells[[x_col]],
      y = window_cells[[y_col]],
      marks = as.factor(window_cells[[ct_col]]),
      window = padded_window
    )
    
    for (local_center_idx in seq_len(nrow(window_cells))) {
      local_center <- window_cells[local_center_idx, ]
      local_center_type <- local_center[[ct_col]]
      
      for (target_type in all_types) {
        if (local_center_type == target_type && sum(window_cells[[ct_col]] == target_type) < 2) {
          auc_val <- 0
        } else if (!(target_type %in% window_cells[[ct_col]])) {
          auc_val <- 0
        } else {
          g_obj <- tryCatch({
            Gcross(window_ppp, i = local_center_type, j = target_type, r = R, correction = "rs")
          }, error = function(e) return(NULL))
          
          if (is.null(g_obj) || all(is.na(g_obj$rs)) || all(is.nan(g_obj$rs))) {
            auc_val <- 0
          } else {
            auc_val <- trapz(R, g_obj$rs)
            if (is.nan(auc_val)) auc_val <- 0
          }
        }
        
        g_auc_full_list[[row_counter]] <- data.frame(
          image = image_name,
          cell_id = curr_img$cell_id[i],
          window_center_id = i,
          interaction_pair = paste(local_center_type, target_type, sep = "→"),
          auc = auc_val
        )
        row_counter <- row_counter + 1
      }
    }
  }
  
  message(paste0("Finished processing image: ", image_name, " (", image_counter, "/", total_images, ")"))
  image_counter <- image_counter + 1
}

# -------------------- Save AUC Matrix --------------------
message("\n--- Saving G-cross AUC matrix to file ---\n")
g_auc_full_df <- do.call(rbind, g_auc_full_list)
g_auc_matrix_full <- g_auc_full_df %>%
  group_by(cell_id, image, window_center_id, interaction_pair) %>%
  summarise(auc = mean(as.numeric(auc), na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(
    names_from = interaction_pair,
    values_from = auc,
    values_fill = list(auc = 0)
  )
write.csv(g_auc_matrix_full, file.path(out_dir, "g_cross_full_matrix_by_window_all_samples.csv"), row.names = FALSE)
message("G-cross AUC matrix saved.")

# -------------------- Save Composition Table --------------------
message("\n--- Saving composition table ---\n")
composition_df <- do.call(rbind, composition_long)
write.csv(composition_df, file.path(out_dir, 'composition_by_window_all_samples.csv'), row.names = FALSE)
message("Composition table saved.")
