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
proj_fd <- 'example'
data_fd <- '/nfs/turbo/umms-ukarvind/shared_data/gcross_data/'
csv_name <- 'subset_data.csv'
csv_file <- paste0(data_fd, '/', proj_fd, '/processed/', csv_name)

source(paste0(code_fd, 'main_code/gcross_functions.R'))

params <- yaml::yaml.load_file(paste0(code_fd, 'param_files/', proj_fd, '.yaml'))
for (param_name in names(params)) {
  assign(param_name, params[[param_name]])
}

R <- 0:60
radius_px <- 60
min_neighbors <- 2
max_pad_local <- 60

# -------------------- Load Data --------------------
tbl <- read.table(csv_file, header = TRUE, sep = ",", comment.char = "%", check.names = TRUE)
img_list <- split(tbl, f = tbl[[im_col]])
image_to_run <- unique(tbl[[im_col]])[1]
curr_img <- img_list[[image_to_run]]

x_vals <- curr_img[[x_col]]
y_vals <- curr_img[[y_col]]
labels <- as.character(curr_img[[ct_col]])

# -------------------- Output Dirs --------------------
out_dir <- paste0(data_fd, '/', proj_fd, '/results/')
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# -------------------- Radius Search --------------------
radius_search <- nn2(
  data = cbind(x_vals, y_vals),
  query = cbind(x_vals, y_vals),
  k = nrow(curr_img),
  searchtype = "radius",
  radius = radius_px
)

# -------------------- Initialize Storage --------------------
all_types <- unique(labels)
g_auc_full_list <- list()
composition_long <- list()
row_counter <- 1

# -------------------- G-cross AUC Calculation --------------------
message("\n--- Running full G-cross AUC extraction for all windows with padding = 60 ---\n")
pb_full <- progress_bar$new(format = "[:bar] :current/:total (:percent) ETA: :eta", total = nrow(curr_img), clear = FALSE, width = 60)

for (i in seq_len(nrow(curr_img))) {
  pb_full$tick()
  neighbors_idx <- radius_search$nn.idx[i, ]
  neighbors_idx <- neighbors_idx[neighbors_idx != 0]
  if (length(neighbors_idx) < min_neighbors) next
  
  window_cells <- curr_img[neighbors_idx, ]
  total_cells <- nrow(window_cells)
  
  # Composition calculation
  comp_table <- table(window_cells[[ct_col]])
  comp_df <- data.frame(
    window_center_id = i,
    cell_type = names(comp_table),
    count = as.integer(comp_table),
    n_neighbors = total_cells,
    composition = round(as.integer(comp_table) / total_cells, 5)
  )
  composition_long[[i]] <- comp_df
  
  # Define padded observation window
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
        window_center_id = i,
        interaction_pair = paste(local_center_type, target_type, sep = "→"),
        auc = auc_val
      )
      row_counter <- row_counter + 1
    }
  }
}

# -------------------- Save AUC Matrix --------------------
g_auc_full_df <- do.call(rbind, g_auc_full_list)
g_auc_matrix_full <- g_auc_full_df %>%
  group_by(window_center_id, interaction_pair) %>%
  summarise(auc = mean(as.numeric(auc), na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(
    names_from = interaction_pair,
    values_from = auc,
    values_fill = list(auc = 0)
  )
write.csv(g_auc_matrix_full, paste0(out_dir, "g_cross_full_matrix_by_window.csv"), row.names = FALSE)

# -------------------- Save Composition Table --------------------
composition_df <- do.call(rbind, composition_long)
write.csv(composition_df, paste0(out_dir, 'composition_by_window.csv'), row.names = FALSE)

# -------------------- Window Plot Function --------------------
plot_window_gg <- function(i, curr_img, radius_search, x_col, y_col, ct_col, radius_px) {
  center_cell <- curr_img[i, ]
  center_x <- center_cell[[x_col]]
  center_y <- center_cell[[y_col]]
  center_type <- center_cell[[ct_col]]
  
  neighbors_idx <- radius_search$nn.idx[i, ]
  neighbors_idx <- neighbors_idx[neighbors_idx != 0]
  window_cells <- curr_img[neighbors_idx, ]
  window_cells$role <- ifelse(rownames(window_cells) == rownames(center_cell), "center", "neighbor")
  
  ggplot(window_cells, aes_string(x = x_col, y = y_col)) +
    geom_point(aes(color = .data[[ct_col]], shape = role), size = 3, alpha = 0.8) +
    scale_shape_manual(values = c(center = 8, neighbor = 16)) +
    geom_point(data = center_cell, aes_string(x = x_col, y = y_col), color = "black", size = 5, shape = 4, stroke = 1.2) +
    ggforce::geom_circle(aes(x0 = center_x, y0 = center_y, r = radius_px), inherit.aes = FALSE, color = "black", linetype = "dashed") +
    theme_minimal() +
    labs(title = paste("Window around cell", i, "-", center_type),
         subtitle = paste("Radius:", radius_px, "px"), color = "Cell Type")
}

# Example: visualize window 1
plot_window_gg(i = 1000, curr_img = curr_img, radius_search = radius_search, x_col = x_col, y_col = y_col, ct_col = ct_col, radius_px = radius_px)
