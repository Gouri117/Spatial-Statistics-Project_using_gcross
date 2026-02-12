####################################################################################
## overview: general script to plot gcross aucs
## created: 4/21/25

## prior requirements:
## 1. 1a_run_gcross.R has been run on mfIHC data

## script inputs:
## 1. processed data file: {path_to_data_folder}/subset_data.csv
## 2. proj_folder/output/1a_run_gcross/gcross_aucs_{imwise/patwise}_{aucobs/ratio}.Rdata

## script outputs:
## 1. auc plots for each rad: proj_folder/output/1b_auc_plot/{imwise/patwise}_rad{rad}.png

####################################################################################
## environment setup ##
library(ggplot2)
library(tidyverse)
library(grid)
library(ggpubr)
library(rstatix)
library(yaml)

####################################################################################
## path/folder setup -- set these variables! ##

## NOTE: directories should be same as what was set in 1a_run_gcross.R
code_fd <- '/nfs/turbo/umms-ukarvind/gouri/gcross_pkg/'         # path to the gcross_pkg code folder; eg: /nfs/turbo/umms-ukarvind/{your_fd}/gcross_pkg/
proj_fd <- 'example'  # project folder/name; eg: stanford_crc

data_fd <- '/nfs/turbo/umms-ukarvind/shared_data/gcross_data/'         # path to general data folder; eg: /nfs/turbo/umms-ukarvind/shared_data/gcross_data/ 
csv_name <- 'subset_data.csv' # name of csv data file; eg: subset_data.csv
csv_file <- paste0(data_fd, '/', proj_fd, '/processed/', csv_name)  # path to csv file; NOTE: if folder is structured as described in protocol, no need to edit thiss

## NOTE: and set this one variable
auc_type <- 'standard' # either 'standard' or 'ratio'

####################################################################################
## load and prepare data ##

# load the data table
tbl <- read.table(csv_file, header=TRUE,sep=",", comment.char = "%", check.names=TRUE) # the '#' character is read as special, so change comment.char to something else

# load the parameter file
params <- yaml.load_file(paste0(code_fd, 'param_files/', proj_fd, '.yaml'))
for (param_name in names(params)) {
  assign(param_name, params[[param_name]])
}

# load the gcross auc data
out_bp <- paste0(data_fd, proj_fd, '/output/')
load(paste0(out_bp, '1a_run_gcross/gcross_aucs_', analysis_type, '_', auc_type, '.RData')) # load the auc data

# define the main sub_col to which the analysis will be conducted on - either by images or patients
if (analysis_type == 'imwise') {
  sub_col <- im_col # colname to subset by
  scale_pad <- 0 # don't need any padding if just considering single image
} else if (analysis_type == 'patwise') { # otherwise, patient-wise analysis
  sub_col <- pat_col # colname to subset by
  if (coords_type == 'global') {
    im_col <- pat_col # colname to subset by
  } else if (has_global_coords == TRUE) {
    global_coords <- convert_to_global_coords(tbl, im_col, x_col, y_col, stamp_coord_type, xy_coord_type, stamp_coord_loc)
    tbl[[x_col]] <- global_coords[[x_col]]
    tbl[[y_col]] <- global_coords[[y_col]]
    im_col <- pat_col
  }
} else {
  stop('invalid analysis type')
}

# subset metadata/response columns of interest
if (length(other_cols) == 0) {
  pat_resp <- unique(tbl[,c(sub_col, resp_col)]) 
} else {
  pat_resp <- unique(tbl[,c(sub_col, resp_col, other_cols)]) 
}

# if response column is numeric, convert to character
if (is.numeric(pat_resp[[resp_col]])) {
  pat_resp[[resp_col]] <- as.character(pat_resp[[resp_col]])
}

cell_pairs <- names(gcross_aucs) # list of the cell pairs
rComp <- colnames(gcross_aucs[[1]])[-1] # get list of radii at which aucs were computed at; character list of rads

####################################################################################
## plotting functions ##

plot_bxp <- function(df_merged, resp_col, out_dir, rad, save_append='') { 
  
  title <- paste0('Gcross Infiltration for Cell Pairs')
  
  bxp <- ggboxplot(df_merged, x='cell_pair', y='AUC', color=resp_col) + 
    ylab('AUC') +
    labs(title=paste(title, ' at rad=',rad, sep='')) +
    theme(text=element_text(size=10),
          axis.text.x = element_text(angle=90, vjust=0.5),
          plot.title = element_text(hjust = 0.5))
  
  # perform statistical test
  stat.test <- df_merged %>%
    group_by(cell_pair) %>%
    wilcox_test(formula(paste("AUC ~", resp_col)), detailed= TRUE) %>%
    adjust_pvalue(method = "fdr") %>%
    add_significance("p.adj") %>%
    add_xy_position(x = 'cell_pair')
  print(stat.test)

  # write the significance test results
  stat.test_sub <- stat.test[, !names(stat.test) %in% c('y.position', 'x', 'groups', 'x', 'xmin', 'xmax')]
  # write.csv(stat.test_sub, paste(out_dir, analysis_type, '_rad', rad, save_append, '.csv', sep=''))

  bxp <- bxp + stat_pvalue_manual(stat.test, label='p.adj.signif', size=2.5) +
    scale_y_continuous(expand=expansion(mult=c(0, .05))) # add 10% spaces btwn labels and plot border
  
  ggsave(paste(out_dir, analysis_type, '_rad', rad, save_append, '.png', sep=''), height=6, width=8)
  
  # get significant pairs of gcross aucs
  signif_pairs_test <- stat.test |>
    filter(p.adj.signif != 'ns')
  print(signif_pairs_test)
  signif_pairs <- as.character(signif_pairs_test$cell_pair)
  
  signif_df_merged <- df_merged[df_merged$cell_pair %in% signif_pairs, ]
  signif_df_merged$cell_pair <- as.factor(signif_df_merged$cell_pair)
  
  if (length(signif_pairs) != 0) {
    cat(length(signif_pairs), ' signif_pairs found\n')
    
    bxp <- ggboxplot(df_merged, x='cell_pair', y='AUC', color=resp_col) + 
      ylab('AUC ') +
      labs(title=paste(title, ' at rad=',rad, sep='')) +
      theme(text=element_text(size=10),
            axis.text.x = element_text(angle=90, vjust=0.5),
            plot.title = element_text(hjust = 0.5))
    
    ggsave(paste(out_dir, analysis_type, '_rad', rad, '_gcross_sig.png', sep=''), height=6, width=8)
    signif_pairs_test <- signif_pairs_test %>% select(-groups)
    # write.csv(signif_pairs_test, paste0(out_dir, analysis_type, '_rad', rad, '_gcross_sig_stat_test.csv'))
    
  } else {
    cat('no signif_pairs found\n')
  }
  
  return (bxp)
  
}

####################################################################################
## generate plots ##

# create output directory if not present
out_dir <- paste0(out_bp, '1b_auc_plot/')
if (!file.exists(out_dir)) { # create out dir if does not exist
  dir.create(out_dir)
}

# for each distance (in rComp), generate individual AUC plot and save
for (r in 1:length(rComp)) {
  full_df <- data.frame() # contains all the aucs for all cell pairs at the radius r only
  rad <- rComp[r]
  for (i in 1:length(gcross_aucs)) { # loop through the cell pairs
    cell_pair <- names(gcross_aucs)[i]
    if (nrow(gcross_aucs[[cell_pair]]) == 0) { # just in case empty
      next
    }
    cell_pair_df <- data.frame(sub_col = gcross_aucs[[cell_pair]][[sub_col]], # the cases - either images of patients
                     AUC = as.numeric(gcross_aucs[[cell_pair]][[rad]]),
                     cell_pair = cell_pair)
    full_df <- rbind(full_df, cell_pair_df)
  }
  colnames(full_df)[1] <- sub_col # rename this column to be either the image or patient 
  
  # now merge with patient response information
  df_merged <- merge(full_df, pat_resp, by=sub_col) # merge by either image or patient

  # VERY IMPORTANT! - convert cell_pair to factor variable, otherwise p-vals are not in right position
  df_merged$cell_pair <- as.factor(df_merged$cell_pair)
  
  bxp <- plot_bxp(df_merged, resp_col, out_dir, rad, save_append = '')

}

load("/nfs/turbo/umms-ukarvind/shared_data/gcross_data/example/output/1a_run_gcross/gcross_aucs_imwise_standard.RData")

library(tidyverse)
library(ggplot2)

# Set the radius you're interested in
r_selected <- "31"  
sub_col <- "sample"  # or "image" or "ID" — first column name of each df
gcross_aucs <- gcross_aucs  # make sure this is already loaded

# Initialize list to hold AUC rows
auc_dfs <- list()

for (pair in names(gcross_aucs)) {
  df <- gcross_aucs[[pair]]
  
  # Skip if radius column is missing or empty
  if (!(r_selected %in% colnames(df)) || nrow(df) == 0) next
  
  # Extract columns safely
  sample_ids <- df[[sub_col]]
  auc_vals <- as.numeric(df[[r_selected]])
  
  # Check that sample and AUC have the same number of entries
  if (length(sample_ids) != length(auc_vals) || length(sample_ids) == 0) next
  
  this_df <- data.frame(
    sample = sample_ids,
    AUC = auc_vals,
    cell_pair = pair
  )
  
  auc_dfs[[pair]] <- this_df
}

# Combine all valid data.frames
full_auc_df <- bind_rows(auc_dfs)

# Plot raw AUC distributions
ggplot(full_auc_df, aes(x = cell_pair, y = AUC)) +
  geom_boxplot(outlier.shape = 21, fill = "skyblue", color = "black", alpha = 0.7) +
  labs(title = paste("G-cross AUCs at r =", r_selected),
       x = "Cell Type Pair", y = "AUC") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


# Choose any image or case from your data
image_id <- names(gcross_data)[2]  # or specify e.g. "Charville_c001..."
pair <- "Macrophage,Stroma"

# Extract G(r) data
g_data <- gcross_data[[image_id]][[pair]]

r_vals <- g_data$r
obs_vals <- g_data$km
theo_vals <- g_data$theo

# Plot
plot(r_vals, obs_vals,
     type = "l", lwd = 2, col = "blue",
     ylim = range(c(obs_vals, theo_vals), na.rm = TRUE),
     xlab = "Radius (µm)", ylab = "G(r)",
     main = paste("Gcross Curve for", pair, "\nSample:", image_id))

lines(r_vals, theo_vals, col = "red", lty = 2, lwd = 2)
legend("bottomright", legend = c("Observed", "CSR (Theoretical)"),
       col = c("blue", "red"), lty = c(1, 2), bty = "n")

