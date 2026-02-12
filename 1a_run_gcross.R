####################################################################################
## overview: general script to run gcross on dataset
## created: 4/21/25

## script inputs:
## 1. processed data file: {path_to_data_folder}/subset_data.csv

## script outputs:
## 1. gcross_data: proj_folder/output/1a_run_gcross/gcross_data_{imwise/patwise}.RData
## 2. gcross_aucs: proj_folder/output/1a_run_gcross/gcross_aucs_{imwise/patwise}_{aucobs/ratio}.RData

####################################################################################
## environment setup ##
library(pracma)
library(spatstat)
library(ggplot2)
library(dplyr)
library(tidyverse)
library(yaml)

####################################################################################
## path/folder setup -- set these variables! ##

code_fd <- '/nfs/turbo/umms-ukarvind/gouri/gcross_pkg/'         # path to the gcross_pkg code folder; eg: /nfs/turbo/umms-ukarvind/{your_fd}/gcross_pkg/
proj_fd <- 'example'  # project folder/name; eg: stanford_crc

data_fd <- '/nfs/turbo/umms-ukarvind/shared_data/gcross_data/'         # path to general data folder; eg: /nfs/turbo/umms-ukarvind/shared_data/gcross_data/ 
csv_name <- 'subset_data.csv' # name of csv data file; eg: subset_data.csv
csv_file <- paste0(data_fd, '/', proj_fd, '/processed/', csv_name)  # path to csv file; NOTE: if folder is structured as described in protocol, no need to edit this

source(paste0(code_fd,'main_code/gcross_functions.R')) # NOTE: no need to edit this

####################################################################################
## gcross paramater setup -- set these variables! ##

# define the gcross radius range of computation (gcross will be computed at each of these values) 
# NOTE: a larger vector/range will just take longer to compute
R <- 0:31

# define the radii at which you want to evaluate the gcross AUC (gcross AUC will be computed at each of these values)
# NOTE: values must be smaller than max(R), as defined above
rcomp <- c(0, 31) 

## these parameters most likely do not need to be modified (max_pad, scale_pad, minW):

# amount to pad the xmax and ymax for defining the obs window; 
# NOTE: if Error in `marks<-.ppp`(`*tmp*`, value = factor(unlist(lb))) : number of points != number of marks
# then increase this max_pad value
max_pad <- 500 

scale_pad <- 2  # scale factor to separate out the images within a single patient (only important for patwise)
minW <- -1      # min window x dim; -1 used in NatComm paper; can fix this if points lying outside error

####################################################################################
## load and prepare data ##

tbl <- read.table(csv_file, header=TRUE,sep=",", comment.char = "%", check.names=TRUE) # the '#' character is read as special, so change comment.char to something else

# load the parameter file
params <- yaml.load_file(paste0(code_fd, 'param_files/', proj_fd, '.yaml'))
for (param_name in names(params)) {
  assign(param_name, params[[param_name]])
}

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

img_list <- split(tbl,f=tbl[[im_col]]) # large list, where each entry is the image name, and sublist is image properties (ie. X, Y, etc.)
access_list <- unique(tbl[, c(im_col, pat_col)]) # mapping between the images and patients
cases <- as.character(unique(access_list[, sub_col])) # will either be unique list of images or unique list of patients - must be character!

####################################################################################
## run gcross ##

gcross_data <- list()
cat('running Gcross...\n')

for(i in 1:length(cases)){ # loop through the cases - either unique images or unique patients
  these_ims <- access_list[access_list[[sub_col]]==cases[i], im_col] # get list of images corresponding to the case
  
  xx <- yy <- lb <- ow <- list() # empty lists to hold x-coords, y-coords, phenotype label, obs window
  all_xmax <- all_ymax <- list() # max values across all points in all images
  
  for(j in 1:length(these_ims)){ # loop through the images corresponding to this case
    print(c(i,j))
    
    curr_imname <- these_ims[j] # name of this image
    curr_img <- img_list[[curr_imname]] # sublist of properties of this image
    
    all_xmax[j] <- max(curr_img[[x_col]]) + max_pad # get the max coordinates in image
    all_ymax[j] <- max(curr_img[[y_col]]) + max_pad
    
    xx[[j]] <- curr_img[[x_col]] + (all_xmax[[j]]*scale_pad*(j-1)) # set coordinates of cells in image
    yy[[j]] <- curr_img[[y_col]] + (all_ymax[[j]]*scale_pad*(j-1)) 
    lb[[j]] <- as.character(curr_img[[ct_col]]) # label the cell type
    
    if (length(xx[[j]]) < 3) { # cw will be NULL if less than 3 points in set
      next
    }
    
    cw <- convexhull.xy(xx[[j]], yy[[j]]) # the convex hull of the current set of coordinates is determined; returns window
    ow[[j]] <- cw$bdry[[1]] # set the observation window
  }
  
  xmax <- max(unlist(all_xmax)) # get the max coordinate out of all the images
  ymax <- max(unlist(all_ymax))
  ww <- owin(poly=ow) # creates the observation windows
  
  # convert the x,y coordinates into the form of spatial point patterns
  pp <- as.ppp(cbind(unlist(xx), unlist(yy)),
               W = c(minW, xmax*2*length(these_ims), minW, ymax*2*length(these_ims))) # As used in NatComm paper
  
  marks(pp) <- factor(unlist(lb)) # label the points as phenotypes, as a factor
  plot(pp) # plot for visualization; check clusters are separated for images within a patient
  pp$window <- ww # set the pp window to be based on the set of points observed
  
  gcross <- alltypes(pp, fun='Gcross', r=R) # run gcross on each combo of cell types; can also just replace this with Kcross if want to
  gfuns <- gcross$fns
  
  ctypes <- as.character(unique(colnames(gcross$which))) # list of phenotypes present in this image
  combs <- expand.grid(ctypes, ctypes, stringsAsFactors=FALSE) # all the pairwise interactions 
  cross_nlist <- (interaction(combs$Var2, combs$Var1, sep=",")) # pairwise interactions converted to single factor in form of i, j where i=center, j=other
  
  names(gfuns) <- cross_nlist # list of gfuns is pairwise interactions
  gcross_data[[cases[i]]] <- gfuns # add to gcross_data list, where sublist is name of case - either im or pat name
  
}

out_dir <- paste(data_fd, proj_fd, '/output/1a_run_gcross/', sep='')
if (!file.exists(out_dir)) { # create out dir if does not exist
  dir.create(out_dir)
}

save(gcross_data, file=paste(out_dir, "gcross_data_", analysis_type, ".RData", sep=''))
cat('Gcross complete & output saved!\n')

####################################################################################
## compute gcross AUCs ##

# get all the cell types
cell_types <- unique(tbl[[ct_col]])
combs <- expand.grid(cell_types, cell_types, stringsAsFactors=FALSE) # all the pairwise interactions 

if (length(cell_pairs) == 0) { # if no interested cell pairs provided, list all pairs
  combs2 <- combs |> 
    dplyr::filter(Var1 != '') |>
    dplyr::filter(Var2 != '') |>
    dplyr::filter(Var1 != Var2) |> # also filter out same relationships
    dplyr::filter(Var1 != 'Other') |>
    dplyr::filter(Var2 != 'Other')
  cell_pairs <- as.character(interaction(combs2$Var2, combs2$Var1, sep=","))
}

# if there were spaces in cell type names, this will fix it to match the gcross data
renamed_list <- lapply(cell_pairs, function(x) {
   gsub(" ", ".", x) })
cell_pairs <- unlist(renamed_list)

## NOTE: this will automatically compute both un-normalized G-Cross AUCs (_standard)
## and normalized G-Cross AUCs (_ratio)
for (j in c('_standard', '_ratio')){ # compute un-normalized and normalized aucs
  gcross_aucs <- list()
  for (k in 1:length(cell_pairs)) { # loop through all cell pairs
    cell_pair <- as.character(cell_pairs[k]) # current cell pair
    
    # create matrix to hold aucs, for each case (im or pat), holds the auc computed at each rComp; indexed by cell_pair
    aucs <- matrix(nrow=length(gcross_data), ncol=length(rcomp)+1)
    colnames(aucs) <- c(sub_col, as.character(rcomp)) # c(im_col or pat_col, rComp)
    
    for (m in 1:length(gcross_data)) { # loop through each case - images or patients
      if (!(cell_pair %in% names(gcross_data[[m]]))) { # if pair not present, then skip
        next
      }
      
      id <- names(gcross_data)[m] # this will either be patient or image id
      this_gcross <- gcross_data[[id]][[cell_pair]] # this case + this cell pair
      
      # corresponding gcross obs and theo values
      gcross_obs <- this_gcross$km # change to iso for kcross, km for gcross; change to bord for inhom
      gcross_theo <- this_gcross$theo
      
      these_aucs <- numeric() # list of aucs for this cell pair
      for (rad in rcomp){ # compute the aucs at each distance of interest
        auc_obs <- trapz(0:rad, (gcross_obs[1:(rad+1)]))
        auc_theo <- trapz(0:rad, (gcross_theo[1:(rad+1)]))
        
        if (j == '_aucobs') {
          these_aucs <- append(these_aucs, auc_obs)
        } else {
          # auc needs to be the ratio between the observed and theoretical - represents on average deviation
          these_aucs <- append(these_aucs, auc_obs/auc_theo)
        }
      }
      
      aucs[m,] <- append(id, these_aucs) # append case id and aucs for each rComp
    }
    
    # save the aucs computed for this case into the aucs matrix
    aucs <- as.data.frame(aucs)
    aucs <- na.omit(aucs)
    gcross_aucs[[cell_pair]] <- aucs
  }
  
  save(gcross_aucs, file=paste(out_dir, "gcross_aucs_", analysis_type, j, ".RData", sep=''))
  cat('Gcross AUCs complete & output saved!\n')
}

# load data
load("/nfs/turbo/umms-ukarvind/shared_data/gcross_data/example/output/1a_run_gcross/gcross_aucs_imwise_standard.RData")
ls()  # Lists all variables in your environment
dim(aucs)
