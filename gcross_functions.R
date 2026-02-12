###
## overview: helper functions for running gcross
## created: 10/2024
##

### convert local image coordinates to global scan coordinates
convert_to_global_coords <- function(tbl, im_col, x_col, y_col, stamp_coord_type, xy_coord_type, stamp_coord_loc) {
  
  convert_to_microns <- function(value, conversion_factor) {  # helper function for conversion
    return(value * conversion_factor)  
  }
  pixel_to_micron <- 0.5 # 1 pixel = 0.5 microns
  
  # Initialize vectors to store the modified coordinates
  global_x <- numeric(nrow(tbl))
  global_y <- numeric(nrow(tbl))
  
  # Iterate over each unique image in the table
  unique_images <- unique(tbl[[im_col]])
  
  for (img in unique_images) {
    
    # Subset the table for the current image
    img_subset <- subset(tbl, tbl[[im_col]] == img)
    
    # Extract stamp coordinates from the image filename (e.g., 29_Scan1_[11836,54920])
    # Use regular expression to extract the values within the square brackets
    stamp_coords <- as.numeric(unlist(strsplit(img, "[^0-9]+"))[3:4])
    
    # Convert stamp coordinates to microns if the stamp_coord_type is in pixels
    if (stamp_coord_type == "pixels") {
      stamp_coords <- convert_to_microns(stamp_coords, pixel_to_micron)
    }
    
    # Convert local xy coordinates to microns if they are in pixels
    if (xy_coord_type == "pixels") {
      img_subset[[x_col]] <- convert_to_microns(img_subset[[x_col]], pixel_to_micron)
      img_subset[[y_col]] <- convert_to_microns(img_subset[[y_col]], pixel_to_micron)
    }
    
    # Handle different reference locations
    if (stamp_coord_loc == 'TL') {  # Top-Left
      global_x[tbl[[im_col]] == img] <- img_subset[[x_col]] + stamp_coords[1]
      global_y[tbl[[im_col]] == img] <- img_subset[[y_col]] + stamp_coords[2]
    }
    
    if (stamp_coord_loc == 'TR') {  # Top-Right
      global_x[tbl[[im_col]] == img] <- stamp_coords[1] - img_subset[[x_col]]
      global_y[tbl[[im_col]] == img] <- img_subset[[y_col]] + stamp_coords[2]
    }
    
    if (stamp_coord_loc == 'C') {  # Center
      global_x[tbl[[im_col]] == img] <- img_subset[[x_col]] + (stamp_coords[1] - img_subset[[x_col]] / 2)
      global_y[tbl[[im_col]] == img] <- img_subset[[y_col]] + (stamp_coords[2] - img_subset[[y_col]] / 2)
    }
  }
  
  # Return a dataframe with column names as x_col and y_col
  global_coords <- data.frame(global_x = global_x, global_y = global_y)
  names(global_coords) <- c(x_col, y_col)  # Rename columns based on input strings
  
  return(global_coords)
}
