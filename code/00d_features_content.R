# ==============================================================================
#            ENVIRONMENTAL DATA - FEATURES - CONTENT & STATUS
# ==============================================================================

#' Script by A. Schickele, 2022
#' PLEASE READ THE FOLLOWING INSCRUCTIONS CAREFULLY:
#' 
#' This script adds variable to the features.grd .gri files for prototyping new
#' models and variable selections
#' 
#' To avoid modifying the C4-CCM paper data and Blue-Cloud data, this script can
#' also reverse the process by re-creating the original features.grd .gri files
#' 
#' It is of primary importance to keep a reversible features.grd .gri files
#' because it is called throughout the pipeline, including the database creation
#' via the X0 table.
#' 
#' When changing the features.grd .gri files, open the 00c.build_omics.data.R
#' script and update the X0 table
#' 
#' This is highly prototypical and will be fixed when the C4-CCM paper will be
#' published

# ============================ LOADING THE CODE ================================
# --- Loading the code
bluecloud_dir <- "/home/aschickele/workspace/bluecloud"
data_dir <- paste0(bluecloud_dir, "/data")

setwd(bluecloud_dir)
source("./code/00a_config.R")
library(raster)

# --- Open feature file and check status
# The original features file has 68 layers
features_cur <- raster::stack(paste0(data.wd,"/features"))

if(nlayers(features_cur) == 68){
  message("--- The features file is currently in its original state")
  names(features_cur)
} else {
  message("--- The features file is currently in a modified state")
  names(features_cur)
}

# ======================= REVERSE FEATURES BACK TO ORIGINAL ====================
# Replacing the current features file by the saved features C4 CCM file
# Otherwise, select the first 68 layers

features_original <- raster::stack(paste0(data_dir,"/features_C4_CCM"))
raster::writeRaster(features_original, paste0(data_dir,"/features"), overwrite = TRUE)

# ==================== ADD NEW LAYERS TO FEATURES FOR PROTOTYPING ==============
# New layers are available in the "home/aschickele/workspace/CMEMS product" dir
# The creation of new layers is external to the bluecloud project

supp_dir <- "/home/aschickele/workspace/CMEMS product"
supp_raster <- list.files(supp_dir) %>% 
  grep(pattern = "gri", value = TRUE) %>% 
  substr(1, nchar(.)-4)

# --- Eventually select rasters to add
message("--- Available supplementary rasters are :")
supp_raster

ID <- NULL
supp_raster <- supp_raster[ID]

# --- Add supplementary rasters
for(i in 1:length(supp_raster)){
  message(paste("--- Stacking", supp_raster[i], "to features_cur ---"))
  tmp <- stack(paste0(supp_dir,"/", supp_raster[i]))
  features_cur <- stack(features_cur, tmp)
}

message("--- The features file is currently in a modified state")
names(features_cur)

# --- Save
raster::writeRaster(features_cur, paste0(data_dir,"/features_new"), overwrite = TRUE)

#  Have to write new, delete and rename... because it doesn't want to overwrite
# Acts as a supplementary security step... kind of !
unlink(paste0(data_dir,"/features.gri"))
unlink(paste0(data_dir,"/features.grd"))
file.rename(paste0(data_dir,"/features_new.gri"), paste0(data_dir,"/features.gri"))
file.rename(paste0(data_dir,"/features_new.grd"), paste0(data_dir,"/features.grd"))


# --- END
