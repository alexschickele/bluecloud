
# --- R Package
# Used in service 1 (R pipeline) and service 2 (R Shiny application)
library(raster) 
library(virtualspecies) 
library(abind) 
library(feather) 
library(reticulate) 
library(RColorBrewer) 
library(parallel) 
library(mvrsquared) 
library(tidyverse) 
library(RSQLite) 
library(RPostgreSQL) 
library(pastecs) 

library(FactoMineR)
library(car)
library(ade4)

# Used only in environmental and genomic data building
# library(ncdf4)
# library(oce)
# library(castr)
# library(vroom)

# --- Input / Output directories
bluecloud.wd <- bluecloud_dir
data.wd <- data_dir

# --- Custom functions
source_python(paste0(bluecloud.wd,"/function/mbtr_function.py"))
source(paste0(bluecloud.wd,"/function/bivarRasterPlot.R"))

# --- Other custom arguments
Sys.setenv(HDF5_USE_FILE_LOCKING="FALSE") # to be able to open .nc from complex

# --- Data specific parameters
DEPTH <- c("SUR") # selected depth for genomic data
FILTER <- c("GGMM","GGZZ") # selected size filter for genomic data

# --- Model specific parameters
NBOOST <- 5000 # maximum number of boosting rounds
N_FOLD <- 10 # number of cross validation runs
MAX_CLUSTER <- 16 # maximum number of clusters for parallel computing

NBOOTSTRAP <- 100 # number of bootstrap rounds for script 03a
