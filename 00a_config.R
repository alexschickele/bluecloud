#' TO DO:
config <- function(data.wd = "/home/aschickele/complex",
                   bluecloud.wd = "/home/aschickele/workspace/bluecloud descriptor",
                   CLUSTER_SELEC = list(MIN_STATIONS = 80, MIN_GENES = 5, MAX_GENES = 25, KEGG_p = "00190"),
                   ENV_METRIC = c("mean","sd","med","mad","dist","bathy")){

ls()
rm(list = setdiff(ls(), lsf.str()))

# --- R Package
library(ncdf4)
library(raster)
library(abind)
library(virtualspecies)
library(feather)
library(reticulate)
library(RColorBrewer)
library(parallel)
library(mvrsquared)
library(tidyverse)
library(oce)
library(castr)
library(vroom)
library(RSQLite)

# --- Input / Output directories
data.wd <<- data.wd
bluecloud.wd <<- bluecloud.wd

# --- Custom functions
source_python(paste0(bluecloud.wd,"/function/mbtr_function.py"))
source(file = paste0(bluecloud.wd,"/function/bivarRasterPlot.R"))

# --- Other custom arguments
Sys.setenv(HDF5_USE_FILE_LOCKING="FALSE") # to be able to open .nc from complex

# --- Data specific parameters
DEPTH <<- c("SUR") # selected depth for genomic data
FILTER <<- c("GGMM") # selected size filter for genomic data

CLUSTER_SELEC <<- CLUSTER_SELEC # for cluster selection
ENV_METRIC <<- ENV_METRIC

# --- Model specific parameters
# Hyperparameters to test in the model (03 and 04)
HYPERPARAMETERS <<- data.frame(LEARNING_RATE = c(1e-2, 1e-2, 1e-2, 1e-2),
                              N_Q = c(5, 10, 20, 50),
                              MEAN_LEAF = c(30, 40, 50, 60))

NBOOST <<- 3000 # maximum number of boosting rounds
N_FOLD <<- 5 # number of k-fold cross validation runs
MAX_CLUSTER <<- 20 # maximum number of clusters for parallel computing

NBOOTSTRAP <<- 20 # number of bootstrap rounds for script 05b

} # end function
