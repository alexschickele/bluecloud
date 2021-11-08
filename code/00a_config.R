ls()
rm(list = setdiff(ls(), lsf.str()))

# --- R Package
# used in step 1 to 3 = background dependencies of shiny
library(raster) # in shiny
library(abind) # in shiny
library(feather) # in shiny
library(reticulate) # in shiny
library(RColorBrewer) # in shiny
library(parallel) # in shiny
library(mvrsquared) #in shiny
library(tidyverse) # in shiny
library(RSQLite) # in shiny

# used only in step 0
library(ncdf4)
library(virtualspecies)
library(oce)
library(castr)
library(vroom)

# --- Input / Output directories
data.wd <- "/home/aschickele/complex"
bluecloud.wd <- "/home/aschickele/workspace/bluecloud descriptor"

# --- Custom functions
source_python(paste0(bluecloud.wd,"/function/mbtr_function.py"))
source(file = paste0(bluecloud.wd,"/function/bivarRasterPlot.R"))

# --- Other custom arguments
Sys.setenv(HDF5_USE_FILE_LOCKING="FALSE") # to be able to open .nc from complex

# --- Data specific parameters
DEPTH <- c("SUR") # selected depth for genomic data
FILTER <- c("GGMM") # selected size filter for genomic data

# CLUSTER_SELEC <- list(MIN_STATIONS = 80, MIN_GENES = 5, MAX_GENES = 25) # for cluster selection
# KEGG_p <-  "00190" # for cluster selection
# ENV_METRIC <- c("mean","sd","med","mad","dist","bathy")

# --- Model specific parameters
# Hyperparameters to test in the model (03 and 04)
# HYPERPARAMETERS <- data.frame(LEARNING_RATE = c(1e-2, 1e-2, 1e-2, 1e-2),
#                               N_Q = c(5, 10, 20, 50),
#                               MEAN_LEAF = c(30, 40, 50, 60))

NBOOST <- 30 # maximum number of boosting rounds
N_FOLD <- 5 # number of k-fold cross validation runs
MAX_CLUSTER <- 20 # maximum number of clusters for parallel computing

NBOOTSTRAP <- 20 # number of bootstrap rounds for script 05b
