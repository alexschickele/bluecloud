#' TO DO:
#' - check conflict with the "DEPTH" parameters

ls()
rm(list=ls())

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

# --- Input / Output directories
data.wd <- "/home/aschickele/complex"
bluecloud.wd <- "/home/aschickele/workspace/bluecloud descriptor"

# --- Custom functions
source_python(paste0(bluecloud.wd,"/function/mbtr_function.py"))
source(file = paste0(bluecloud.wd,"/function/bivarRasterPlot.R"))

# --- Other custom arguments
Sys.setenv(HDF5_USE_FILE_LOCKING="FALSE") # to be able to open .nc from complex

# --- Data specific parameters
DEPTH <- data.frame(top=0, bottom=30) #selected depth for WOA data

MIN.GENE <- 3 # minimum of genes per cluster represented in the data summary
MIN.STATION <- 100 # minimum presence of genes in the data summary

DEPTH <- "SUR" # selected depth for genomic data
FILTER <- "QQSS" # selected size filter for genomic data
CLUSTER <- "CC_842178" # selected protein family

# --- Model specific parameters
# Hyperparameters to test in the model (03 and 04)
HYPERPARAMETERS <- data.frame(LEARNING_RATE = c(6e-2, 3e-2, 1e-2, 6e-3, 3e-3, 1e-3),
                              N_Q = c(5, 10, 20, 33, 50, 100),
                              MEAN_LEAF = c(3, 6, 10, 15, 20, 25))

NBOOST <- 1000 # maximum number of boosting rounds
N_FOLD <- 3 # number of k-fold cross validation runs
MAX_CLUSTER <- 36 # maximum number of clusters for parallel computing

NBOOTSTRAP <- 5 # number of bootstrap rounds for script 05b

# --- END
