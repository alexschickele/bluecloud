#' TO DO:

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
library(tidyverse)
library(oce)
library(castr)

# --- Input / Output directories
data.wd <- "/home/aschickele/complex"
bluecloud.wd <- "/home/aschickele/workspace/bluecloud descriptor"

# --- Custom functions
source_python(paste0(bluecloud.wd,"/function/mbtr_function.py"))
source(file = paste0(bluecloud.wd,"/function/bivarRasterPlot.R"))

# --- Other custom arguments
Sys.setenv(HDF5_USE_FILE_LOCKING="FALSE") # to be able to open .nc from complex

# --- Data specific parameters
MIN.GENE <- 3 # minimum of genes per cluster represented in the data summary
MIN.STATION <- 10 # minimum presence of genes in the data summary

DEPTH <- c("SUR") # selected depth for genomic data
FILTER <- c("QQSS") # selected size filter for genomic data
CLUSTER <- "CC_995588" # selected protein family

# --- Model specific parameters
# Hyperparameters to test in the model (03 and 04)
HYPERPARAMETERS <- data.frame(LEARNING_RATE = c(1e-4, 1e-4, 3e-4, 1e-4),
                              N_Q = c(10, 20, 50, 100),
                              MEAN_LEAF = c(10,20,30,50))

NBOOST <- 10000 # maximum number of boosting rounds
N_FOLD <- 3 # number of k-fold cross validation runs
MAX_CLUSTER <- 36 # maximum number of clusters for parallel computing

NBOOTSTRAP <- 5 # number of bootstrap rounds for script 05b

# --- END
