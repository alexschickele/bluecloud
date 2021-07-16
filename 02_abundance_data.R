
#' Here we build the relative abundance dataset with the TARA OCEAN data
#' The output are
#' 1. A relative abundance datasset : Ytargets * Nobs
#' 2. The corresponding environmental features : Xfeatures * Nobs

input.wd <- "~/complex/data"
output.wd <- "~/workspace/bluecloud descriptor"

# --- Loading R packages
library(raster)
library(feather)

# --- Input parameters


# --- Loading data
lonlat <- read.csv(paste0(input.wd,"/SMAGs_Env.csv"), sep=';', header = TRUE)
relabs <- read.table(paste0(input.wd,"/FF_metaG_Unknown_annot_subset.txt"))
