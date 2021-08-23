
#' TO DO LIST:
#' - implement nice maps with uncertainty in HCL color figures
#' i.e. vary the luminosity in function of the standard deviation
#' - print in PDF ?

ls()
rm(list=ls())

input.wd <- "~/workspace/bluecloud descriptor"
output.wd <- "~/workspace/bluecloud descriptor"

# --- Loading R packages
library(reticulate)
library(feather)
library(raster)
library(virtualspecies)

