
#' Master function for the BLUECLOUD project descriptor 2 modelling procedure
#' This script successively calls R and Python scripts to
#' - prepare the environmental and relative abundance data
#' - run hyperparameter search for the mbtr model under mse loss
#' - evaluate the selected model
#' - project the results on a global map and corresponding figures
#' 
#' TO DO LIST :
#' - put all the scripts into functions and call them in the master script

ls()
rm(list=ls())

setwd("~/workspace/bluecloud descriptor")

source(file = "01_environmental_data.R")

env_data(input.wd="~/complex/share/WOA/DATA",
         output.wd="~/workspace/bluecloud descriptor")


