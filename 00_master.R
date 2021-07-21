ls()
rm(list=ls())

setwd("~/workspace/bluecloud descriptor")

source(file = "01_environmental_data.R")

env_data(input.wd="~/complex/share/WOA/DATA",
         output.wd="~/workspace/bluecloud descriptor")


