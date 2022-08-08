message("load global functions")

library(shiny)
library(shinyjs)
library(shinythemes)
library(leaflet)
library(raster)
library(tidyverse)

bluecloud_dir <- bluecloud.dir <- "/home/aschickele/workspace/bluecloud/app"

# --- Loading data tables and functions
setwd(bluecloud_dir)
source("./app_function.R")
load("./app_data.RData")

# --- Defining the enzyme - KO correspondance
plot_list <- list(RUBISCO = "01601|01602",
                  PEPC = "1595",
                  GOT = "14454|14455",
                  PEPCK = "01610",
                  MDH_NAD = "00024|00025|00026",
                  MDH_NADP = "00051",
                  MDC_NADP = "00029",
                  MDC_NAD = "00028",
                  GPT_GGAT = "00814|14272",
                  PEPDK = "1006")

# right CRS for leaflet
crs(proj$proj) <- sp::CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")

