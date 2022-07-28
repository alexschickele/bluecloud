message("load global functions")

library(shiny)
library(shinybusy)
library(shinyjs)

bluecloud_dir <- "/home/aschickele/workspace/bluecloud"
data_dir <- "/home/aschickele/workspace/bluecloud/data"

source("00a_config.R")
source("03a_bootstrap_predict.R") # for plot functions

# --- Loading data tables
load(paste0(data_dir,"/shiny_data.RData"))
CC_desc_e <- query$CC_desc[query$e$vr,] %>% inner_join(query$nn_ca)

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

