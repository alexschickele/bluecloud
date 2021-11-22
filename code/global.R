message("load global functions")

library(shiny)
library(shinybusy)

source("00a_config.R")
source("01_query_data.R")
source("02a_model_param.R")
source("02b_model_eval.R")
source("03a_bootstrap_predict.R")

bluecloud_dir <- "/home/aschickele/workspace/bluecloud descriptor"
data_dir <- "/home/aschickele/workspace/bluecloud descriptor/data"