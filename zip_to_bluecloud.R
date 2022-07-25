
library(rmd2jupyter)
rmd2jupyter("~/workspace/bluecloud/service_1.Rmd")

zip(zipfile = "~/workspace/bluecloud/bluecloud_wp3_d2_s1.zip",
    files = c("service_1.ipynb",
              "./data/X.feather",
              "./data/Y.feather",
              "./function",
              "./code/00a_config.R",
              "./code/01a_query_BC.R",
              "./code/02a_model_param.R",
              "./code/02b_model_eval.R",
              "./code/03a_bootstrap_predict.R"))
