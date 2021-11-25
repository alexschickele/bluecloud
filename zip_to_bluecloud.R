
library(rmd2jupyter)
rmd2jupyter("~/workspace/bluecloud/Service_1.Rmd")

zip(zipfile = "~/workspace/bluecloud/bluecloud_wp3_d2_s1.zip",
    files = c("Service_1.ipynb",
              "./function",
              "./code/00a_config.R",
              "./code/01_query_data.R",
              "./code/02a_model_param.R",
              "./code/02b_model_eval.R",
              "./code/03a_bootstrap_predict.R"))
