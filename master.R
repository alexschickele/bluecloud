#' Master script for running the R pipeline, open and save results as well as
#' performing global analysis on several model runs

# ============================ LOADING THE CODE ================================
# --- Loading the code
bluecloud_dir <- "/home/aschickele/workspace/bluecloud"
data_dir <- paste0(bluecloud_dir, "/data")

setwd(bluecloud_dir)
source("./code/00a_config.R")
source("./code/01_query_data.R")
source("./code/02a_model_param.R")
source("./code/02b_model_eval.R")
source("./code/03a_bootstrap_predict.R")

MAX_CLUSTER <- 20

# =========================== DEFINE PARAMETERS ================================
kegg_p0 = c("00190","00195","00710","00680","00910","00920")
kegg_m0 = list(paste0("M00",c(144:160,416,417)),
               paste0("M00",c(157,161:163)),
               paste0("M00",165:172),
               paste0("M00",c(567,357,356,563,358,608,174,346,345,344,378,935,422)),
               paste0("M00",c(175,531,530,529,528,804)),
               paste0("M00",c(176,596,595,"021")))

for(i in 1:length(kegg_p0)){
  kegg_p = kegg_p0[i]
  kegg_m = kegg_m0[[i]]
  cc_id = NULL
  cluster_selec = c(30,2,30)
  env_metric = c("mean", "sd")
  relative = TRUE
  
  output_dir <- paste0(kegg_p, "_", cc_id, "_", 
                       paste(cluster_selec, collapse = "_"), "_",
                       paste(env_metric, collapse = "_"), "_",
                       if(relative==TRUE){"rel"} else {"abs"})
  
  # ================= LOAD IF MODEL EXISTS // CREATE DIR IF NOT ==================
  if (length(list.files(paste0(bluecloud_dir, "/output/", output_dir)))!=0){
    message(">>> This model has already been run ... loading files <<<")
    old_files <- list.files(paste0(bluecloud_dir, "/output/", output_dir))
    old_files <- old_files[-grep("pdf|RData", old_files)]
    for(i in old_files){
      file.copy(from = paste0(bluecloud_dir, "/output/", output_dir, "/", i), to = data_dir, overwrite = TRUE)
    }
    load(file = paste0(bluecloud_dir, "/output/", output_dir, "/output.RData"))
  } else {
    message(">>> This model has never been run ... create directory <<<")
    dir.create(paste0(bluecloud_dir, "/output/", output_dir))
  }
  
  # ===================== RUN THE PIPELINE AND SAVE ==============================
  # With the data present in /data if partially run
  
  # 1. Query data, synthetic cluster plot & save ---------------------------------
  query <- query_data(bluecloud.wd = bluecloud_dir,
                      CC_id = cc_id,
                      KEGG_m = kegg_m,
                      CLUSTER_SELEC = list(N_CLUSTERS = cluster_selec[1], MIN_GENES = cluster_selec[2], MAX_GENES = cluster_selec[3]),
                      ENV_METRIC = c(env_metric, "dist","bathy"),
                      relative = TRUE)
  
  file.copy(from = paste0(data_dir, "/X.feather"), to = paste0(bluecloud_dir, "/output/", output_dir))
  file.copy(from = paste0(data_dir, "/Y.feather"), to = paste0(bluecloud_dir, "/output/", output_dir))
  file.copy(from = paste0(data_dir, "/CC_desc.feather"), to = paste0(bluecloud_dir, "/output/", output_dir))
  
  png(filename = paste0(bluecloud_dir, "/output/", output_dir, "/CC_module_plot.png"))
  image(query$CC_module, col = "black", axes = F, main = "Which KEGG module is related to which CC ?")
  axis(side = 1, at = seq(0,1,length.out = nrow(query$CC_module)), labels = rownames(query$CC_module), las = 2, cex.axis = 0.6)
  axis(side = 2, at = seq(0,1,length.out = ncol(query$CC_module)), labels = colnames(query$CC_module), las = 2, cex.axis = 0.6)
  abline(v = seq(0,1,length.out = nrow(query$CC_module)), h  = seq(0,1,length.out = ncol(query$CC_module)), lty = "dotted")
  box()
  dev.off()
  
  # 2. Run model & save ----------------------------------------------------------
  run <- model_run(bluecloud.wd = bluecloud_dir,
                   HYPERPARAMETERS = data.frame(LEARNING_RATE = c(1e-3, 1e-3, 1e-3, 1e-3),
                                                N_Q = c(10, 10, 10, 10),
                                                MEAN_LEAF = c(20, 30, 40, 50)),
                   verbose = TRUE)
  
  file.copy(from = paste0(data_dir, "/HYPERPARAMETERS.feather"), to = paste0(bluecloud_dir, "/output/", output_dir), overwrite = TRUE)
  file.copy(from = paste0(data_dir, "/Station_FOLD.feather"), to = paste0(bluecloud_dir, "/output/", output_dir), overwrite = TRUE)
  file.copy(from = paste0(data_dir, "/Station_ID.feather"), to = paste0(bluecloud_dir, "/output/", output_dir), overwrite = TRUE)
  file.copy(from = paste0(data_dir, "/m"), to = paste0(bluecloud_dir, "/output/", output_dir), overwrite = TRUE)
  for(n in 1:N_FOLD){
    file.copy(from = paste0(data_dir, "/", n,"_X_val.feather"), to = paste0(bluecloud_dir, "/output/", output_dir), overwrite = TRUE)
  }
  
  # 3. Evaluate model ------------------------------------------------------------
  eval <- model_eval(bluecloud.wd = bluecloud_dir,
                     by_target = TRUE,
                     var_importance = FALSE)
  
  # 4. Do projections & save -----------------------------------------------------
  proj <- model_proj(bluecloud.wd = bluecloud_dir,
                     data.wd = data_dir,
                     ENV_METRIC = c(env_metric, "dist","bathy"))
  
  save(query, eval, proj, file = paste0(bluecloud_dir, "/output/", output_dir, "/output.RData"))
  
  # 5. Do plot outputs & save ----------------------------------------------------
  pdf(paste0(bluecloud_dir, "/output/", output_dir, "/", output_dir, ".pdf"))
  par(mfrow = c(4,2))
  legend_proj(col_matrix = proj$col_matrix, cutx = proj$cutx, cuty = proj$cuty)
  map_proj(proj = proj$proj, 
           col = proj$col, 
           targetID = seq(1:nlayers(proj$proj))[order(match(query$CC_desc$CC,rownames(query$CC_module)))], 
           targetNAME = query$CC_desc$CC)
  cor_proj(y_hat_m = proj$y_hat_m[,order(match(query$CC_desc$CC,rownames(query$CC_module)))])
  dev.off()
} # global run loop



# --- END