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

MAX_CLUSTER <- 10

# =========================== DEFINE PARAMETERS ================================
kegg_p0 = c("00190","00195","00196","00710","00680","00910","00920")

for(i in 1:length(kegg_p0)){
  kegg_p = kegg_p0[i]
  cc_id = NULL
  cluster_selec = c(20,2,25)
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
  
  # 1. Query data & save ---------------------------------------------------------
  query <- query_data(bluecloud.wd = bluecloud_dir,
                      CC_id = cc_id,
                      KEGG_p = kegg_p,
                      CLUSTER_SELEC = list(N_CLUSTERS = cluster_selec[1], MIN_GENES = cluster_selec[2], MAX_GENES = cluster_selec[3]),
                      ENV_METRIC = c(env_metric, "dist","bathy"),
                      relative = relative)
  
  file.copy(from = paste0(data_dir, "/X.feather"), to = paste0(bluecloud_dir, "/output/", output_dir))
  file.copy(from = paste0(data_dir, "/Y.feather"), to = paste0(bluecloud_dir, "/output/", output_dir))
  file.copy(from = paste0(data_dir, "/CC_desc.feather"), to = paste0(bluecloud_dir, "/output/", output_dir))
  
  # 2. Run model & save ----------------------------------------------------------
  run <- model_run(bluecloud.wd = bluecloud_dir,
                   HYPERPARAMETERS = data.frame(LEARNING_RATE = c(1e-3, 1e-3),
                                                N_Q = c(10, 10),
                                                MEAN_LEAF = c(25, 50)),
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
  map_proj(proj = proj$proj, col = proj$col, targetID = seq(1:nlayers(proj$proj)))
  cor_proj(y_hat_m = proj$y_hat_m)
  dev.off()
} # global run loop



# --- END