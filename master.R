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
# kegg_p0 = c("00190","00195","00710","00680","00910","00920")
# kegg_m0 = list(paste0("M00",c(144:160,416,417)),
#                paste0("M00",c(157,161:163)),
#                paste0("M00",165:172),
#                paste0("M00",c(567,357,356,563,358,608,174,346,345,344,378,935,422)),
#                paste0("M00",c(175,531,530,529,528,804)),
#                paste0("M00",c(176,596,595,"021")))

kegg_p0 = c("00710","00710","00710","00710","00710","00710")
kegg_m0 = list(paste0("M00",165:172),
               paste0("M00",165:172),
               paste0("M00",165:172),
               paste0("M00",165:172),
               paste0("M00",165:172),
               paste0("M00",165:172))
cluster_selec0 = list(c(30,5,0),
                     c(30,5,0.5),
                     c(30,5,1),
                     c(30,3,0),
                     c(30,3,0.5),
                     c(30,3,1))

for(p in 1:length(kegg_p0)){
  kegg_p = kegg_p0[p]
  kegg_m = kegg_m0[[p]]
  cc_id = NULL
  cluster_selec = cluster_selec0[[p]]
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
                      CLUSTER_SELEC = list(N_CLUSTERS = cluster_selec[1], MIN_GENES = cluster_selec[2], EXCLUSIVITY_R = cluster_selec[3]),
                      ENV_METRIC = c(env_metric, "dist","bathy"),
                      relative = TRUE)
  
  file.copy(from = paste0(data_dir, "/X.feather"), to = paste0(bluecloud_dir, "/output/", output_dir))
  file.copy(from = paste0(data_dir, "/Y.feather"), to = paste0(bluecloud_dir, "/output/", output_dir))
  file.copy(from = paste0(data_dir, "/CC_desc.feather"), to = paste0(bluecloud_dir, "/output/", output_dir))
  
  # Relation between KEGG_modules and the selected CC
  png(filename = paste0(bluecloud_dir, "/output/", output_dir, "/CC_module_plot.png"))
  image(query$CC_module, col = "black", axes = F, main = "Which KEGG module is related to which CC ?")
  axis(side = 1, at = seq(0,1,length.out = nrow(query$CC_module)), labels = rownames(query$CC_module), las = 2, cex.axis = 0.6)
  axis(side = 2, at = seq(0,1,length.out = ncol(query$CC_module)), labels = colnames(query$CC_module), las = 2, cex.axis = 0.6)
  abline(v = seq(0,1,length.out = nrow(query$CC_module)), h  = seq(0,1,length.out = ncol(query$CC_module)), lty = "dotted")
  box()
  dev.off()
  
  # PCA highlighting the selected CC among ubiquitous, abundant, exclusive etc...
  pca_res <- PCA(X = query$CC_desc[,c(2,8:13)], graph = FALSE)
  ind_col <- rep("black", nrow(query$CC_desc))
  ind_col[query$e$vr[1:cluster_selec[1]]] <- "red"
  plot.PCA(x = pca_res, choix = "ind", habillage = "ind", col.hab = ind_col)
  plot.PCA(x = pca_res, choix = "var")
  
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
  # --- New version with clustering
  y_hat_m_rescaled <- apply(proj$y_hat_m, 2, function(x){x/max(x, na.rm = TRUE)})
  colnames(y_hat_m_rescaled) <- query$CC_desc$CC[query$e$vr[1:cluster_selec[1]]]
  tree <- dist(t(y_hat_m_rescaled)) %>% 
    hclust(method = "ward.D2")
  group <- cutree(tree, k = 5)

  pdf(paste0(bluecloud_dir, "/output/", output_dir, "/", output_dir, ".pdf"))
  par(mfrow = c(2,1),  mar = c(3,11,3,11))
  barplot(rev(tree$height), ylab = "Euclidiian distance", xlab = "Number of clusters", main = "Connected Component (CC) cutoff", col = "black")
  abline(h = tree$height[length(tree$height)-4], col = "red")
  plot(tree, cex = 0.6, ylab= "Euclidian distance", main = "Connected Component (CC) dendrogram")
  abline(h = tree$height[length(tree$height)-4], col = "red")
  cor_proj(y_hat_m = proj$y_hat_m[,tree$order],
           targetNAME = tree$labels[tree$order])
  legend_proj(col_matrix = proj$col_matrix)
  
  par(mfrow = c(4,2), mar=c(2,2,2,0))
  for(n in 1:cluster_selec[1]){
    map_proj(proj = proj$proj, 
             col = proj$col, 
             targetID = tree$order[n], 
             targetNAME = tree$labels)
    map_scale <- max(proj$y_hat_m[,tree$order[n]], na.rm = TRUE)
    barplot(c(map_scale, rep(0, 9)), ylim = c(0,1), border = NA, col = "black",
            main = "Description :")
    abline(v = c(0.2,1.2))
    title(main = "Scale :", adj = 0)
    text(x = rep(2,6), y = c(0.8,0.6,0.5,0.4,0.3,0.1), cex = 0.8, adj = c(0,1),
         labels = c(paste("H. Clustering group:", group[tree$order][n]),
                    paste("N. Modules:", query$CC_desc$max_mod[query$e$vr[1:cluster_selec[1]]][tree$order[n]]),
                    paste("Modules:", query$CC_desc$kegg_module[query$e$vr[1:cluster_selec[1]]][tree$order[n]]),
                    paste("KOs:", query$CC_des$kegg_ko[query$e$vr[1:cluster_selec[1]]][tree$order[n]]),
                    paste("Unknown rate:", query$CC_desc$unknown_rate[query$e$vr[1:cluster_selec[1]]][tree$order[n]]),
                    paste("Class:", query$CC_desc$class[query$e$vr[1:cluster_selec[1]]][tree$order[n]])))
  }
  dev.off()
  
  
} # global run loop



# --- END