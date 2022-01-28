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
  tree_cut <- 5
  group <- cutree(tree, k = tree_cut)

  pdf(paste0(bluecloud_dir, "/output/", output_dir, "/", output_dir, ".pdf"))
  par(mfrow = c(2,1),  mar = c(3,11,3,11))
  barplot(rev(tree$height), ylab = "Euclidian distance", xlab = "Number of clusters", main = "Connected Component (CC) cutoff",
          col = c(rep("black", tree_cut), rep("gray75", cluster_selec[1]-tree_cut)))
  abline(h = tree$height[length(tree$height)-tree_cut], col = "red")
  plot(tree, cex = 0.6, ylab= "Euclidian distance", main = "Connected Component (CC) dendrogram")
  abline(h = tree$height[length(tree$height)-(tree_cut)], col = "red")
  cor_proj(y_hat_m = proj$y_hat_m[,tree$order], targetNAME = tree$labels[tree$order])
  abline(h = match(1:tree_cut, group[tree$order])-0.5, v = match(1:tree_cut, group[tree$order])-0.5)
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
  
  # 6. Multivariate analysis on map clusters for interpretation ----------------
  # Create table
  factor_raw <- list(query$CC_desc$kegg_module[query$e$vr[1:cluster_selec[1]]],
                  query$CC_desc$kegg_ko[query$e$vr[1:cluster_selec[1]]],
                  query$CC_desc$class[query$e$vr[1:cluster_selec[1]]])
  factor_names <- lapply(factor_raw, function(x){x <- gsub(x, pattern = " ", replacement = "") %>% 
                                                       strsplit(split = ",") %>% 
                                                       unlist() %>% unique()
                                                  x <- x[-which(x == "-" | x == "NA")]})
  df <- data.frame(group = as.factor(group), unknown_rate = query$CC_desc$unknown_rate[query$e$vr[1:cluster_selec[1]]])
  
  for(j in 1:length(factor_raw)){
    tmp <- matrix("FALSE", ncol = length(factor_names[[j]]), nrow = cluster_selec[1],
                  dimnames = list(query$CC_desc$CC[query$e$vr[1:cluster_selec[1]]], factor_names[[j]]))
    for(k in 1:nrow(tmp)){
      for(l in 1:ncol(tmp)){
        if(str_detect(factor_raw[[j]][k], as.character(factor_names[[j]][l])) == TRUE){tmp[k,l] <- "TRUE"}
      } # factor_names
    } # CC names
    df <- cbind(df, tmp)
  }# j variable groups

  res_mfa <- MFA(base = df, group = c(1,1,16,35,14), type = c("n", "s" ,"n", "n", "n"), 
                 name.group = c("group","unknown_rate","Module","KO","Class"), graph = FALSE)
  
  # Layout
  pdf(paste0(bluecloud_dir, "/output/", output_dir, "/MFA.pdf"))
  layout(mat = matrix(c(1,1,2,3,4,5), ncol = 3, nrow = 2, byrow = TRUE), widths = c(2,2,1), heights = c(2,1))
  # Plot of individuals
  par(mar = c(4,4,2,1))
  plot(res_mfa$ind$coord[,1], res_mfa$ind$coord[,2], pch = 20, main = "Individuals",
       xlab = paste0("Dimension 1 (", round(res_mfa$eig[1,2],2),"%)"), 
       ylab = paste0("Dimension 2 (",round(res_mfa$eig[2,2],2) ,"%)"),
       col = brewer.pal(tree_cut, "Set2")[group])
  text(res_mfa$ind$coord[,1], res_mfa$ind$coord[,2], labels = names(group), 
       col = brewer.pal(tree_cut, "Set2")[group], pos = 3, offset = 0.2, cex = 0.5)
  grid()
  abline(h = 0, v = 0, lty = "longdash")
  
  # Legend
  par(mar = c(0,0,0,0))
  plot.new()
  legend(0, 0.7, legend = paste("Group", seq(1:5)), fill = brewer.pal(tree_cut, "Set2"), )
  
  # Dimension 1 : contribution
  par(mar = c(6,4,2,2))
  barplot(sort(res_mfa$quali.var$contrib[,1], decreasing = TRUE)[1:20], las = 2,
          cex.names = 0.5, main  = "Quali.var contribution to dim. 1", ylab = "(%)",
          col = rep(brewer.pal(4, "Set2"), res_mfa$call$group.mod[-2])[order(res_mfa$quali.var$contrib[,1], decreasing = TRUE)[1:20]])
  abline(h = mean(res_mfa$quali.var$contrib[,1]), lty = "longdash")
  
  # Dimension 2 : contribution
  barplot(sort(res_mfa$quali.var$contrib[,2], decreasing = TRUE)[1:20], las = 2,
          cex.names = 0.5, main  = "Quali.var contribution to dim. 2", ylab = "(%)",
          col = rep(brewer.pal(4, "Set2"), res_mfa$call$group.mod[-2])[order(res_mfa$quali.var$contrib[,2], decreasing = TRUE)[1:20]])
  abline(h = mean(res_mfa$quali.var$contrib[,2]), lty = "longdash")
  
  # Legend
  par(mar = c(0,0,0,0))
  plot.new()
  legend(0, 0.7, legend = res_mfa$call$name.group[-2], fill = brewer.pal(4, "Set2"))
  dev.off()
  
} # global run loop

toto <- order(res_mfa$quali.var$contrib[,1], decreasing = TRUE)[1:20]


# --- END