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
kegg_p0 = c("CCMlight", "Glu", "N", "C4ko","C4ko", "C4ko","Nglu")
kegg_m0 = list(paste0("K",c("01601","01602",
                            "01743","01674","18245","18246",
                            "00028","00029","01610")),
               paste0("K",c("00265","00264","00284")),
               paste0("K",c("00367","10534","00372","00360","00366","17877", #NR + NiR
                            "01948","00611","01755", # Urea in and out to citrate
                            "01915","00265","00264","00284")), # GS I to III)
               paste0("K",c("01595","00051","00028","00029","00814","14272","01006","14454","14455","00024","00025","00026","01610")),
               paste0("K",c("01595","00051","00028","00029","00814","14272","01006","14454","14455","00024","00025","00026","01610")),
               paste0("K",c("01595","00051","00028","00029","00814","14272","01006","14454","14455","00024","00025","00026","01610")),
               paste0("K",c("00367","10534","00372","00360","00366","17877", #NR + NiR
                            "01915","00265","00264","00284"))) # GS I to III

cluster_selec0 = list(c(50,1,1),
                      c(50,1,0),
                      c(50,1,0),
                      c(50,1,1),
                      c(50,3,1),
                      c(50,1,0.5),
                      c(50,1,0.5))

for(p in 1:length(kegg_p0)){
  kegg_p = kegg_p0[p]
  kegg_m = kegg_m0[[p]]
  cc_id = NULL
  cluster_selec = cluster_selec0[[p]]
  env_metric = c("mean","sd")
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
                      ENV_METRIC = c(env_metric),
                      relative = TRUE)
  
  file.copy(from = paste0(data_dir, "/X.feather"), to = paste0(bluecloud_dir, "/output/", output_dir))
  file.copy(from = paste0(data_dir, "/Y.feather"), to = paste0(bluecloud_dir, "/output/", output_dir))
  file.copy(from = paste0(data_dir, "/CC_desc.feather"), to = paste0(bluecloud_dir, "/output/", output_dir))
  
  # Number of genes, unknown and escoufier selection
  # if(is.null(cc_id)==TRUE){
  #   pdf(paste0(bluecloud_dir, "/output/", output_dir, "/CC_desc.pdf"))
  #   layout(mat = matrix(c(1,2), ncol = 2, nrow = 1, byrow = TRUE), widths = c(2,1))
  #   par(mar = c(2,6,3,1))
  #   barplot(query$CC_desc$n_genes[query$e$vr], horiz = TRUE, xlim = c(0,50),
  #           col = rep(c("#536f8c","gray80"), each = c(min(length(query$e$vr),cluster_selec[1]), length(query$e$vr)-min(length(query$e$vr),cluster_selec[1]))), 
  #           main = "Genes per connected components", names.arg = query$CC_desc$CC[query$e$vr], las = 2, axes = FALSE, cex.names = 0.6)
  #   barplot(query$CC_desc$n_genes[query$e$vr]*query$CC_desc$unknown_rate[query$e$vr]/100, add = TRUE, horiz = TRUE,
  #           col = rep(c("#28425c","gray50"), each = c(min(length(query$e$vr),cluster_selec[1]), length(query$e$vr)-min(length(query$e$vr),cluster_selec[1]))))
  #   par(new = TRUE)
  #   abline(v = c(2,4,6,8,seq(10,50,10)), lty = "dotted")
  #   abline(v = 0, h = 1.2*30, lwd = 2, lty = c("dashed","solid"), col = c("#536f8c","black"))
  #   par(mar = c(2,0,3,1))
  #   plot(x = query$e$RV, y = 1:length(query$e$RV), ylab = "", main = "Escoufier rel. imp.", type = 'l', lwd = 2, axes = FALSE)
  #   axis(side = 1, at = seq(0.4,1,0.1), labels = seq(0.4,1,0.1))
  #   abline(v = c(seq(0.4,1,0.1)), lty = "dotted")
  #   abline(v = 0.4, h = 30.2, lwd = 2, lty = c("dashed","solid"), col=c("#536f8c","black"))
  #   dev.off()
  # }

  # Assign points to the nearest escoufier selected point (correspondent analysis)
  res_ca <- CA(query$Y0, graph = FALSE, ncp = 50)
  plot(res_ca$col$coord[,1], res_ca$col$coord[,2], pch = 19,
       xlab = paste0("Dimension 1 (", round(res_ca$eig[1,2],2),"%)"), 
       ylab = paste0("Dimension 2 (",round(res_ca$eig[2,2],2) ,"%)"))
  points(res_ca$col$coord[1:cluster_selec[1],1], res_ca$col$coord[1:cluster_selec[1],2], pch = 19, col = "red")
  ndim_ca <- which(res_ca$eig[,3]>80)[1]
  dist_ca <- dist(res_ca$col$coord[,1:ndim_ca]) %>% as.matrix()
  dist_ca <- dist_ca[,1:cluster_selec[1]]

  nn_ca <- data.frame(CC = dimnames(dist_ca)[[1]],
                      pos_CC = 1:nrow(dist_ca),
                      nn_CC = apply(dist_ca, 1, function(x){x = names(which(x == min(x)))}),
                      pos_nn_CC = apply(dist_ca, 1, function(x){x = which(x == min(x))}))

  query[["nn_ca"]] <- nn_ca
  
  # 2. Run model & save ----------------------------------------------------------
  run <- model_run(bluecloud.wd = bluecloud_dir,
                   HYPERPARAMETERS = data.frame(LEARNING_RATE = c(1e-2, 1e-2, 1e-2, 1e-2),
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
                     var_importance = TRUE)
  
  # 4. Do projections & save -----------------------------------------------------
  proj <- model_proj(bluecloud.wd = bluecloud_dir,
                     data.wd = data_dir,
                     ENV_METRIC = c(env_metric, "dist","bathy"))
  
  save(query, eval, proj, file = paste0(bluecloud_dir, "/output/", output_dir, "/output.RData"))

} # global run loop    
  
  # 5. Do plot outputs & save ----------------------------------------------------
  # --- Create rescaled map values
  y_hat_m_rescaled <- apply(proj$y_hat_m, 2, function(x){x/max(x, na.rm = TRUE)})
  if(is.null(cc_id) == TRUE){colnames(y_hat_m_rescaled) <- query$CC_desc$CC[query$e$vr[1:min(length(query$e$vr),cluster_selec[1])]]
  } else {colnames(y_hat_m_rescaled) <- query$CC_desc$Genes[query$e$vr[1:min(length(query$e$vr),cluster_selec[1])]]}
  
  # --- Do spatial clustering ---
  tree <- hclust(as.dist(1-cor(y_hat_m_rescaled, use = "pairwise.complete.obs")), method = "ward.D2")
  tree_cut <- 8
  group <- cutree(tree, k = tree_cut)
  
  # --- Start PDF with clustering related plots
  pdf(paste0(bluecloud_dir, "/output/", output_dir, "/", output_dir, ".pdf"))
  par(mfrow = c(2,1),  mar = c(3,11,3,11))
  barplot(rev(tree$height), ylab = "Euclidian distance", xlab = "Number of clusters", main = "Connected Component (CC) cutoff",
          col = c(rep("black", tree_cut), rep("gray75", cluster_selec[1]-tree_cut)))
  abline(h = tree$height[length(tree$height)-tree_cut], col = "red")
  plot(tree, cex = 0.6, ylab= "Euclidian distance", main = "Connected Component (CC) dendrogram")
  abline(h = tree$height[length(tree$height)-(tree_cut-1)], col = "red")
  cor_proj(y_hat_m = y_hat_m_rescaled[,tree$order], targetNAME = tree$labels[tree$order])
  abline(h = match(1:tree_cut, group[tree$order])-0.5, v = match(1:tree_cut, group[tree$order])-0.5)
  legend_proj(col_matrix = proj$col_matrix)
  
  # --- Continue with map, labels and related nn_CC / nn_KO
  par(mfrow = c(4,2), mar=c(2,2,2,0))
  for(n in 1:min(length(query$e$vr),cluster_selec[1])){
    map_proj(proj = proj$proj, 
             col = proj$col, 
             targetID = tree$order[n], 
             targetNAME = tree$labels)
    map_scale <- max(proj$y_hat_m[,tree$order[n]], na.rm = TRUE)
    barplot(c(map_scale, rep(0, 9)), ylim = c(0,1), border = NA, col = "black",
            main = "Description :")
    abline(v = c(0.2,1.2))
    title(main = "Scale :", adj = 0)
    
    linked_CC <- query$CC_desc$CC[query$e$vr][which(query$nn_ca$nn_CC == tree$labels[tree$order[n]])][-1] %>% 
      paste(collapse = ",")
    linked_KO <- query$CC_desc$kegg_ko[query$e$vr][which(query$nn_ca$nn_CC == tree$labels[tree$order[n]])][-1] %>% 
      strsplit(",") %>% unlist() %>% unique() %>% paste(collapse = ",")
    
    text(x = rep(2,7), y = c(0.9,0.7,0.6,0.5,0.4,0.2,0.1), cex = 0.8, adj = c(0,1),
         labels = c(paste("H. Clustering group:", group[tree$order][n]),
                    paste("KOs:", query$CC_des$kegg_ko[query$e$vr[1:min(length(query$e$vr),cluster_selec[1])]][tree$order[n]]),
                    paste("Unknown rate:", query$CC_desc$unknown_rate[query$e$vr[1:min(length(query$e$vr),cluster_selec[1])]][tree$order[n]]),
                    paste("nn_CCs:", linked_CC),
                    paste("nn_KOs:", linked_KO),
                    paste("Class:", query$CC_desc$class[query$e$vr[1:min(length(query$e$vr),cluster_selec[1])]][tree$order[n]]),
                    paste("Genus:", query$CC_desc$genus[query$e$vr[1:min(length(query$e$vr),cluster_selec[1])]][tree$order[n]])))
  }
  dev.off()
  
  # 6. Multivariate analysis on map clusters for interpretation ----------------
  # Create table
  factor_raw <- list(query$CC_desc$kegg_module[query$e$vr[1:min(length(query$e$vr),cluster_selec[1])]],
                  query$CC_desc$kegg_ko[query$e$vr[1:min(length(query$e$vr),cluster_selec[1])]],
                  query$CC_desc$class[query$e$vr[1:min(length(query$e$vr),cluster_selec[1])]])
  factor_names <- lapply(factor_raw, function(x){x <- gsub(x, pattern = " ", replacement = "") %>% 
                                                       strsplit(split = ",") %>% 
                                                       unlist() %>% unique()
                                                 if(length(which(x == "-" | x == "NA")) > 0){
                                                    x <- x[-which(x == "-" | x == "NA")]
                                                 } else {x <-  x}
                                                })
  df <- data.frame(group = as.factor(group), unknown_rate = query$CC_desc$unknown_rate[query$e$vr[1:min(length(query$e$vr),cluster_selec[1])]])
  factor_names[[2]] <- kegg_m[paste0("ko:",kegg_m)%in%factor_names[[2]]]
  # Binary variables
  for(j in 1:length(factor_raw)){
    tmp <- matrix("FALSE", ncol = length(factor_names[[j]]), nrow = min(length(query$e$vr),cluster_selec[1]),
                  dimnames = list(query$CC_desc$CC[query$e$vr[1:min(length(query$e$vr),cluster_selec[1])]], factor_names[[j]]))
    for(k in 1:nrow(tmp)){
      for(l in 1:ncol(tmp)){
        if(str_detect(factor_raw[[j]][k], as.character(factor_names[[j]][l])) == TRUE){tmp[k,l] <- "TRUE"}
      } # factor_names
    } # CC names
    df <- cbind(df, tmp)
  }# j variable groups

  res_mfa <- MFA(base = df, group = c(1,1,length(factor_names[[1]]),length(factor_names[[2]]),length(factor_names[[3]])), type = c("n", "s" , "n", "n", "n"), 
                 name.group = c("group","unknown_rate","Module","KO","Class"), num.group.sup = c(1,2,3) ,graph = FALSE)

  # Layout
  pdf(paste0(bluecloud_dir, "/output/", output_dir, "/MFA.pdf"))
  layout(mat = matrix(c(1,1,2,3,4,5), ncol = 3, nrow = 2, byrow = TRUE), widths = c(2,2,1), heights = c(2,1))
  # Plot of individuals
  par(mar = c(4,4,2,1))
  dataEllipse(res_mfa$ind$coord[,1], res_mfa$ind$coord[,2], levels = 0.8, groups = as.factor(group),
              pch = 20, lwd = 1, main = "Individuals", plot.points = TRUE, robust = TRUE,
              xlab = paste0("Dimension 1 (", round(res_mfa$eig[1,2],2),"%)"), 
              ylab = paste0("Dimension 2 (",round(res_mfa$eig[2,2],2) ,"%)"),
              col = brewer.pal(tree_cut, "Set2"))
  points(res_mfa$ind$coord[,1], res_mfa$ind$coord[,2], pch = 20, 
       col = brewer.pal(tree_cut, "Set2")[group])
  text(res_mfa$ind$coord[,1], res_mfa$ind$coord[,2], labels = names(group), 
       col = brewer.pal(tree_cut, "Set2")[group], pos = 3, offset = 0.2, cex = 0.5)
  points(res_mfa$quali.var$coord[seq(2,nrow(res_mfa$quali.var$coord),2),1], res_mfa$quali.var$coord[seq(2,nrow(res_mfa$quali.var$coord),2),2], pch = 20, 
         col = "black")
  text(res_mfa$quali.var$coord[seq(2,nrow(res_mfa$quali.var$coord),2),1], res_mfa$quali.var$coord[seq(2,nrow(res_mfa$quali.var$coord),2),2], 
       labels = rownames(res_mfa$quali.var$coord)[seq(2,nrow(res_mfa$quali.var$coord),2)], 
       col = "black", pos = 3, offset = 0.2, cex = 0.5)
  grid()
  abline(h = 0, v = 0, lty = "longdash")

  # Legend
  par(mar = c(0,0,0,0))
  plot.new()
  legend(0, 0.7, legend = paste("Group", seq(1:tree_cut)), fill = brewer.pal(tree_cut, "Set2"), )
  
  # Dimension 1 : contribution
  par(mar = c(6,4,2,2))
  barplot(sort(res_mfa$quali.var$contrib[,1], decreasing = TRUE)[1:20], las = 2,
          cex.names = 0.5, main  = "Quali.var contribution to dim. 1", ylab = "(%)",
          col = rep(brewer.pal(3, "Set2")[1:2], res_mfa$call$group.mod[-c(1,2,3)])[order(res_mfa$quali.var$contrib[,1], decreasing = TRUE)[1:20]])
  abline(h = mean(res_mfa$quali.var$contrib[,1]), lty = "longdash")
  
  # Dimension 2 : contribution
  barplot(sort(res_mfa$quali.var$contrib[,2], decreasing = TRUE)[1:20], las = 2,
          cex.names = 0.5, main  = "Quali.var contribution to dim. 2", ylab = "(%)",
          col = rep(brewer.pal(3, "Set2")[1:2], res_mfa$call$group.mod[-c(1,2,3)])[order(res_mfa$quali.var$contrib[,2], decreasing = TRUE)[1:20]])
  abline(h = mean(res_mfa$quali.var$contrib[,2]), lty = "longdash")
  
  # Legend
  par(mar = c(0,0,0,0))
  plot.new()
  legend(0, 0.7, legend = res_mfa$call$name.group[-c(1,2,3)], fill = brewer.pal(3, "Set2")[1:2])
  dev.off()
  
  # 7. Specific correlations to test -------------------------------------------
  # Defining the different maps to plot
  plot_list <- list(RUBISCO = "01601|01602",
                    C4 = "00028|00029|01610",
                    PEPCK = "01610",
                    ME_NADP = "00029",
                    ME_NAD = "00028")

  plot_list <- list(NR = "00367|00372|00360|10534",
                    NiR  = "00366|17877",
                    Glu1 = "264",
                    Glu2 = "265",
                    Glu3 = "284")
  
  plot_list <- list(PPC = "1595",
                    GOT = "14454|14455",
                    PEPCK = "01610",
                    MDH = "00024|00025|00026",
                    MD = "00051",
                    ME_NADP = "00029",
                    ME_NAD = "00028",
                    GPT = "00814|14272",
                    PPDK = "1006")

  CC_desc_e <- query$CC_desc[query$e$vr,] %>% inner_join(query$nn_ca)
  r0 <- stack(paste0(data.wd,"/features"))[[1]]
  
  proj_data <- y_hat_m_rescaled
  proj_out <- NULL
  taxo_out <- NULL

  par(mfrow = c(5,2))
  for(j in 1:length(plot_list)){
    # Extract spatial data
    id <- CC_desc_e$pos_nn_CC[which(str_detect(CC_desc_e$kegg_ko, plot_list[[j]])==TRUE)]
    print(paste("---", names(plot_list[j]), "// Nb of CC :", length(id), "---"))
    
    tmp <- apply(proj_data[,id], 1, mean)
    tmp_proj <- setValues(r0, tmp)
    
    par(mar = c(2,1,2,2))
    plot(tmp_proj, main = paste(names(plot_list[j]), "// Nb of CC :", length(id)), col = colorRampPalette(col = rev(brewer.pal(10, "Spectral")))(100))
    
    # Taxonomic proportions
    df <- matrix(0, nrow = length(id), ncol = length(factor_names[[3]]), dimnames = list(CC_desc_e$CC[id], factor_names[[3]]))
    for(k in 1:dim(df)[[1]]){
      for(l in 1:dim(df)[[2]]){
        if(str_detect(factor_raw[[3]][id[k]], factor_names[[3]][l])==TRUE){df[k,l] <- df[k,l]+sum(proj_data[,id[k]], na.rm = TRUE)}
      }
    }
    df <- apply(df,2,sum)
    df <- df/sum(df)
    
    par(mar=c(5,2,1,0))
    barplot(df, las = 2, cex.names = 0.6, col = "gray20", ylim = c(0,1))
    abline(h = c(0.2,0.4,0.6,0.8), lty = "dotted")
    box()
    
    #  save for later
    proj_out <- cbind(proj_out, tmp)
    taxo_out <- cbind(taxo_out, df)
  }
  
  # naming saved outputs
  colnames(proj_out) <- colnames(taxo_out) <- names(plot_list)
  
  # Supplementary plots
  library(vegan)
  map_similarity <- as.dist(cor(proj_out, use = "pairwise.complete.obs"))
  taxo_similarity <- 1-vegdist(t(taxo_out), "bray")
  mantel(taxo_similarity, map_similarity, method = "pearson", permutations = 1e+5)
  print(map_similarity)

  
  


# --- END