#' Master script for running the R pipeline, open and save results as well as
#' performing global analysis on several model runs

# ==============================================================================
#                           MODELLING PIPELINE
# ==============================================================================

# ============================ LOADING THE CODE ================================
# --- Loading the code
bluecloud_dir <- "/home/aschickele/workspace/bluecloud"
data_dir <- paste0(bluecloud_dir, "/data")

setwd(bluecloud_dir)
source("./code/00a_config.R")
source("./code/01a_query_data.R")
source("./code/02a_model_param.R")
source("./code/02b_model_eval.R")
source("./code/03a_bootstrap_predict.R")

MAX_CLUSTER <- 20

# =========================== DEFINE PARAMETERS ================================
kegg_p0 = c("C4","C4_RUBISCO")
kegg_m0 = list(paste0("K",c("01595","00051","00028","00029","00814","14272","01006","14454","14455","00024","00025","00026","01610")),
               paste0("K",c("01601","01602","01595","00051","00028","00029","00814","14272","01006","14454","14455","00024","00025","00026","01610")))

cluster_selec0 = list(c(50,1,1),
                      c(50,1,1))

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
  
  # Assign points to the nearest escoufier selected point (correspondent analysis)
  res_ca <- CA(query$Y0, graph = FALSE, ncp = 50)
  plot(res_ca$col$coord[,1], res_ca$col$coord[,2], pch = 19,
       xlab = paste0("Dimension 1 (", round(res_ca$eig[1,2],2),"%)"), 
       ylab = paste0("Dimension 2 (",round(res_ca$eig[2,2],2) ,"%)"))
  points(res_ca$col$coord[1:cluster_selec[1],1], res_ca$col$coord[1:cluster_selec[1],2], pch = 19, col = "red")
  ndim_ca <- which(res_ca$eig[,3]>80)[1]
  dist_ca <- dist(res_ca$col$coord[,1:ndim_ca]) %>% as.matrix()
  dist_ca <- dist_ca[,1:cluster_selec[1]]
  
  # Distance plot
  tmp <- dist_ca[1:50,] %>% apply(2, min)
  hist(tmp[50:length(tmp)], breaks = 50, main = "distance to the nearest Escoufier neighboor", col = "gray50")
  abline(h = c(1,5), col = c("red","black"))
  grid(col  = "black")
  
  # Quantify dominance of CC
  tmp <- apply(query$Y0,2,sum)

  # Concatenate in dataframe
  nn_ca <- data.frame(CC = dimnames(dist_ca)[[1]],
                      pos_CC = 1:nrow(dist_ca),
                      nn_CC = apply(dist_ca, 1, function(x){x = names(which(x == min(x)))}),
                      pos_nn_CC = apply(dist_ca, 1, function(x){x = which(x == min(x))}),
                      sum_CC = tmp,
                      scale_CC = tmp/tmp[apply(dist_ca, 1, function(x){x = which(x == min(x))})])
  
  # Add it to query
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
                     ENV_METRIC = c(env_metric))
  
  save(query, eval, proj, file = paste0(bluecloud_dir, "/output/", output_dir, "/output.RData"))

} # global run loop    
  

# ==============================================================================
#                           STANDARD OUTPUTS
# ==============================================================================

# 1. Do plot outputs & save ----------------------------------------------------
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

# ==============================================================================
#                           SPECIFIC OUTPUTS
# ==============================================================================

# ============================ INITIALIZE PARAMETERS ===========================
# --- Create factor table
factor_raw <- list(query$CC_desc$kegg_module[query$e$vr],
                   query$CC_desc$kegg_ko[query$e$vr],
                   query$CC_desc$class[query$e$vr],
                   query$CC_desc$mag[query$e$vr],
                   query$CC_desc$phylum[query$e$vr])
factor_names <- lapply(factor_raw, function(x){x <- gsub(x, pattern = " ", replacement = "") %>% 
                                                     strsplit(split = ",") %>% 
                                                     unlist() %>% unique()
                                               if(length(which(x == "-" | x == "NA")) > 0){
                                                  x <- x[-which(x == "-" | x == "NA")]
                                               } else {x <-  x}
                                              })
factor_names[[2]] <- kegg_m[paste0("ko:",kegg_m)%in%factor_names[[2]]]

# --- Defining the different maps to plot
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

# --- Supplementary parameters parameters
CC_desc_e <- query$CC_desc[query$e$vr,] %>% inner_join(query$nn_ca)
r0 <- stack(paste0(data.wd,"/features"))[[1]]
scaled <- FALSE

proj_data <- apply(proj$y_hat, c(2,3), function(x){x = x/sum(x, na.rm = TRUE)}) 

# --- Colors
pal <- colorRampPalette(col = rev(brewer.pal(10,"Spectral")))(100)

# =================== BUILDING FUNCTIONAL DATA & PROJ ==========================
# 1. Initialize parameters -----------------------------------------------------
func_data <- NULL
func_r_pal <- list()

col_matrix <- colmat(pal = colorRampPalette(rev(brewer.pal(10,"Spectral")))(100), value = 0,
                     xlab = "Coef. Variation (%)", ylab = "Relative Abundance")

# 2. Building data, rasters and profiles ---------------------------------------
for(j in 1:length(plot_list)){
  # --- Extract nearest neighbor data
  id <- CC_desc_e$pos_nn_CC[which(str_detect(CC_desc_e$kegg_ko, plot_list[[j]])==TRUE)]
  if(scaled == TRUE){scale_CC <- query$nn_ca$sum_CC[which(str_detect(CC_desc_e$kegg_ko, plot_list[[j]])==TRUE)]} else {scale_CC <- 1}
  
  # --- Building functional data
  tmp <- apply(proj_data[,id,],c(1,3), function(x){x = x*scale_CC}) # re-scale by raw data if necessary
  tmp <- apply(tmp, c(2,3), sum) # matrix transposed for some reasons...
  
  # --- Re-scaling the data between 0 and 1 now
  tmp <- apply(tmp, 2, function(x) (x = x/max(x, na.rm = TRUE)))
  tmp[tmp<0] <- 1e-10 # Negative values (i.e. NA later) are model artefact. 
                      # Rescaled to infinite small positive to be considered as
                      # 0 when using raster::cut() in bivarmap

  # --- Save
  func_data <- abind(func_data, tmp, along = 3)
  
  # --- Building functional raster and corresponding 3D color palette
  r_m <- setValues(r0, apply(tmp, 1, function(x) (x = mean(x, na.rm = TRUE))))
  r_cv <- setValues(r0, apply(tmp, 1, function(x) (x = cv(x, na.rm = TRUE))))
  r_cv[r_cv > 100] <- 100
  
  tmp <- bivar_map(rasterx = r_cv, rastery = r_m, colormatrix = col_matrix, cutx = seq(0,100,1), cuty = seq(0,1,0.01))
  
  if(j==1){func_r <- tmp[[1]]} 
  else {func_r <- stack(func_r, tmp[[1]])}
  func_r_pal[[j]] <- tmp[[2]]
}
# --- Rename data matrix
dimnames(func_data)[[3]] <- names(plot_list)

# --- Synchronise NA and rename raster
features <- stack(paste0(data.wd,"/features"))
func_r <- synchroniseNA(stack(features[[1]], func_r))[[-1]]
names(func_r) <- names(plot_list)

# ============== BUILDING TAXONOMIC DATA =======================================
# 1. Initialize parameters -----------------------------------------------------
taxo_data <- NULL

# 2. Building data, rasters and profiles ---------------------------------------
for(j in 1:length(factor_names[[3]])){
  # Extract taxonomic data
  id <- CC_desc_e$pos_nn_CC[which(str_detect(CC_desc_e$class, factor_names[[3]][[j]])==TRUE)]
  if(scaled == TRUE){scale_CC <- query$nn_ca$sum_CC[which(str_detect(CC_desc_e$class, factor_names[[3]][[j]])==TRUE)]} else {scale_CC <- 1}
  
  # Building taxonomic data
  if(length(id) > 1){
    tmp <- apply(proj_data[,id,], c(1,3), function(x){x = x*scale_CC})
    tmp <- apply(tmp, c(2,3), sum) # matrix transposed for some reasons...
  } else {
    tmp <- proj_data[,id,]*scale_CC
  } # if length ID > 1

  taxo_data <- abind(taxo_data, tmp, along = 3)
}
dimnames(taxo_data)[[3]] <- factor_names[[3]]

# ============== BUILDING MAG DATA =============================================
# 1. Initialize parameters -----------------------------------------------------
# Not including boostrap based uncertainty (no need)
mag_data <- NULL
proj_data0 <- apply(proj_data, c(1,2), mean)

# 2. Building data, rasters and profiles ---------------------------------------
for(j in 1:length(factor_names[[4]])){
  # Extract mag data
  id <- CC_desc_e$pos_nn_CC[which(str_detect(CC_desc_e$mag, factor_names[[4]][[j]])==TRUE)]
  if(scaled == TRUE){scale_CC <- query$nn_ca$sum_CC[which(str_detect(CC_desc_e$mag, factor_names[[4]][[j]])==TRUE)]} else {scale_CC <- 1}
  
  # Building mag data
  if(length(id)> 1){
    tmp <- apply(proj_data0[,id], c(1), function(x){x = x*scale_CC})
    tmp <- apply(tmp, c(2), sum) # matrix transposed for some reasons...
  } else {
    tmp <- proj_data0[,id]*scale_CC
  } # if length ID > 1
  
  mag_data <- abind(mag_data, tmp, along = 2)
}
dimnames(mag_data)[[2]] <- factor_names[[4]]

# =========================== GRAPHICAL OUTPUTS ================================
# 1. Functional maps -----------------------------------------------------------
r <- func_r

if(scaled == TRUE){pdf(paste0(bluecloud_dir,"/output/", output_dir, "/Functional_map_scaled.pdf"))
}else{pdf(paste0(bluecloud_dir,"/output/", output_dir, "/Functional_map.pdf"))}

par(mfrow = c(4,3), mar = c(3,2,3,2))
for(j in 1:length(plot_list)){
  map_proj(proj = r, 
           col = func_r_pal, 
           targetID = j, 
           targetNAME = names(plot_list))
}
legend_proj(col_matrix = col_matrix)
dev.off()

# 2. Correlation plots ---------------------------------------------------------
# --- Load additional libraries
library(vegan)
library(corrplot)

if(scaled == TRUE){pdf(paste0(bluecloud_dir,"/output/", output_dir, "/Correlations_scaled.pdf"))
}else{pdf(paste0(bluecloud_dir,"/output/", output_dir, "/Correlations.pdf"))}

par(mfrow = c(2,2), mar = c(2,5,5,3))
# --- Functional correlation
func_similarity <- as.dist(cor(apply(func_data, c(1,3), mean), use = "pairwise.complete.obs"))
corrplot(as.matrix(func_similarity), type = 'lower', diag = FALSE,
          tl.cex = 1, tl.col =  "black", order = 'FPC')

# --- Functional vs Taxonomic correlations
cor_mat <- cor(apply(func_data, c(1,3), mean), apply(taxo_data, c(1,3), mean), use = 'pairwise.complete.obs')
func_hclust <- hclust(dist(cor_mat), method = 'ward.D2')
taxo_hclust <- hclust(dist(t(cor_mat)), method = 'ward.D2')
cor_mat <- cor_mat[func_hclust$order, taxo_hclust$order]

corrplot(as.matrix(cor_mat), diag = FALSE,
         title = "Functional vs taxo composition similarity", tl.cex = 0.8, tl.col =  "black")

par(mar = c(10,5,5,2))
plot(as.dendrogram(func_hclust), main = "Functional dendrogram")
plot(as.dendrogram(taxo_hclust), main = "Taxonomic dendrogram")

dev.off()

# --- END
  
  