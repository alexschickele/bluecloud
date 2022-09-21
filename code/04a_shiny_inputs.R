
# ==============================================================================
#             PREPARATION OF THE RDATA TO BE USED AS INPUT OF SHINY APP
# ==============================================================================

# This script takes different sections of the master script and saves the
# outputs in a specific Rdata used for shiny only

# ============================ LOADING THE CODE ================================
# --- Loading the code
bluecloud_dir <- "/home/aschickele/workspace/bluecloud"
data_dir <- paste0(bluecloud_dir, "/data")

setwd(bluecloud_dir)
source("./code/00a_config.R")

# =========================== DEFINE PARAMETERS ================================
kegg_p0 = c("C4_RUBISCO_clean")
kegg_m0 = list(paste0("K",c("01601","01602","01595","00051","00028","00029","00814","14272","01006","14454","14455","00024","00025","00026","01610")))
cluster_selec0 = list(c(50,1,1))

kegg_p = kegg_p0[1]
kegg_m = kegg_m0[[1]]
exclude = "New_MAST-4|Oomycota"
cluster_selec = cluster_selec0[[1]]
env_metric = c("mean","sd")
relative = TRUE

output_dir <- paste0(kegg_p, "_", exclude, "_", 
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
  message(">>> This model has never been run ... try again <<<")
}

# ==============================================================================
#       Enzyme level raster, color matrix and corresponding CC IDs
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
scaled <- TRUE

proj_data <- apply(proj$y_hat, c(2,3), function(x){x = x/sum(x, na.rm = TRUE)}) 

# --- Colors
pal <- colorRampPalette(col = rev(brewer.pal(10,"Spectral")))(100)

# ======================== UNSCALED OUTPUT DATA ================================
# 1. Initialize parameters -----------------------------------------------------
func_r_pal_unscaled <- list()
col_matrix <- colmat(pal = colorRampPalette(rev(brewer.pal(10,"Spectral")))(100), value = 0,
                     xlab = "Coef. Variation (%)", ylab = "Relative Abundance")

# 2. Building data, rasters and profiles ---------------------------------------
for(j in 1:length(plot_list)){
  # --- Extract nearest neighbor data
  id <- CC_desc_e$pos_nn_CC[which(str_detect(CC_desc_e$kegg_ko, plot_list[[j]])==TRUE)]
  scale_CC <- 1
  
  # --- Building functional data
  tmp <- apply(proj_data[,id,],c(1,3), function(x){x = x*scale_CC}) # re-scale by raw data if necessary
  tmp <- apply(tmp, c(2,3), sum) # matrix transposed for some reasons...
  
  # --- Re-scaling the data between 0 and 1 now
  tmp <- apply(tmp, 2, function(x) (x = x/max(x, na.rm = TRUE)))
  tmp[tmp<0] <- 1e-10 # Negative values (i.e. NA later) are model artefact. 
                      # Rescaled to infinite small positive to be considered as
                      # 0 when using raster::cut() in bivarmap
  
  # --- Building functional raster and corresponding 3D color palette
  r_m <- setValues(r0, apply(tmp, 1, function(x) (x = mean(x, na.rm = TRUE))))
  r_cv <- setValues(r0, apply(tmp, 1, function(x) (x = cv(x, na.rm = TRUE))))
  r_cv[r_cv > 100] <- 100
  
  tmp <- bivar_map(rasterx = r_cv, rastery = r_m, colormatrix = col_matrix, cutx = seq(0,100,1), cuty = seq(0,1,0.01))
  
  if(j==1){func_r_unscaled <- tmp[[1]]} 
  else {func_r_unscaled <- stack(func_r_unscaled, tmp[[1]])}
  func_r_pal_unscaled[[j]] <- tmp[[2]]
}

# --- Synchronise NA and rename raster
features <- stack(paste0(data.wd,"/features"))
func_r_unscaled <- synchroniseNA(stack(features[[1]], func_r_unscaled))[[-1]]
names(func_r_unscaled) <- names(plot_list)

# ========================= SCALED OUTPUT DATA ================================
# 1. Initialize parameters -----------------------------------------------------
func_r_pal_scaled <- list()
col_matrix <- colmat(pal = colorRampPalette(rev(brewer.pal(10,"Spectral")))(100), value = 0,
                     xlab = "Coef. Variation (%)", ylab = "Relative Abundance")

# 2. Building data, rasters and profiles ---------------------------------------
for(j in 1:length(plot_list)){
  # --- Extract nearest neighbor data
  id <- CC_desc_e$pos_nn_CC[which(str_detect(CC_desc_e$kegg_ko, plot_list[[j]])==TRUE)]
  scale_CC <- query$nn_ca$sum_CC[which(str_detect(CC_desc_e$kegg_ko, plot_list[[j]])==TRUE)]
  
  # --- Building functional data
  tmp <- apply(proj_data[,id,],c(1,3), function(x){x = x*scale_CC}) # re-scale by raw data if necessary
  tmp <- apply(tmp, c(2,3), sum) # matrix transposed for some reasons...
  
  # --- Re-scaling the data between 0 and 1 now
  tmp <- apply(tmp, 2, function(x) (x = x/max(x, na.rm = TRUE)))
  tmp[tmp<0] <- 1e-10 # Negative values (i.e. NA later) are model artefact. 
  # Rescaled to infinite small positive to be considered as
  # 0 when using raster::cut() in bivarmap
  
  # --- Building functional raster and corresponding 3D color palette
  r_m <- setValues(r0, apply(tmp, 1, function(x) (x = mean(x, na.rm = TRUE))))
  r_cv <- setValues(r0, apply(tmp, 1, function(x) (x = cv(x, na.rm = TRUE))))
  r_cv[r_cv > 100] <- 100
  
  tmp <- bivar_map(rasterx = r_cv, rastery = r_m, colormatrix = col_matrix, cutx = seq(0,100,1), cuty = seq(0,1,0.01))
  
  if(j==1){func_r_scaled <- tmp[[1]]} 
  else {func_r_scaled <- stack(func_r_scaled, tmp[[1]])}
  func_r_pal_scaled[[j]] <- tmp[[2]]
}

# --- Synchronise NA and rename raster
features <- stack(paste0(data.wd,"/features"))
func_r_scaled <- synchroniseNA(stack(features[[1]], func_r_scaled))[[-1]]
names(func_r_scaled) <- names(plot_list)


# =============== DIFFERENCE SCALED - UNSCALED OUTPUT DATA =====================
# 1. Initialize parameters -----------------------------------------------------
func_r_pal_diff <- list()
col_matrix <- colmat(pal = colorRampPalette(rev(brewer.pal(10,"Spectral")))(100), value = 0,
                     xlab = "Coef. Variation (%)", ylab = "Relative Abundance")

# 2. Building data, rasters and profiles ---------------------------------------
for(j in 1:length(plot_list)){
  # --- Extract nearest neighbor data
  id <- CC_desc_e$pos_nn_CC[which(str_detect(CC_desc_e$kegg_ko, plot_list[[j]])==TRUE)]
  scale_CC <- query$nn_ca$sum_CC[which(str_detect(CC_desc_e$kegg_ko, plot_list[[j]])==TRUE)]
  
  # --- Building functional data
  tmp_scaled <- apply(proj_data[,id,],c(1,3), function(x){x = x*scale_CC}) %>% # re-scale by raw data if necessary
    apply(c(2,3), sum) %>% 
    apply(2, function(x) (x = x/max(abs(x), na.rm = TRUE)))
  
  tmp_unscaled <- apply(proj_data[,id,],c(1,3), function(x){x = x*1}) %>% # re-scale by raw data if necessary
    apply(c(2,3), sum) %>% 
    apply(2, function(x) (x = x/max(abs(x), na.rm = TRUE)))
  
  tmp <- tmp_scaled-tmp_unscaled
  
  # --- Building functional raster and corresponding 3D color palette
  r_m <- setValues(r0, apply(tmp, 1, function(x) (x = mean(x, na.rm = TRUE))))
  r_cv <- setValues(r0, apply(tmp, 1, function(x) (x = cv(x, na.rm = TRUE))))
  r_cv[r_cv > 100] <- 100
  
  tmp <- bivar_map(rasterx = r_cv, rastery = r_m, colormatrix = col_matrix, cutx = seq(0,100,1), cuty = seq(-1,1,0.02))
  
  if(j==1){func_r_diff <- tmp[[1]]} 
  else {func_r_diff <- stack(func_r_diff, tmp[[1]])}
  func_r_pal_diff[[j]] <- tmp[[2]]
}

# --- Synchronise NA and rename raster
features <- stack(paste0(data.wd,"/features"))
func_r_diff <- synchroniseNA(stack(features[[1]], func_r_diff))[[-1]]
names(func_r_diff) <- names(plot_list)

# ==============================================================================
#       CC level raster, color matrix and corresponding CC IDs
# ==============================================================================



# ==============================================================================
#                           Save all in a RData
# ==============================================================================

proj0 <- list(proj = proj$proj, col_matrix = proj$col_matrix, col = proj$col)
proj <- proj0

save(query, eval, proj, CC_desc_e,
     func_r_scaled, func_r_pal_scaled,
     func_r_unscaled, func_r_pal_unscaled,
     func_r_diff, func_r_pal_diff,
     file = paste0(data_dir, "/shiny_data.RData"))



