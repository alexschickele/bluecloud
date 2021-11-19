#' @concept re-trains MBTR on all data according to a bootstrap procedure and the
#' hyperparameters selected with script 03_ and 04_
#' 
#' @source MBTR python library
#' 
#' @param bluecloud.wd path to the bluecloud descriptor file
#' @param HYPERPARAMETERS dataframe of hyperparameter to test in the model
#' @param NBOOST number of maximum boosting round to perform
#' @param NBOOTSTRAP number of bootstrap rounds to perform
#' @param MAX_CLUSTER maximum CPU clustering for parallel computing
#' 
#' @return .mbtr object per bootstrap rounds
#' @return .pdf in /graphic containing mapped projections per target and
#' interaction matrix between targets
#' 
#' TO DO LIST : 
#' - do more pretty graphics
#' - think about a synthetic way to present the graphics

model_proj <- function(bluecloud.wd = "/home/jovyan/bluecloud",
                       data.wd = "/home/jovyan/dataspace/PlanktonGenomic_datasets/",
                       ENV_METRIC = c("mean","sd","med","mad","dist","bathy")){
  
  HYPERPARAMETERS <- read_feather(paste0(bluecloud.wd,"/data/HYPERPARAMETERS.feather"))
  
  # --- Load files
  features <- stack(paste0(data.wd,"/features"))
  features <- features[[grep(paste(ENV_METRIC, collapse = "|"), names(features))]]
  Y0 <- as.data.frame(read_feather(paste0(bluecloud.wd,"/data/Y.feather")))
  X0 <- as.data.frame(read_feather(paste0(bluecloud.wd,"/data/X.feather")))
  N <- nrow(X0)
  X <- as.data.frame(getValues(features))
  
  # ================ PART 1 : generating predictions with bootstrap ==============
  # --- Create bootstrap X_tr and Y_tr
  for(b in 1:NBOOTSTRAP){
    id <- sample(seq(1:nrow(X0)), replace = TRUE)
    X_tr <- X0[id,]
    Y_tr <- Y0[id,]
    
    write_feather(X_tr, paste0(bluecloud.wd,"/data/",b,"_X_tr.feather"))
    write_feather(Y_tr, paste0(bluecloud.wd,"/data/",b,"_Y_tr.feather"))
  } # bootstrap loop
  
  # --- Re-fitting models on all data
  b <- 1:NBOOTSTRAP
  m <- mcmapply(FUN=mbtr_fit, 
                path=paste0(bluecloud.wd, "/data/", b),
                hp_id = as.character(b),
                loss_type='mse',
                n_boosts = as.integer(HYPERPARAMETERS$n_boost),
                min_leaf= HYPERPARAMETERS$MEAN_LEAF,
                learning_rate=HYPERPARAMETERS$LEARNING_RATE,
                lambda_weights=HYPERPARAMETERS$LEARNING_RATE/100,
                lambda_leaves=HYPERPARAMETERS$LEARNING_RATE*0,
                n_q= as.integer(HYPERPARAMETERS$N_Q),
                early_stopping_rounds = as.integer(HYPERPARAMETERS$n_boost),
                SIMPLIFY = FALSE,
                USE.NAMES = FALSE,
                mc.cores = min(c(MAX_CLUSTER, NBOOTSTRAP)))
  
  # --- Reload best_hp models from python because of "previous session invalidity"
  m <- list()
  for(b in 1:NBOOTSTRAP){
    m0 <- py_load_object(paste0(bluecloud.wd,"/data/",b,"_",b,"_m"), pickle = "pickle")
    m0 <- m0[[1]]
    m <- append(m, list(m0))
  } #b loop
  
  # --- Predicting values
  cat(paste("---", Sys.time(), "// Projecting maps // This may take some time --- \n"))
  y_hat <- mclapply(m,
                function(a_boot) mbtr_predict(a_boot, X),
                mc.cores = min(c(MAX_CLUSTER, NBOOTSTRAP))) %>% 
    abind(along = 3)
  cat(paste("---", Sys.time(), "// Done --- \n"))
  
  # ==================  PART 2 : calculate projections ============================
  # --- Initializing raster and projections data
  r0 <- raster(res=res(features), ext=extent(features))
  
  y_hat_m <- apply(y_hat, c(1,2), mean)
  y_hat_sd <- apply(y_hat, c(1,2), sd)
  
  cutx <- seq(0,0.5,0.05)
  cuty <- seq(0,1,0.1)
  
  # --- Rasterizing bivariate projections
  col_matrix <- colmat(pal = brewer.pal(10, "RdYlBu"), saturation = 0,
                       xlab = "Standard deviation", ylab = "Relative Abundance")
  
  for(i in 1:ncol(y_hat_m)){
    r_m <- setValues(r0, y_hat_m[,i])
    r_sd <- setValues(r0, y_hat_sd[,i])
    tmp <- bivar_map(rasterx = r_sd, rastery = r_m, colormatrix = col_matrix, cutx = cutx, cuty = cuty)
    
    if(i == 1){
      proj <- tmp[[1]]
      col <- list(tmp[[2]])
    } else {
      proj <- stack(proj, tmp[[1]])
      col[[i]] <- tmp[[2]]
    }
  } # i target loop
  proj <- synchroniseNA(stack(features[[1]], proj))[[-1]]
  names(proj) <- paste(1:ncol(y_hat_m))
  
  return(list(proj = proj, col_matrix = col_matrix, col = col, cutx = cutx, cuty = cuty, y_hat_m = y_hat_m))
} # end function

# ==================  PART 3 : plotting projections ============================
# --- Plot the legend ---
legend_proj <- function(col_matrix, cutx, cuty){
  par(mar=c(5,5,1,5))
  colmat_plot(col_matrix, xlab = "Standard deviation", ylab = "Relative Abundance")
  axis(side = 1, at = seq(0,1,0.1), labels = cutx)
  axis(side = 2, at = seq(0,1,0.1), labels = cuty)
}

# --- Plot the correlation ---
cor_proj <- function(y_hat_m){
  par(mar = c(2,2,3,5))
  proj_cor <- cor(y_hat_m) %>% 
    raster(xmn = 0.5, ymn = 0.5, xmx = ncol(y_hat_m)+0.5, ymx = ncol(y_hat_m)+0.5)
  proj_cor <- flip(proj_cor, direction = "y") %>%
    plot(col = brewer.pal(10, "RdBu"), main = "pair-wise correlation", axes = FALSE)
  
  axis(side = 1, at = seq(0,ncol(y_hat_m),1), labels = seq(0,ncol(y_hat_m),1))
  axis(side = 2, at = seq(0,ncol(y_hat_m),1), labels = seq(0,ncol(y_hat_m),1))
  grid(ncol(y_hat_m), lwd = 2, col = "white", lty = "solid")
}

# --- Plot the maps ---
map_proj <- function(proj, col, targetID = seq(1:nlayers(proj))){
  par(mar=c(2,2,2,1))
  for(i in targetID){
    plot(proj[[i]], col = col[[i]], legend = FALSE, ylim=c(-90, 90), 
         main= paste("Spatial relative abundance of target cluster n°",i))
    abline(h=c(66, 23, 0, -23, -66), lty = c("dotted","dotted","dashed","dotted","dotted"), lwd = c(1,1,1,1,1), col = "black")
  }
}





