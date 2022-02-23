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

model_proj <- function(bluecloud.wd = bluecloud_dir,
                       data.wd = data_dir,
                       ENV_METRIC = c("mean","sd","dist","bathy")){
  
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
                function(a_boot) mbtr_predict(model = a_boot, X_pred = X, n_boosts = HYPERPARAMETERS$n_boost),
                mc.cores = min(c(MAX_CLUSTER, NBOOTSTRAP))) %>% 
    abind(along = 3)
  cat(paste("---", Sys.time(), "// Done --- \n"))
  
  # --- Back Hellinger like transformation
  # y_hat <- sqrt(y_hat, na.rm = TRUE)
  # cat(paste("--- Rescale post Hellinger done --- \n"))
  
  # ==================  PART 2 : calculate projections ============================
  # --- Initializing raster and projections data
  r0 <- raster(res=res(features), ext=extent(features))
  
  y_hat_m <- apply(y_hat, c(1,2), mean)
  y_hat_m[y_hat_m<0] <- 1e-10 #Negative values (i.e. NA later) are model artefact. 
                              #Rescaled to infinite small positive to be considered as
                              # 0 when using raster::cut() in bivarmap
  y_hat_cv <- apply(y_hat, c(1,2), cv)
  y_hat_cv[y_hat_cv > 100] <- 100  #CV over 100% is rescaled to 100% for interpretation purpose
  
  cutx <- seq(0,100,1)
  cuty <- seq(0,1,0.01)
  
  # --- Rasterizing bivariate projections
  custom_pal <- rev(brewer.pal(10,"Spectral"))
  col_matrix <- colmat(pal = colorRampPalette(custom_pal)(100), value = 0,
                       xlab = "Coef. Variation (%)", ylab = "Relative Abundance")
  y_hat_m_rescaled <- apply(y_hat_m, 2, function(x){x/max(x, na.rm = TRUE)})
  
  for(i in 1:ncol(y_hat_m)){
    r_m <- setValues(r0, y_hat_m_rescaled[,i])
    r_cv <- setValues(r0, y_hat_cv[,i])
    tmp <- bivar_map(rasterx = r_cv, rastery = r_m, colormatrix = col_matrix, cutx = cutx, cuty = cuty)
    
    if(i == 1){
      # proj <- tmp
      proj <- tmp[[1]]
      col <- list(tmp[[2]])
    } else {
      # proj <- stack(proj,  tmp)
      proj <- stack(proj, tmp[[1]])
      col[[i]] <- tmp[[2]]
    }
  } # i target loop
  
  # --- SynchroniseNA with coastline etc...
  proj <- synchroniseNA(stack(features[[1]], proj))[[-1]]
  names(proj) <- paste(1:ncol(y_hat_m))
  
  for(i in 1:ncol(y_hat_m)){
    y_hat_m[which(is.na(getValues(proj[[i]]))),i] <- NA
  }
  
  return(list(proj = proj, col_matrix = col_matrix, col = col, y_hat_m = y_hat_m))
} # end function

# ==================  PART 3 : plotting projections ============================
# --- Plot the legend ---
legend_proj <- function(col_matrix){
  colmat_plot(col_matrix, xlab = "Coef. Variation", ylab = "Relative Abundance")
  axis(side = 1, at = seq(0,1,0.2), labels = seq(0,100,20))
  axis(side = 2, at = seq(0,1,0.1), labels = seq(0,1,0.1))
}

# --- Plot the correlation ---
cor_proj <- function(y_hat_m, targetNAME = seq(1:ncol(y_hat_m))){
  proj_cor <- cor(y_hat_m, use = "pairwise.complete.obs") %>% 
    raster(xmn = 0.5, ymn = 0.5, xmx = ncol(y_hat_m)+0.5, ymx = ncol(y_hat_m)+0.5)
  proj_cor <- flip(proj_cor, direction = "y") %>%
    plot(col = colorRampPalette(col = rev(brewer.pal(9, "RdBu")))(100), main = "pair-wise correlation", axes = FALSE)
  
  axis(side = 1, at = seq(1,ncol(y_hat_m),1), labels = targetNAME, las = 2, cex.axis = 0.4)
  axis(side = 2, at = seq(1,ncol(y_hat_m),1), labels = targetNAME, las = 2, cex.axis = 0.4)
}

# --- Plot the maps ---
map_proj <- function(proj, col, targetID = seq(1:nlayers(proj)), targetNAME = seq(1:nlayers(proj))){
  # par(mar=c(2,2,2,1))
  for(i in targetID){
    plot(proj[[i]], col = col[[i]], legend = FALSE, ylim=c(-90, 90), 
         main= paste("Spatial relative abundance:",targetNAME[i]))
    abline(h=c(66, 23, 0, -23, -66), lty = c("dotted","dotted","dashed","dotted","dotted"), lwd = c(1,1,1,1,1), col = "black")
  }
}






