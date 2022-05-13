
# ==============================================================================
# MA VERSION QUI MARCHE PAS :'(
# ==============================================================================

bluecloud.wd = bluecloud_dir
data.wd = data_dir
ENV_METRIC = c(env_metric)

HYPERPARAMETERS <- read_feather(paste0(bluecloud.wd,"/data/HYPERPARAMETERS.feather"))

# --- Load files
features <- stack(paste0(data.wd,"/features"))
features <- features[[grep(paste(ENV_METRIC, collapse = "|"), names(features))]]
Y0 <- as.data.frame(read_feather(paste0(bluecloud.wd,"/data/Y.feather")))
X0 <- as.data.frame(read_feather(paste0(bluecloud.wd,"/data/X.feather")))
# N <- nrow(X0)
# X <- as.data.frame(getValues(features))

# --- Create X_tr and Y_tr
b <- 1 # i.e. kept from bootstrap script, set to 1
X_tr <- X0
Y_tr <- Y0

write_feather(X_tr, paste0(bluecloud.wd,"/data/",b,"_X_tr.feather"))
write_feather(Y_tr, paste0(bluecloud.wd,"/data/",b,"_Y_tr.feather"))

# --- Re-fitting models on all data
m <- mbtr_fit(path=paste0(bluecloud.wd, "/data/", b),
              hp_id = as.character(b),
              loss_type='mse',
              n_boosts = as.integer(HYPERPARAMETERS$n_boost),
              min_leaf= HYPERPARAMETERS$MEAN_LEAF,
              learning_rate=HYPERPARAMETERS$LEARNING_RATE,
              lambda_weights=HYPERPARAMETERS$LEARNING_RATE/100,
              lambda_leaves=HYPERPARAMETERS$LEARNING_RATE*0,
              n_q= as.integer(HYPERPARAMETERS$N_Q),
              early_stopping_rounds = as.integer(HYPERPARAMETERS$n_boost))

# --- Reload best_hp models from python because of "previous session invalidity"
m <- py_load_object(paste0(bluecloud.wd,"/data/",b,"_",b,"_m"), pickle = "pickle")
m <- m[[1]]

# --- Defining function for multivariate PDP adapted to MBTR
mpartial <- function(object, train, pred.var, grid.resolution=10
                     , quantiles=FALSE, cores=10, ...) {
  #' @param object a fitted model for which a predict() method exists
  #' @param train training set on which the model was fitted
  #' @param pred.var name (or index) of the variable of the training set for which the pdp is to be computed
  #' @param grid.resolution number of points along the values of `pred.var` at which the pdp is to be computed
  #' @param quantiles when TRUE, choose those points as quantiles of the original variable
  #' @param ... ignored
  #'
  #' @return A data.frame with the values of `pred.var` and the predicted value of all response variables
  # define the grid of pred.var values at which to compute the pdp
  x <- select(train, pred.var)
  if (quantiles) {
    grid <- quantile(x, probs=seq(0,1,length.out=grid.resolution))
  } else {
    grid <- seq(min(x, na.rm=T), max(x, na.rm=T), length.out=grid.resolution)
  }
  
  # compute the pdp for each variable
  yhat <- parallel::mclapply(grid, function(val) {
    X <- train
    X[,all_of(pred.var)] <- val # Not sure about the 'all_of'
    this_yhat <- mbtr_predict(model = object, X_pred = X, n_boosts = HYPERPARAMETERS$n_boost)
    this_yhat <- apply(this_yhat, 2, mean)
  }, mc.cores=cores)
  
  # combine the results for all grid values
  yhat <- data.frame(x=grid, do.call(rbind, yhat))
  # add nice names
  names(yhat) <- c(pred.var, paste0(names(yhat)[-1], "hat"))
  
  return(yhat)
}

# --- Defining the different maps to plot
plot_list <- list(PPC = "1595",
                  GOT = "14454|14455",
                  PEPCK = "01610",
                  MDH_NAD = "00024|00025|00026",
                  MDH_NADP = "00051",
                  MDC_NADP = "00029",
                  MDC_NAD = "00028",
                  GPT_GGAT = "00814|14272",
                  PPDK = "1006")

# --- Supplementary parameters parameters
CC_desc_e <- query$CC_desc[query$e$vr,] %>% inner_join(query$nn_ca)
r0 <- stack(paste0(data.wd,"/features"))[[1]]
scaled <- TRUE

# --- Build PDP
all_pdp <- list()

for(i in 1:dim(X_tr)[2]){
  all_pdp[[names(features)[i]]] <- list()
  # --- Calculating PDPs
  mpdp <- mpartial(object = m, pred.var = i, grid.resolution = 20, train = X_tr, cores = 1)
  all_pdp[[i]][["grid"]] <- mpdp[,1]
  mpdp <- mpdp[,-1]
  mpdp <- apply(mpdp, 2, function(x){x = x/sum(x, na.rm = TRUE)}) # sum columns = 1
  colnames(mpdp) <- names(Y_tr)
  cat(paste(Sys.time(), " --- mpdp computed ---", names(features)[i], "\n"))
  
  # --- Aggregating by functions
  for(j in 1:length(plot_list)){
    # Extract nearest neighbor data
    id <- CC_desc_e$pos_nn_CC[which(str_detect(CC_desc_e$kegg_ko, plot_list[[j]])==TRUE)]
    if(scaled == TRUE){scale_CC <- query$nn_ca$sum_CC[which(str_detect(CC_desc_e$kegg_ko, plot_list[[j]])==TRUE)]} else {scale_CC <- 1}
    
    tmp <- apply(mpdp[,id],1, function(x){x = x*scale_CC}) %>% t() # select and scale patterns corresponding to enzyme
    tmp <- apply(tmp, 1, sum)                                      # aggregate all patterns by sum
    all_pdp[[i]][["pdp"]] <- cbind(all_pdp[[i]][["pdp"]], tmp)     # all enzymes together
    colnames(all_pdp[[i]][["pdp"]])[j] <- names(plot_list)[j]      # pretty names
  } # j aggregating loop
  cat(paste(Sys.time(), " --- mpdp aggregated ---", names(features)[i], "\n"))
} # i th feature
  
  
# --- Plot
if(scaled == TRUE){pdf(paste0(bluecloud_dir,"/output/", output_dir, "/PDP_scaled.pdf"))
}else{pdf(paste0(bluecloud_dir,"/output/", output_dir, "/PDP.pdf"))}
par(mfrow = c(3,3), mar = c(4,2,4,1))

for(i in 1:dim(X_tr)[2]){
  for(j in 1:length(plot_list)){
    # Get maximum per enzyme across all pdp
    max_j <- lapply(all_pdp, function(x){max(x[["pdp"]][,j])}) %>% unlist()
    
    # Rescaled at max = 1 for not yet defined reasons
    plot(all_pdp[[i]][["grid"]], all_pdp[[i]][["pdp"]][,j]/max(max_j), type = 'l', lwd = 3, ylim = c(0,1),
         xlab = names(features)[i], main = names(plot_list)[j])
    grid(col = "gray20")
  }
} # i th feature
dev.off()

# --- Plot between 0 and 1
if(scaled == TRUE){pdf(paste0(bluecloud_dir,"/output/", output_dir, "/PDP01_scaled.pdf"))
}else{pdf(paste0(bluecloud_dir,"/output/", output_dir, "/PDP01.pdf"))}
par(mfrow = c(3,3), mar = c(4,2,4,1))

for(i in 1:dim(X_tr)[2]){
  for(j in 1:length(plot_list)){
    # Get maximum per enzyme across all pdp
    max_j <- lapply(all_pdp, function(x){max(x[["pdp"]][,j])}) %>% unlist()
    min_j <- lapply(all_pdp, function(x){min(x[["pdp"]][,j])}) %>% unlist()
    
    # Rescaled at max = 1 for not yet defined reasons
    tmp <- (all_pdp[[i]][["pdp"]][,j]-min(min_j))/(max(max_j)-min(min_j))
    plot(all_pdp[[i]][["grid"]], tmp, type = 'l', lwd = 3, ylim = c(0,1),
         xlab = names(features)[i], main = names(plot_list)[j])
    grid(col = "gray20")
  }
} # i th feature
dev.off()

# --- Plot between 0 and 1 all together to better compare variables
pal = brewer.pal(9, "Paired")
if(scaled == TRUE){pdf(paste0(bluecloud_dir,"/output/", output_dir, "/PDP01_synthetic_scaled.pdf"))
}else{pdf(paste0(bluecloud_dir,"/output/", output_dir, "/PDP01_synthetic.pdf"))}
par(mfrow = c(3,3), mar = c(2,2,4,6))

for(i in 1:dim(X_tr)[2]){
  for(j in 1:length(plot_list)){
    # Get maximum per enzyme across all pdp
    max_j <- lapply(all_pdp, function(x){max(x[["pdp"]][,j])}) %>% unlist()
    min_j <- lapply(all_pdp, function(x){min(x[["pdp"]][,j])}) %>% unlist()
    
    # Rescaled at max = 1 for not yet defined reasons
    tmp <- (all_pdp[[i]][["pdp"]][,j]-min(min_j))/(max(max_j)-min(min_j))
    
    if(j == 1){
      plot(all_pdp[[i]][["grid"]], tmp, type = 'l', lwd = 1, ylim = c(0,1), col = pal[j],
           xlab = "", main = names(features)[i])
      mtext(side = 4, at = tail(tmp, 1), text = names(plot_list)[j], col = pal[j], padj = 0.5, las = 1, cex = 0.6)
      grid(col = "gray20")
    } else {
      lines(all_pdp[[i]][["grid"]], tmp, lwd = 1, col = pal[j])
      mtext(side = 4, at = tail(tmp, 1), text = names(plot_list)[j], col = pal[j], padj = 0.5, las = 1, cex = 0.6)
    } # end if
    
  }
} # i th feature
dev.off()

# --- Plot features for comparison
pal <- colorRampPalette(col = rev(brewer.pal(10,"Spectral")))(100)

pdf(paste0(bluecloud_dir,"/output/", output_dir, "/features.pdf"))
par(mfrow=c(3,3))

for(j in 1:nlayers(features)){
  plot(features[[j]], col= pal, main = names(features)[j])
}
dev.off()

