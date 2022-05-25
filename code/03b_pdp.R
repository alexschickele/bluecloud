
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
X0 <- as.data.frame(read_feather(paste0(bluecloud.wd,"/data/X.feather")))
Y0 <- as.data.frame(read_feather(paste0(bluecloud.wd,"/data/Y.feather")))

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

# --- Building the Partial Dependence Plot (PDPs)
all_pdp <- list() # final storage structure

for(i in 1:dim(X_tr)[2]){
  all_pdp[[names(features)[i]]] <- list() # sub-structure by environmental variable
  mpdp <- NULL
  
  # --- Calculating PDPs by bootstrap
  for(b in 1:NBOOTSTRAP){
    # --- Reload model and train data from bootstrap predictions
    X_tr <- X0 # to avoid different train, hence predict data set for PDPs (i.e., different X axis)
    m <- proj$m[[b]]
    
    # --- Calculating PDPs
    mpdp0 <- mpartial(object = m, pred.var = i, grid.resolution = 20, train = X_tr, cores = min(c(MAX_CLUSTER, NBOOTSTRAP)))
    all_pdp[[i]][["grid"]] <- mpdp0[,1]
    mpdp0 <- mpdp0[,-1]
    mpdp0 <- apply(mpdp0, 2, function(x){x = x/sum(x, na.rm = TRUE)}) # sum columns = 1
    colnames(mpdp0) <- names(Y0)
    
    mpdp <- abind(mpdp, mpdp0, along = 3)
  } # b bootstrap loop
  cat(paste(Sys.time(), " --- mpdp computed ---", names(features)[i], "\n"))

  # --- Aggregating by functions
  for(j in 1:length(plot_list)){
    # Extract nearest neighbor data
    id <- CC_desc_e$pos_nn_CC[which(str_detect(CC_desc_e$kegg_ko, plot_list[[j]])==TRUE)]
    if(scaled == TRUE){scale_CC <- query$nn_ca$sum_CC[which(str_detect(CC_desc_e$kegg_ko, plot_list[[j]])==TRUE)]} else {scale_CC <- 1}
    
    tmp <- apply(mpdp[,id,], c(1,3), function(x){x = x*scale_CC}) %>% aperm(c(2,1,3)) # select and scale patterns corresponding to enzyme
    tmp <- apply(tmp, c(1,3), sum)                                                    # aggregate all patterns by sum
    
    all_pdp[[i]][["pdp"]] <- cbind(all_pdp[[i]][["pdp"]], apply(tmp, 1, mean))        # all enzymes together (MEAN)
    all_pdp[[i]][["cv"]] <- cbind(all_pdp[[i]][["cv"]], apply(tmp, 1, cv))            # all enzymes together (CV)
    
    colnames(all_pdp[[i]][["pdp"]])[j] <- colnames(all_pdp[[i]][["cv"]])[j] <- names(plot_list)[j]  # pretty names
  } # j aggregating loop
  cat(paste(Sys.time(), " --- mpdp aggregated ---", names(features)[i], "\n"))
} # i th feature
  
  
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
    
    # Computing 1-CV-based confidence interval
    tmp_cv <- tmp*all_pdp[[i]][["cv"]][,j]/100
    
    if(j == 1){
      plot(all_pdp[[i]][["grid"]], tmp, type = 'l', lwd = 1, ylim = c(0,1), col = pal[j],
           xlab = "", main = names(features)[i])
      polygon(x = c(all_pdp[[i]][["grid"]], rev(all_pdp[[i]][["grid"]])),
              y = c(tmp-tmp_cv, rev(tmp+tmp_cv)),
              col = alpha(pal[j], 0.3), border = NA)
      mtext(side = 4, at = tail(tmp, 1), text = names(plot_list)[j], col = pal[j], padj = 0.5, las = 1, cex = 0.6)
      grid(col = "gray20")
    } else {
      lines(all_pdp[[i]][["grid"]], tmp, lwd = 1, col = pal[j])
      polygon(x = c(all_pdp[[i]][["grid"]], rev(all_pdp[[i]][["grid"]])),
              y = c(tmp-tmp_cv, rev(tmp+tmp_cv)),
              col = alpha(pal[j], 0.3), border = NA)
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

