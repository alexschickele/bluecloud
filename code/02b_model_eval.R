#' @concept evaluate each validation fold with R2, MSE and RMSE and computes variable
#' importance
#' 
#' @source MBTR python library
#' 
#' @param bluecloud.wd path to the bluecloud folder containing the master script
#' @param by_target TRUE or FALSE, performs the above mentionned metrics by target as well
#' @param var_importance TRUE or FALSE, performs a variable importance output
#' 
#' @return an eval object with all above mentionned metrics

model_eval <- function(bluecloud.wd = bluecloud_dir,
                       by_target = TRUE,
                       var_importance = FALSE){
  
  # --- 1. Loading data
  # Here we load the MBTR models producing the minimum average loss for further performance evaluation
  setwd(paste0(bluecloud.wd,"/data/"))
  m <- py_load_object("m", pickle = "pickle")
  HYPERPARAMETERS <- read_feather(paste0(bluecloud.wd,"/data/HYPERPARAMETERS.feather"))
  # We also load the station ID, target and feature tables
  ID <- read_feather(paste0(bluecloud.wd,"/data/Station_FOLD.feather"))
  Y0 <- as.data.frame(read_feather(paste0(bluecloud.wd,"/data/Y.feather")))
  X0 <- as.data.frame(read_feather(paste0(bluecloud.wd,"/data/X.feather")))
  
  # --- 2. Initializing outputs
  y_hat <- r2_tar <- r2cor_tar <- rmse_tar <- NULL
  
  # --- 3. Predicting on the validation data
  # The models corresponding to the selected hyperparameter set and each n-cross-validation runs
  # were evaluated by predicting on the corresponding n-test sets, not seen before by the models.
  for(cv in 1:N_FOLD){
    # 3.1. Loading test data and models
    X_val <- as.data.frame(read_feather(paste0(bluecloud.wd,"/data/", cv, "_X_val.feather")))
    m0 <- m[[cv]][[1]]
    
    # 3.2. Do predictions
    y_hat <- rbind(y_hat, mbtr_predict(model = m0, X_pred = X_val, n_boosts = HYPERPARAMETERS$n_boost))
  } # k-fold cv  loop
  
  # 3.3. Concatenate the y_hat corresponding to each prediction on test folds and 
  # order them back with the pre-fold station order to match Y again before evaluation
  y_hat <- y_hat[order(ID$Station),]
  
  # --- 4. Testing the correlation structure integrity
  # During the modelling process, between the different targets
  dist_y_hat <- cor(y_hat) %>% as.dist()
  dist_Y0 <- cor(as.matrix(Y0)) %>% as.dist()
  message("--- Mantel test between Y0 and y_hat to test correlation structure")
  mantel <- mantel.rtest(dist_y_hat, dist_Y0)
  message(paste("--- Mantel correlation:", round(mantel$obs, 3)))
  message(paste("--- p-value:", mantel$pvalue))
  
  # --- 5. Evaluate the model performance
  # 5.1. Global model evaluation metrics
  r2 <- calc_rsquared(as.matrix(Y0), y_hat)
  rmse <- sqrt(mean(as.matrix((Y0-y_hat)^2), na.rm=TRUE))
  
  print(paste("--- model multidimensional R-squarred is :", round(r2, 2), "---"))
  print(paste("--- model multidimensional RMSE is :", round(rmse, 2), "---"))
  
  # 5.2. Target by target evaluation metrics
  if(by_target==TRUE){
    for(t in 1:ncol(Y0)){
      # numeric R2
      r2_tar <- c(r2_tar, calc_rsquared(as.matrix(Y0[,t]), as.matrix(y_hat[,t])))
      # Correlation R2
      r2cor_tar <- c(r2cor_tar, cor(as.matrix(Y0[,t]), as.matrix(y_hat[,t])))
      rmse_tar <- c(rmse_tar, sqrt(mean(as.matrix((Y0[,t]-y_hat[,t])^2), na.rm=TRUE)))
    }
    print("--- R2 by target :")
    print(round(r2_tar,2))
    print(round(r2cor_tar,2))
    print("--- RMSE by target :")
    print(round(rmse_tar,2))
  }
  
  # --- 6. Plot predictions vs truth for each target
  par(mar = c(4,4,5,0), mfrow = c(2,1))
  for(t in 1:ncol(Y0)){
    plot(y = Y0[,t], x = as.numeric(sort(ID$Station)), 
         type = 'p', xlim = c(0,153), ylab = "Abundance", xlab = "Station number", 
         main = paste("Target n°", t, "(", names(Y0[t]), ")"), cex.main = 1,
         pch = 16, axes = FALSE)
    segments(y0 = Y0[,t], x0 = as.numeric(sort(ID$Station)),
             y1 = y_hat[,t], x1 = as.numeric(sort(ID$Station)),
             col = "black")
    points(y = y_hat[,t], x = as.numeric(sort(ID$Station)), pch = 16, col = "gray50")
    axis(side = 1, at = seq(1:153), cex.axis = 0.5, las = 2)
    axis(side = 2)
    axis(side = 3, at = c(15, 48, 78, 115, 148), labels = c("Med.", "S. Ind.", "S. Alt.", "Pac.", "N. Alt"), tick = FALSE)
    abline(v = c(30.5, 66.5, 89.5, 140.5), lty = "dotted")
    box()
  }
  
  # --- 7. Calculating variable importance
  # Performed according to the loss evolution corresponding to each tree split, 
  # related to the corresponding environmental variable used for the split
  # 7.1. Initialize object
  var_count <- matrix(0, ncol = ncol(X0), nrow=N_FOLD)
  colnames(var_count) <- colnames(X0)
  
  # 7.2. Iteratively extract the loss corresponding to each tree split x environmental variable
  if(var_importance == TRUE){
    for(cv in 1:N_FOLD){
      m0 <- m[[cv]][[1]]
      n_tree <- length(m0$trees)
      var_count0 <- mcmapply(FUN = function(t, cv, nvar){
                                    tmp <- matrix(0, ncol = nvar, nrow = 1)
                                    n_node <-length(m0$trees[[t]]$g$nodes$`_nodes`)
                                    for(n in 1:n_node){
                                      var_nb <- m0$trees[[t]]$g$nodes$`_nodes`[[n]]$variable+1 #py index start at 0, R at 1
                                      if(!is.null(var_nb)){
                                        dloss <- m0$trees[[t]]$g$nodes$`_nodes`[[n]]$loss*(-1)
                                        tmp[var_nb] <- tmp[var_nb]+dloss
                                      }
                                      return(tmp)}
                                    }, # end if ; for ; function
                             t = 1:n_tree, cv = cv, nvar = ncol(X0),
                             SIMPLIFY = FALSE,
                             USE.NAMES = FALSE,
                             mc.cores = 20)
      var_count0 <- Reduce(`+`, var_count0)
      var_count[cv,] <- var_count0
    } # fold loop
    
    # 7.3. Plotting the variable importance
    var_imp_order <- t(apply(var_count, 1, function(x) (x*100)/sum(x, na.rm = TRUE))) %>% 
      apply(2, median) %>% 
      order(decreasing = TRUE)
    
    var_imp <- t(apply(var_count, 1, function(x) (x*100)/sum(x, na.rm = TRUE)))[,var_imp_order]
    
    boxplot(var_imp, axes = FALSE, ylab="variable importance (%)", xlab = "", col = "gray50")
    abline(h=(seq(0,100,10)), lty="dotted", col="black")
    axis(side = 2, at = seq(0,100,10), labels = seq(0,100,10))
    axis(side = 1, at = 1:ncol(X0), labels = colnames(X0)[var_imp_order], las = 2, 
         cex.axis = if(ncol(X0) > 20) {0.7} else {1})
    box()

  } # var imp
  return(list(r2 = r2, rmse = rmse, r2_tar = r2_tar, r_cor_tar = r2cor_tar, rmse_tar = rmse_tar, var_count = var_count, y_hat = y_hat, mantel = mantel))

} # end function
