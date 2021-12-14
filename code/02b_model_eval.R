#' @concept evaluate each validation fold with R2, MSE and RMSE and computes variable
#' importance
#' 
#' @source MBTR python library
#' 
#' @param bluecloud.wd path to the bluecloud descriptor file
#' @param HYPERPARAMETERS dataframe of hyperparameter to test in the model
#' @param N_FOLD number of cross validation fold to perform
#' @param MAX_CLUSTER maximum CPU clustering for parallel computing
#' 
#' @return a .pdf in /graphic with the Y and Y_hat per station and variable importance

model_eval <- function(bluecloud.wd = bluecloud_dir,
                       by_target = TRUE,
                       var_importance = FALSE){
  
  # --- Loading data
  setwd(paste0(bluecloud.wd,"/data/"))
  m <- py_load_object("m", pickle = "pickle")
  HYPERPARAMETERS <- read_feather(paste0(bluecloud.wd,"/data/HYPERPARAMETERS.feather"))
  ID <- read_feather(paste0(bluecloud.wd,"/data/Station_FOLD.feather"))
  Y0 <- as.data.frame(read_feather(paste0(bluecloud.wd,"/data/Y.feather")))
  X0 <- as.data.frame(read_feather(paste0(bluecloud.wd,"/data/X.feather")))
  
  # --- Initializing outputs
  y_hat <- r2_tar <- r2cor_tar <- rmse_tar <- NULL
  
  # --- Predicting on the validation data
  for(cv in 1:N_FOLD){
    # --- Loading test data and models
    X_val <- as.data.frame(read_feather(paste0(bluecloud.wd,"/data/", cv, "_X_val.feather")))
    m0 <- m[[cv]][[1]]
    
    # --- Do predictions
    y_hat <- rbind(y_hat, mbtr_predict(model = m0, X_pred = X_val, n_boosts = HYPERPARAMETERS$n_boost))
  } # k-fold cv  loop
  
  # --- Put y_hat_all back with the pre-fold station order
  # We performed a random sorting to create the folds, we therefore need to match Y again before evaluation
  y_hat <- y_hat[order(ID$Station),]
  
  # --- Global model evaluation
  r2 <- calc_rsquared(as.matrix(Y0), y_hat)
  rmse <- sqrt(mean(as.matrix((Y0-y_hat)^2), na.rm=TRUE))
  
  print(paste("--- model multidimensional R-squarred is :", round(r2, 2), "---"))
  print(paste("--- model multidimensional RMSE is :", round(rmse, 2), "---"))
  
  # --- Target by target evaluation
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
  
  # --- Plot predictions vs truth
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

  # --- Calculating variable importance
  if(var_importance == TRUE){
    var_count <- matrix(0, ncol = ncol(X0), nrow=N_FOLD)
    colnames(var_count) <- colnames(X0)
    
    for(cv in 1:N_FOLD){
      m0 <- m[[cv]][[1]]
      n_tree <- length(m0$trees)
      
      for(t in 1:n_tree){
        n_node <- length(m0$trees[[t]]$g$nodes$`_nodes`)
        
        for(n in 1:n_node){
          var_nb <- m0$trees[[t]]$g$nodes$`_nodes`[[n]]$variable+1 #py index start at 0, R at 1
          
          if(!is.null(var_nb)){
            dloss <- m0$trees[[t]]$g$nodes$`_nodes`[[n]]$loss*(-1)
            var_count[cv,var_nb] <- var_count[cv,var_nb]+dloss
          }
        } # node loop
      } # tree loop
    } # fold loop
    
    var_imp <- t(apply(var_count, 1, function(x) (x*100)/sum(x, na.rm = TRUE))) %>% 
      apply(2, mean)
    
    # --- Plotting variable importance
    par(mar = c(6, 4, 1, 1))
    plot(x=seq(1,ncol(X0)), y=var_imp, type = 'h',
         lwd=10, lend = 2, ylim=c(0,100), axes = FALSE,
         ylab="variable importance (%)", xlab = "")
    abline(h=(seq(0,100,10)), lty="dotted", col="black")
    axis(side = 2, at = seq(0,100,10), labels = seq(0,100,10))
    axis(side = 1, at = 1:ncol(X0), labels = colnames(X0), las = 2, 
         cex.axis = if(ncol(X0) > 20) {0.7} else {1})
    box()

  } # var imp


} # end function
