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

model_eval <- function(bluecloud.wd = "/home/jovyan/bluecloud",
                       by_target = FALSE,
                       var_importance = FALSE){
  
  # --- Loading data
  setwd(paste0(bluecloud.wd,"/data/"))
  m <- py_load_object("m", pickle = "pickle")
  Y0 <- as.data.frame(read_feather(paste0(bluecloud.wd,"/data/Y.feather")))
  X0 <- as.data.frame(read_feather(paste0(bluecloud.wd,"/data/X.feather")))
  
  # --- Initializing outputs
  r2 <- mse <- NULL
  r2_tar <- mse_tar <- matrix(NA, ncol = ncol(Y0), nrow = N_FOLD)
  
  pal <- brewer.pal(ncol(Y0), "Spectral")
  par(bg="black", col="white", col.axis = "white", col.lab="white",col.main="white",
      mfrow =c(ceiling(N_FOLD^0.5),ceiling(N_FOLD^0.5)), mar = c(2,5,2,1))
  
  # --- Evaluating model and plotting relative abundance
  for(cv in 1:N_FOLD){
    # --- Loading test data and models
    X_val <- as.data.frame(read_feather(paste0(bluecloud.wd,"/data/", cv, "_X_val.feather")))
    Y_val <- as.data.frame(read_feather(paste0(bluecloud.wd,"/data/", cv, "_Y_val.feather")))
    m0 <- m[[cv]][[1]]
    
    # --- Do predictions on test set
    y_hat <- mbtr_predict(m0, X_val)
    
    # --- Evaluate model fit
    r2 <- c(r2,calc_rsquared(as.matrix(Y_val), y_hat))
    se <- (Y_val-y_hat)^2
    mse <- c(mse, mean(as.matrix(se), na.rm=TRUE))
    
    # --- Evaluate individual target fit
    for(t in 1:ncol(Y0)){
      r2_tar[cv,t] <- calc_rsquared(as.matrix(Y_val)[,t], y_hat[,t])
      mse_tar[cv,t] <- mean(as.matrix(Y_val[,t]-y_hat[,t])^2, na.rm = TRUE)
    }
  
    # --- Plot prediction against test set
    plot(Y_val[,1], type='l', ylim = c(0,1), ylab = "relative abundance", xlab = "obs", col="black",
         main = paste("fold n°", cv))
    for(i in 1:ncol(Y_val)){
      lines(y_hat[,i], col=pal[i], lwd=1)
      lines(Y_val[,i], lty="dotted", col=pal[i], lwd=2)
    } # i target
    legend(x=nrow(y_hat)-0.2*nrow(y_hat), 1, legend = seq(1:ncol(y_hat)),
           fill = brewer.pal(ncol(y_hat), "Spectral"),
           title = "tar. nb. :", border="white", box.col = "white")
  } # k-fold cv  loop
  
  # --- Print model fit ---
  print(paste("--- model multidimensional R-squarred is :", round(mean(r2),2), "+/-", round(sd(r2),2), "---"))
  print(paste("--- model multidimensional MSE is :", round(mean(mse),2), "+/-", round(sd(mse),2), "---"))
  print(paste("--- model multidimensional RMSE is :", round(mean(sqrt(mse)),2), "+/-", round(sd(sqrt(mse)),2), "---"))
  
  # --- Print individual target fit
  if(by_target == TRUE){
    cat(paste("--- R2 for target", 1:ncol(Y0),":", round(apply(r2_tar, 2, mean),2),"+/-", round(apply(r2_tar, 2, sd),2), "--- \n"))
    cat(paste("--- RMSE for target", 1:ncol(Y0),":", round(apply(sqrt(mse_tar), 2, mean),2),"+/-", round(apply(sqrt(mse_tar), 2, sd),2), "--- \n"))
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
            var_count[cv,var_nb] <- var_count[cv,var_nb]+1
          }
        } # node loop
      } # tree loop
    } # fold loop
    
    var_imp <- t(apply(var_count, 1, function(x) (x*100)/sum(x, na.rm = TRUE)))
    
    # --- Plotting variable importance
    pal <- rep(brewer.pal(ncol(X0), "Spectral"), each = N_FOLD)
    
    plot(x=rep(seq(1,ncol(X0)), each = N_FOLD), y=var_imp, col = pal,
         pch=18, cex=2, ylim=c(0,100),
         ylab="variable importance (%)", xlab="variable")
    abline(h=(seq(0,100,20)), lty="dotted", col="white")
    legend(x=ncol(X0)-0.2*ncol(X0), 100, legend = colnames(X0),
           fill = brewer.pal(ncol(X0), "Spectral"),
           title = "variables :", border="white", box.col = "white")
  }

  
  # --- Clean up temporary files
  # data_file <- list.files(paste0(bluecloud.wd,"/data/"))
  # for(cv in 1:9){
  #   rem_file <- data_file[grep(paste0(cv,"_"), data_file)]
  #   file.remove(rem_file)
  # }

} # end function
