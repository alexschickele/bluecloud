#' @concept train Multivariate Boosted Regression Tree (MBTR) model and find the
#' best set of hyperparameters
#' 
#' @source MBTR python library
#' 
#' @param bluecloud.wd path to the bluecloud folder containing the master script
#' @param HYPERPARAMETERS dataframe of hyperparameter to test in the model
#' @param verbose TRUE or FALSE to silence or not the output in the console
#' 
#' @return a .mbtr object containing the trained model per hyperparameter and
#' cross validation fold

model_run <- function(bluecloud.wd = bluecloud_dir,
                      HYPERPARAMETERS = data.frame(LEARNING_RATE = c(1e-2, 1e-2, 1e-2, 1e-2),
                                                   N_Q = c(10, 10, 10, 10),
                                                   MEAN_LEAF = c(20, 30, 40, 50)),
                      verbose = TRUE){
  
  # --- 1. Custom functions
  # 1.1. Function to split a vector (e.g. station ID) in equal folds
  kfold <- function(x,n) split(x, cut(seq_along(x), n, labels = FALSE))

  # --- 2. Load data
  X0 <- read_feather(paste0(bluecloud.wd,"/data/X.feather"))
  N <- nrow(X0)
  Y0 <- read_feather(paste0(bluecloud.wd,"/data/Y.feather"))

  # --- 3. Perform k-fold cross validation splits
  # 3.1. Initialize the station names
  ID <- read_feather(paste0(bluecloud.wd,"/data/Station_ID.feather"))
  names(ID) <- c("Station","Longitude","Latitude")
  # 3.2. Perform splits according to the station number names
  # One can also do random splits by randomizing the order of "id"
  id <- 1:N
  FOLDS <- kfold(id,N_FOLD)

  # 3.3. Save the new station ID order
  write_feather(ID[id,1], path = paste0(bluecloud.wd,"/data/Station_FOLD.feather"))
  
  # 3.4. Create the different train and validation sets and save them in the data folder
  for (cv in 1:N_FOLD){
    X_tr <- as.data.frame(X0[(unlist(FOLDS[-cv])),])
    Y_tr <- as.data.frame(Y0[(unlist(FOLDS[-cv])),])
    
    write_feather(X_tr, paste0(bluecloud.wd,"/data/",cv,"_X_tr.feather"))
    write_feather(Y_tr, paste0(bluecloud.wd,"/data/",cv,"_Y_tr.feather"))
    
    X_val <- as.data.frame(X0[(FOLDS[[cv]]),])
    Y_val <- as.data.frame(Y0[(FOLDS[[cv]]),])
    
    write_feather(X_val, paste0(bluecloud.wd,"/data/",cv,"_X_val.feather"))
    write_feather(Y_val, paste0(bluecloud.wd,"/data/",cv,"_Y_val.feather"))
  }
  
  # --- 4. Parallel model fitting
  # 4.1. Load the cross validation and hyperparameter order
  cv <- rep(seq(1:N_FOLD),nrow(HYPERPARAMETERS))
  hp <- rep(seq(1:nrow(HYPERPARAMETERS)), each = N_FOLD)
  
  # 4.2. Fit a model to each cv x hp combination
  m <- mcmapply(FUN=mbtr_fit, 
                path=paste0(bluecloud.wd, "/data/", cv),
                hp_id = as.character(hp),
                loss_type='mse',
                n_boosts = as.integer(NBOOST),
                min_leaf= HYPERPARAMETERS$MEAN_LEAF[hp],
                learning_rate=HYPERPARAMETERS$LEARNING_RATE[hp],
                lambda_weights=HYPERPARAMETERS$LEARNING_RATE[hp]/100,
                lambda_leaves=0,
                n_q= as.integer(HYPERPARAMETERS$N_Q[hp]),
                val_path = paste0(bluecloud.wd,"/data/", cv),
                early_stopping_rounds = as.integer(max(50,(1/HYPERPARAMETERS$LEARNING_RATE[hp]))),
                SIMPLIFY = FALSE,
                USE.NAMES = FALSE,
                mc.cores = min(c(MAX_CLUSTER, N_FOLD*nrow(HYPERPARAMETERS))))
  
  # --- 5. Hyperparameter selection
  # 5.1. Plotting the loss per boosting round and cross validation split
  if(verbose == TRUE){
    pal <- rep(brewer.pal(nrow(HYPERPARAMETERS), "Spectral"), each = N_FOLD)
    plot(unlist(m[[1]][[2]]), type='l', ylim = c(0,0.2), xlim=c(0,NBOOST), ylab = "Loss", xlab = "Number of boost rounds")
    for (hp in 1:nrow(HYPERPARAMETERS)){
      for (cv in 1:N_FOLD){
        v <- unlist(m[[(hp-1)*N_FOLD+cv]][[2]])
        v <- v/nrow(Y_tr)
        lines(v, type='l', lwd = 2,
              col = pal[(hp-1)*N_FOLD+cv])
      } # k-fold cv loop
    } # hp hyperparameter loop
    grid(lwd = 2)
    abline(h = seq(0,0.2,0.01), v = seq(0,NBOOST, NBOOST/20), lty = "dotted")
    legend(x=NBOOST-0.1*NBOOST, 0.2, legend = seq(1:nrow(HYPERPARAMETERS)),
           fill = brewer.pal(nrow(HYPERPARAMETERS), "Spectral"),
           title = "hyp. nb. :")
  }

  # 5.2. Select the best hyperparameter
  # By detecting the minimum average loss across cross-validation runs
  losses <- lapply(m, function(x) {x[[2]]})
  HYPERPARAMETERS <- cbind(HYPERPARAMETERS,
                       min_loss = NA,
                       n_boost = NA)
  
  for(hp in 1:nrow(HYPERPARAMETERS)){
    extract_loss <- NULL
    for(cv in 1:N_FOLD){
      tmp <- c(unlist(losses[[(hp-1)*N_FOLD+cv]]), rep(NA, NBOOST-length(losses[[(hp-1)*N_FOLD+cv]])))
      extract_loss <- cbind(extract_loss,tmp)
      extract_loss <- apply(extract_loss, 1, mean)
    }
    HYPERPARAMETERS$min_loss[hp] <- min(extract_loss, na.rm = TRUE)
    HYPERPARAMETERS$n_boost[hp] <- which(extract_loss==min(extract_loss, na.rm = TRUE))
  }
  
  # 5.3. Select the best hyperparameter model object
  best_hp <- which(HYPERPARAMETERS$min_loss==min(HYPERPARAMETERS$min_loss))
  
  if(verbose == TRUE){
    mtext(text = paste("Best set of hyperparameters =", best_hp), side = 3)
    print(paste("Best set of hyperparameters is :", best_hp))
  }

  # 5.4. Reload best_hp models from python to avoid error "previous session invalidity"
  best_m <- list()
  for(cv in 1:N_FOLD){
    m0 <- py_load_object(paste0(bluecloud.wd,"/data/",cv,"_",best_hp,"_m"), pickle = "pickle")
    best_m <- append(best_m, list(m0))
  } #cv loop
  
  # 5.5. Save the corresponding object
  write_feather(HYPERPARAMETERS[best_hp,], paste0(bluecloud.wd,"/data/HYPERPARAMETERS.feather"))
  py_save_object(best_m, paste0(bluecloud.wd,"/data/m"), pickle = "pickle")
  
} # end function
