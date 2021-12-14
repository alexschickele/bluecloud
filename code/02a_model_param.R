#' @concept train Multivariate Boosted Regression Tree (MBTR) model and find the
#' best set of hyperparameters
#' 
#' @source MBTR python library
#' 
#' @param bluecloud.wd path to the bluecloud descriptor file
#' @param HYPERPARAMETERS dataframe of hyperparameter to test in the model
#' @param N_FOLD number of cross validation fold to perform
#' @param NBOOST number of maximum boosting round to perform
#' @param MAX_CLUSTER maximum CPU clustering for parallel computing
#' 
#' @return a .mbtr object containing the trained model per hyperparameter and
#' cross validation fold

model_run <- function(bluecloud.wd = bluecloud_dir,
                      HYPERPARAMETERS = data.frame(LEARNING_RATE = c(1e-2, 1e-2, 1e-2, 1e-2),
                                                   N_Q = c(10, 10, 10, 10),
                                                   MEAN_LEAF = c(20, 30, 40, 50)),
                      verbose = TRUE){
  
  # --- Custom functions
  kfold <- function(x,n) split(x, cut(seq_along(x), n, labels = FALSE))
  ma <- function(x, n = 10){stats::filter(x, rep(1 / n, n), sides = 2)}
  
  # --- Load data
  X0 <- read_feather(paste0(bluecloud.wd,"/data/X.feather"))
  N <- nrow(X0)
  Y0 <- read_feather(paste0(bluecloud.wd,"/data/Y.feather"))

  # --- Initialize k-fold cross validation splits
  ID <- read_feather(paste0(bluecloud.wd,"/data/Station_ID.feather"))
  id <- sample(x = seq(1:N), size = N, replace = FALSE)
  FOLDS <- kfold(id,N_FOLD)
  write_feather(ID[id,], path = paste0(bluecloud.wd,"/data/Station_FOLD.feather"))
  
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
  
  # --- Parallel model fitting
  cv <- rep(seq(1:N_FOLD),nrow(HYPERPARAMETERS))
  hp <- rep(seq(1:nrow(HYPERPARAMETERS)), each = N_FOLD)
  
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
  
  # --- Plotting results
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

  # --- Detect minimum loss and save corresponding model in one object
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
  best_hp <- which(HYPERPARAMETERS$min_loss==min(HYPERPARAMETERS$min_loss))
  
  if(verbose == TRUE){
    mtext(text = paste("Best set of hyperparameters =", best_hp), side = 3)
    print(paste("Best set of hyperparameters is :", best_hp))
  }

  # --- Reload best_hp models from python because of "previous session invalidity"
  best_m <- list()
  for(cv in 1:N_FOLD){
    m0 <- py_load_object(paste0(bluecloud.wd,"/data/",cv,"_",best_hp,"_m"), pickle = "pickle")
    best_m <- append(best_m, list(m0))
  } #cv loop
  
  # --- Save only the best models
  write_feather(HYPERPARAMETERS[best_hp,], paste0(bluecloud.wd,"/data/HYPERPARAMETERS.feather"))
  py_save_object(best_m, paste0(bluecloud.wd,"/data/m"), pickle = "pickle")
  
  # # --- Remove temporary files
  # data_file <- list.files(paste0(bluecloud.wd,"/data/"))
  # for(hp in 1:nrow(HYPERPARAMETERS)){
  #   rem_file <- data_file[grep(paste0("_",hp,"_"), data_file)]
  #   file.remove(rem_file)
  # }
  
} # end function
