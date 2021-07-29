
#' This script uses X.feather and Y.feather as input for a multivariate
#' boosted regression tree model.
#' The MBTR model is run from python via the reticulate package
#' See MBTR_installation.Rmd for more information on the custom loss function that
#' was implemented
#' 
#' TO DO LIST :
#' - automatically save the 5 fold models from the best_hp to be able to 
#' evaluate the standard deviation in script 04 and project standard deviation
#' in script 05 :: DONE
#' - cleanup the .feather that are save or not to speed up the process and
#' limit the RAM used
#' - cleanup the parameters and the .feather loaded in the MBTR_functions
#' - put all the script into a function

input.wd <- "~/workspace/bluecloud descriptor"
output.wd <- "~/workspace/bluecloud descriptor"

# --- Loading R packages
library(reticulate)
library(feather)

# --- Custom functions
kfold <- function(x,n) split(x, cut(seq_along(x), n, labels = FALSE))
source_python(paste0(input.wd,"/function/mbtr_function.py"))

# --- Load data
X0 <- read_feather(paste0(input.wd,"/data/X.feather"))
Y0 <- read_feather(paste0(input.wd,"/data/Y.feather"))

# --- Defining HYPERPARAMETERS
HYPERPARAMETERS <- expand.grid(LEARNING_RATE = c(0.01, 0.03, 0.06, 0.1, 0.15, 0.2),
                               LAMBDA_WEIGHTS = c(1e-5, 1e-4, 1e-3, 1e-2, 1e-1),
                               LAMBDA_LEAVES = c(1e-5, 1e-4, 1e-3, 1e-2, 1e-1),
                               RMSE_mean = NA,
                               RMSE_sd = NA)

# --- Parameters
N_fold <- 5

# --- Excluding TEST DATASET
# Used for final model testing, once the hyperparameters are selected
ID_te <- sample(seq(1, nrow(X0)), 0.2*nrow(X0))

X_te <- as.data.frame(X0[sort(ID_te),])
Y_te <- as.data.frame(Y0[sort(ID_te),])

write_feather(X_te, paste0(output.wd,"/data/X_te.feather"))
write_feather(Y_te, paste0(output.wd,"/data/Y_te.feather"))

X <- as.data.frame(X0[-ID_te,])
Y <- as.data.frame(Y0[-ID_te,])
N <- nrow(X)

# ======================= 1. HYPERPARAMETER SEARCH =============================
hp_sample <- sample(seq(1,nrow(HYPERPARAMETERS)), 30)
for (hp in hp_sample){
  cat(paste("---", Sys.time(), "computing combination", hp, "/", nrow(HYPERPARAMETERS), "--- \n",
            "   lr =", HYPERPARAMETERS$LEARNING_RATE[hp], "lw =", 
            HYPERPARAMETERS$LAMBDA_WEIGHTS[hp], "l_leaf =",
            HYPERPARAMETERS$LAMBDA_LEAVES[hp], "\n"))
  
  m <- list()
  
  for (cv in 1:N_fold){
    # --- Train and validation datasets
    tmp <- sample(x = seq(1:N), size = N, replace = FALSE)
    FOLDS <- kfold(tmp,N_fold)
    
    X_tr <- as.data.frame(X[sort(unlist(FOLDS[-cv])),])
    Y_tr <- as.data.frame(Y[sort(unlist(FOLDS[-cv])),])
    
    write_feather(X_tr, paste0(output.wd,"/data/X_tr.feather"))
    write_feather(Y_tr, paste0(output.wd,"/data/Y_tr.feather"))
    
    X_val <- as.data.frame(X[sort(FOLDS[[cv]]),])
    Y_val <- as.data.frame(Y[sort(FOLDS[[cv]]),])
    
    write_feather(X_val, paste0(output.wd,"/data/X_val.feather"))
    write_feather(Y_val, paste0(output.wd,"/data/Y_val.feather"))
    
    # --- HYPERPARAMETER search : fitting and evaluating model
    m0 <- mbtr_fit(path=input.wd, loss_type='mse', min_leaf=10, 
                   learning_rate=HYPERPARAMETERS$LEARNING_RATE[hp],
                   lambda_weights=HYPERPARAMETERS$LAMBDA_WEIGHTS[hp],
                   lambda_leaves=HYPERPARAMETERS$LAMBDA_LEAVES[hp])
    
    m <- append(m, m0[[1]])
    
    # --- Model evaluation
    if (cv == 1){
      rmse <- sqrt(mean((as.matrix(Y_val) - as.matrix(m0[[2]]))^2))
    } else {
      rmse <- c(rmse, sqrt(mean((as.matrix(Y_val) - as.matrix(m0[[2]]))^2)))
    }
  } # cv k-fold cross validation loop
  HYPERPARAMETERS$RMSE_mean[hp] <- mean(rmse)
  HYPERPARAMETERS$RMSE_sd[hp] <- sd(rmse)
  
  if(HYPERPARAMETERS$RMSE_mean[hp] < min(HYPERPARAMETERS$RMSE_mean[-hp], na.rm = TRUE)){
    best_m <- m
  }
} # hp hyperparameter loop

py_save_object(best_m, "m", pickle = "pickle")
write_feather(HYPERPARAMETERS, paste0(output.wd,"/data/HYPERPARAMETERS.feather"))





# --- END

