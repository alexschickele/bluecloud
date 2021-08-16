ls()
rm(list=ls())

input.wd <- "/home/aschickele/workspace/bluecloud descriptor"
output.wd <- "/home/aschickele/workspace/bluecloud descriptor"

# --- Loading R packages
library(reticulate)
library(feather)
library(RColorBrewer)
library(foreach)
library(doParallel)

# --- Custom functions
kfold <- function(x,n) split(x, cut(seq_along(x), n, labels = FALSE))
source_python(paste0(input.wd,"/function/mbtr_function.py"))

# --- Load data
X0 <- read_feather(paste0(input.wd,"/data/X.feather"))
Y0 <- read_feather(paste0(input.wd,"/data/Y.feather"))
N <- nrow(X0)

# --- Defining PARAMETERS/ model HYPERPARAMETERS
HYPERPARAMETERS <- data.frame(LEARNING_RATE = c(1e-1, 3e-2,1e-2),
                              N_Q = rev(c(10, 20, 50)),
                              MEAN_LEAF = rev(c(10, 30, 50)))
NBOOST <- 10
N_FOLD <- 3

# --- Initialize k-fold cross validation splits
id <- sample(x = seq(1:N), size = N, replace = FALSE)
FOLDS <- kfold(id,N_FOLD)

# --- Apply stuff
for (cv in 1:N_FOLD){
  X_tr <- as.data.frame(X0[sort(unlist(FOLDS[-cv])),])
  Y_tr <- as.data.frame(Y0[sort(unlist(FOLDS[-cv])),])
  
  write_feather(X_tr, paste0(output.wd,"/data/",cv,"_X_tr.feather"))
  write_feather(Y_tr, paste0(output.wd,"/data/",cv,"_Y_tr.feather"))
  
  X_val <- as.data.frame(X0[sort(FOLDS[[cv]]),])
  Y_val <- as.data.frame(Y0[sort(FOLDS[[cv]]),])
  
  write_feather(X_val, paste0(output.wd,"/data/",cv,"_X_val.feather"))
  write_feather(Y_val, paste0(output.wd,"/data/",cv,"_Y_val.feather"))
}

library(parallel)
cl <- makeCluster(min(c(12, N_FOLD*nrow(HYPERPARAMETERS))))

cv <- rep(seq(1:N_FOLD),nrow(HYPERPARAMETERS))
hp <- rep(seq(1:N_FOLD), each = nrow(HYPERPARAMETERS))

m <- mcmapply(FUN=mbtr_fit, 
              path=paste0(output.wd, "/data/", cv),
              hp_id = as.character(hp),
              loss_type='mse',
              n_boosts = as.integer(NBOOST),
              min_leaf= HYPERPARAMETERS$MEAN_LEAF[hp],
              learning_rate=HYPERPARAMETERS$LEARNING_RATE[hp],
              lambda_weights=0,
              lambda_leaves=0,
              n_q= as.integer(HYPERPARAMETERS$N_Q[hp]),
              val_path = paste0(output.wd,"/data/", cv),
              early_stopping_rounds = as.integer(0.5/HYPERPARAMETERS$LEARNING_RATE[hp]),
              SIMPLIFY = FALSE,
              USE.NAMES = FALSE)

stopCluster(cl)

m <- list()
for (hp in 1:nrow(HYPERPARAMETERS)){
  for(cv in 1:N_FOLD){
    m0 <- py_load_object(paste0(output.wd,"/data/",cv,"_",hp,"_m"), pickle = "pickle")
    m <- append(m, list(m0))
  } #cv loop
} #hp loop

for(hp in 1:nrow(HYPERPARAMETERS)){
  cat(paste("---", Sys.time(), "fit :", toString(names(HYPERPARAMETERS)), toString(HYPERPARAMETERS[hp,])), "--- \n")
  
  # foreach(cv = 1:N_FOLD, .packages=c("feather", "reticulate"), .verbose = TRUE) %dopar% {
  for (cv in 1:N_FOLD){
    if(hp==1){
      # --- Preparing training and validation data
      X_tr <- as.data.frame(X0[sort(unlist(FOLDS[-cv])),])
      Y_tr <- as.data.frame(Y0[sort(unlist(FOLDS[-cv])),])
      
      write_feather(X_tr, paste0(output.wd,"/data/",cv,"_X_tr.feather"))
      write_feather(Y_tr, paste0(output.wd,"/data/",cv,"_Y_tr.feather"))
      
      X_val <- as.data.frame(X0[sort(FOLDS[[cv]]),])
      Y_val <- as.data.frame(Y0[sort(FOLDS[[cv]]),])
      
      write_feather(X_val, paste0(output.wd,"/data/",cv,"_X_val.feather"))
      write_feather(Y_val, paste0(output.wd,"/data/",cv,"_Y_val.feather"))
    }
    
    # --- Fitting the model
    m0 <- mbtr_fit(path=paste0(output.wd, "/data/", cv),
                   loss_type='mse',
                   n_boosts = as.integer(NBOOST),
                   min_leaf= HYPERPARAMETERS$MEAN_LEAF[hp],
                   learning_rate=HYPERPARAMETERS$LEARNING_RATE[hp],
                   lambda_weights=0,
                   lambda_leaves=0,
                   n_q= as.integer(HYPERPARAMETERS$N_Q[hp]),
                   val_path = paste0(output.wd,"/data/", cv),
                   early_stopping_rounds = as.integer(0.5/HYPERPARAMETERS$LEARNING_RATE[hp]))
    
    # m <- append(m, list(m0))
    
    if(hp==nrow(HYPERPARAMETERS)){
      # --- Deleting the train and validation files
      file.remove(paste0(output.wd,"/data/",cv,"_X_tr.feather"))
      file.remove(paste0(output.wd,"/data/",cv,"_Y_tr.feather"))
      file.remove(paste0(output.wd,"/data/",cv,"_X_val.feather"))
      file.remove(paste0(output.wd,"/data/",cv,"_Y_val.feather"))
    }
  } # k-fold cross validation loop
} # hyperparameter loop

# stopCluster(cl)




zz <- list(1,2,3)
vars1<-c(1,2,3)
vars2<-c(10,20,30)
mult_one<-function(var1,var2)
{
  list(a=var1, b=var2)
}
zz <- mapply(FUN=mult_one,var1=vars1,var2=vars2, SIMPLIFY = FALSE)






#' Testing some custom loss function on the y and y_hat

Y <- matrix(data = c(0,0.3,0.7,0.3,0,0.7,0,0.7,0.3),
            nrow = 3,
            ncol = 3)

Y_HAT_good <- Y*10

Y_HAT_bad <- matrix(data = rep((1/3),10),
                     nrow = 3,
                     ncol = 3)

Y_HAT_bad <- matrix(data = c(0,0.5,0.5,0.4,0,0.6,0,0.7,0.3),
                    nrow = 3,
                    ncol = 3)

sqrt(mean((Y - Y_HAT_bad)^2))
sqrt(mean((Y - Y_HAT_good)^2))

sqrt(mean((Y - Y_HAT_bad/(sum(Y_HAT_bad)/sum(Y)))^2))
sqrt(mean((Y - Y_HAT_good/(sum(Y_HAT_good)/sum(Y)))^2)) # works 

# NOW how to calculate gradient and hessian of this ? :D
g <- Y - Y_HAT_bad/(sum(Y_HAT_bad)/sum(Y))
h <- Y/Y




# PROTOTYPE DE FONCTION

git clone https://github.com/zzd1992/GBDTMO.git

from gbdtmo import load_lib, GBDTMulti, GBDTSingle
import numpy as np
from pandas import read_feather

LIB = load_lib("/home/aschickele/workspace/custom package/GBDTMO/build/gbdtmo.so")


inp_dim, out_dim = 10, 5
params = {"max_depth": 5, "lr": 0.1, 'loss': b"mse", "num_threads": 1}
booster = GBDTMulti(LIB, out_dim=out_dim, params=params)

inp_dim, out_dim = 5, 10
x_train, y_train = np.random.rand(10000, inp_dim), np.random.rand(10000, out_dim)
x_valid, y_valid = np.random.rand(10000, inp_dim), np.random.rand(10000, out_dim)


x_train, y_train = np.random.rand(10000, inp_dim), np.random.rand(10000, out_dim)
x_valid, y_valid = np.random.rand(10000, inp_dim), np.random.rand(10000, out_dim)
booster.set_data((x_train, y_train), (x_valid, y_valid))






