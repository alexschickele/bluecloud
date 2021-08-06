#' This script uses X.feather and Y.feather as input for a multivariate
#' boosted regression tree model.
#' The MBTR model is run from python via the reticulate package and the
#' mbtr_function.py implementation
#' 
#' TO DO LIST :
#' - put all the script into a function
#' - check in mbtr_function.py why the loss does not increase after a while...

ls()
rm(list=ls())

input.wd <- "~/workspace/bluecloud descriptor"
output.wd <- "~/workspace/bluecloud descriptor"

# --- Loading R packages
library(reticulate)
library(feather)
library(RColorBrewer)

# --- Custom functions
kfold <- function(x,n) split(x, cut(seq_along(x), n, labels = FALSE))
source_python(paste0(input.wd,"/function/mbtr_function.py"))

# --- Load data
X0 <- read_feather(paste0(input.wd,"/data/X.feather"))
Y0 <- read_feather(paste0(input.wd,"/data/Y.feather"))
N <- nrow(X0)

# --- Defining PARAMETERS/ model HYPERPARAMETERS
HYPERPARAMETERS <- data.frame(LEARNING_RATE = c(1e-1, 1e-2, 1e-3),
                              N_Q = c(10, 20, 50),
                              MEAN_LEAF = c(3, 5, 10))
NBOOST <- 30
N_FOLD <- 3

# --- Initialize k-fold cross validation splits
id <- sample(x = seq(1:N), size = N, replace = FALSE)
FOLDS <- kfold(id,N_FOLD)
m <- list()

for(hp in 1:nrow(HYPERPARAMETERS)){
  cat(paste("---", Sys.time(), "fit :", toString(names(HYPERPARAMETERS)), toString(HYPERPARAMETERS[hp,])), "--- \n")

  # foreach(cv = 1:N_FOLD, .packages=c("feather")) %dopar% {
  for (cv in 1:N_FOLD){
    # --- Preparing training and validation data
    X_tr <- as.data.frame(X0[sort(unlist(FOLDS[-cv])),])
    Y_tr <- as.data.frame(Y0[sort(unlist(FOLDS[-cv])),])
    
    write_feather(X_tr, paste0(output.wd,"/data/X_tr.feather"))
    write_feather(Y_tr, paste0(output.wd,"/data/Y_tr.feather"))
    
    X_val <- as.data.frame(X0[sort(FOLDS[[cv]]),])
    Y_val <- as.data.frame(Y0[sort(FOLDS[[cv]]),])
    
    write_feather(X_val, paste0(output.wd,"/data/X_val.feather"))
    write_feather(Y_val, paste0(output.wd,"/data/Y_val.feather"))
    
    # --- Fitting the model
    m0 <- mbtr_fit(path=input.wd,
                   n_boosts = as.integer(NBOOST),
                   min_leaf= HYPERPARAMETERS$MEAN_LEAF[hp], 
                   learning_rate=HYPERPARAMETERS$LEARNING_RATE[hp],
                   lambda_weights=0,
                   lambda_leaves=0,
                   n_q= as.integer(HYPERPARAMETERS$N_Q[hp]),
                   val_path = paste0(output.wd,"/data"))
    
    m <- append(m, list(m0))
  } # k-fold cross validation loop
} # hyperparameter loop

# --- Save models
setwd(paste0(output.wd,"/data/"))
py_save_object(m, "m", pickle = "pickle")
write_feather(HYPERPARAMETERS, paste0(output.wd,"/data/HYPERPARAMETERS.feather"))

# --- Plotting results
pal <- rep(brewer.pal(nrow(HYPERPARAMETERS), "Spectral"), each = N_FOLD)

par(bg="black", col="white", col.axis = "white", col.lab="white",col.main="white")
plot(unlist(m[[1]][[2]]), type='l', ylim = c(0,10), ylab = "Loss", xlab = "Number of boost rounds")
for (hp in 1:nrow(HYPERPARAMETERS)){
  for (cv in 1:N_FOLD){
    v <- unlist(m[[(hp-1)*N_FOLD+cv]][[2]])
    lines(v, type='l', lwd = 2,
          col = pal[(hp-1)*N_FOLD+cv])
  } # k-fold cv loop
} # hp hyperparameter loop
legend(x=NBOOST-0.1*NBOOST, 10, legend = seq(1:ncol(HYPERPARAMETERS)),
       fill = brewer.pal(nrow(HYPERPARAMETERS), "Spectral"),
       title = "hyp. nb. :")

# --- END