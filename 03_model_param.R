#' This script uses X.feather and Y.feather as input for a multivariate
#' boosted regression tree model.
#' The MBTR model is run from python via the reticulate package and the
#' mbtr_function.py implementation
#' 
#' TO DO LIST :
#' - put all the script into a function
#' - check in mbtr_function.py why the loss does not increase after a while...

source(file = "/home/aschickele/workspace/bluecloud descriptor/00_config.R")

# --- Custom functions
kfold <- function(x,n) split(x, cut(seq_along(x), n, labels = FALSE))

# --- Load data
X0 <- read_feather(paste0(bluecloud.wd,"/data/X.feather"))
Y0 <- read_feather(paste0(bluecloud.wd,"/data/Y.feather"))
N <- nrow(X0)

# --- Initialize k-fold cross validation splits
id <- sample(x = seq(1:N), size = N, replace = FALSE)
FOLDS <- kfold(id,N_FOLD)

for (cv in 1:N_FOLD){
  X_tr <- as.data.frame(X0[sort(unlist(FOLDS[-cv])),])
  Y_tr <- as.data.frame(Y0[sort(unlist(FOLDS[-cv])),])
  
  write_feather(X_tr, paste0(bluecloud.wd,"/data/",cv,"_X_tr.feather"))
  write_feather(Y_tr, paste0(bluecloud.wd,"/data/",cv,"_Y_tr.feather"))
  
  X_val <- as.data.frame(X0[sort(FOLDS[[cv]]),])
  Y_val <- as.data.frame(Y0[sort(FOLDS[[cv]]),])
  
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
              lambda_weights=0,
              lambda_leaves=0,
              n_q= as.integer(HYPERPARAMETERS$N_Q[hp]),
              val_path = paste0(bluecloud.wd,"/data/", cv),
              early_stopping_rounds = as.integer(100),
              SIMPLIFY = FALSE,
              USE.NAMES = FALSE,
              mc.cores = min(c(MAX_CLUSTER, N_FOLD*nrow(HYPERPARAMETERS))))

# --- Load models from Python and save in one object
m <- list()
for (hp in 1:nrow(HYPERPARAMETERS)){
  for(cv in 1:N_FOLD){
    m0 <- py_load_object(paste0(bluecloud.wd,"/data/",cv,"_",hp,"_m"), pickle = "pickle")
    m <- append(m, list(m0))
  } #cv loop
} #hp loop

py_save_object(m, paste0(bluecloud.wd,"/data/m"), pickle = "pickle")
write_feather(HYPERPARAMETERS, paste0(bluecloud.wd,"/data/HYPERPARAMETERS.feather"))

# --- Remove temporary files
data_file <- list.files(paste0(bluecloud.wd,"/data/"))
for(cv in 1:N_FOLD){
  rem_file <- data_file[grep(paste0(cv,"_"), data_file)]
  file.remove(rem_file)
}

# --- Plotting results
pal <- rep(brewer.pal(nrow(HYPERPARAMETERS), "Spectral"), each = N_FOLD)

par(bg="black", col="white", col.axis = "white", col.lab="white",col.main="white")
plot(unlist(m[[1]][[2]]), type='l', ylim = c(0,10), xlim=c(0,NBOOST), ylab = "Loss", xlab = "Number of boost rounds")
for (hp in 1:nrow(HYPERPARAMETERS)){
  for (cv in 1:N_FOLD){
    v <- unlist(m[[(hp-1)*N_FOLD+cv]][[2]])
    lines(v, type='l', lwd = 2,
          col = pal[(hp-1)*N_FOLD+cv])
  } # k-fold cv loop
} # hp hyperparameter loop
legend(x=NBOOST-0.1*NBOOST, 10, legend = seq(1:nrow(HYPERPARAMETERS)),
       fill = brewer.pal(nrow(HYPERPARAMETERS), "Spectral"),
       title = "hyp. nb. :")

# --- END