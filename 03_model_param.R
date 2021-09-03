#' This script uses X.feather and Y.feather as input for a multivariate
#' boosted regression tree model.
#' The MBTR model is run from python via the reticulate package and the
#' mbtr_function.py implementation
#' 
#' TO DO LIST :
#' - put all the script into a function
#' - check in mbtr_function.py why the loss does not increase after a while...

dev.off()
source(file = "/home/aschickele/workspace/bluecloud descriptor/00_config.R")

# --- Custom functions
kfold <- function(x,n) split(x, cut(seq_along(x), n, labels = FALSE))
ma <- function(x, n = 10){stats::filter(x, rep(1 / n, n), sides = 2)}

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
grid(lwd = 2)
abline(h = seq(0,10,0.5), v = seq(0,NBOOST, NBOOST/20), lty = "dotted")
legend(x=NBOOST-0.1*NBOOST, 10, legend = seq(1:nrow(HYPERPARAMETERS)),
       fill = brewer.pal(nrow(HYPERPARAMETERS), "Spectral"),
       title = "hyp. nb. :")

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
mtext(text = paste("Best set of hyperparameters =", best_hp), side = 3)

# --- Reload best_hp models from python because of "previous session invalidity"
best_m <- list()
for(cv in 1:N_FOLD){
  m0 <- py_load_object(paste0(bluecloud.wd,"/data/",cv,"_",best_hp,"_m"), pickle = "pickle")
  best_m <- append(best_m, list(m0))
} #cv loop

# --- Save only the best models
write_feather(HYPERPARAMETERS[best_hp,], paste0(bluecloud.wd,"/data/HYPERPARAMETERS.feather"))
py_save_object(best_m, paste0(bluecloud.wd,"/data/m"), pickle = "pickle")

# --- Remove temporary files
data_file <- list.files(paste0(bluecloud.wd,"/data/"))
for(hp in 1:nrow(HYPERPARAMETERS)){
  rem_file <- data_file[grep(paste0("_",hp,"_"), data_file)]
  file.remove(rem_file)
}

# --- END