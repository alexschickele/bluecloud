
#' WARNING :
#' - the model predicts on different environmental variables than the ones
#' used for training. It works, therefore the maps are totally wrong... but useful
#' to test the code !
#' 
#' TO DO LIST:
#' - cleanup
#' - put inside a function

source(file = "/home/aschickele/workspace/bluecloud descriptor/00_config.R")
NBOOTSTRAP <- 5
HYPERPARAMETERS <- read_feather(paste0(bluecloud.wd,"/data/HYPERPARAMETERS.feather"))

# --- Load files
features <- stack(paste0(bluecloud.wd,"/data/features"))
Y0 <- as.data.frame(read_feather(paste0(bluecloud.wd,"/data/Y.feather")))
X0 <- as.data.frame(read_feather(paste0(bluecloud.wd,"/data/X.feather")))
N <- nrow(X0)

# --- Build X_pred from raster
X <- as.data.frame(getValues(features))

# --- Train model on all data
for(b in 1:NBOOTSTRAP){
  id <- sample(seq(1:nrow(X0)), replace = TRUE)
  X_tr <- X0[id,]
  Y_tr <- Y0[id,]
  
  write_feather(X_tr, paste0(bluecloud.wd,"/data/",b,"_X_tr.feather"))
  write_feather(Y_tr, paste0(bluecloud.wd,"/data/",b,"_Y_tr.feather"))
  
} # bootstrap loop
b <- 1:NBOOTSTRAP
m <- mcmapply(FUN=mbtr_fit, 
              path=paste0(bluecloud.wd, "/data/", b),
              hp_id = as.character(b),
              loss_type='mse',
              n_boosts = as.integer(HYPERPARAMETERS$n_boost),
              min_leaf= HYPERPARAMETERS$MEAN_LEAF,
              learning_rate=HYPERPARAMETERS$LEARNING_RATE,
              lambda_weights=0,
              lambda_leaves=0,
              n_q= as.integer(HYPERPARAMETERS$N_Q),
              val_path = NULL,
              early_stopping_rounds = as.integer(HYPERPARAMETERS$n_boost),
              SIMPLIFY = FALSE,
              USE.NAMES = FALSE,
              mc.cores = min(c(MAX_CLUSTER, NBOOTSTRAP)))

# --- Predicting values
y_hat <- NULL



cat(paste("---", Sys.time(), "// Projecting maps // This may take some time --- \n"))
for (cv in 1:N_FOLD){
  m0 <- m[[cv]][[1]]
  y_hat <- abind(y_hat, mbtr_predict(m0, X), along = 3)
}
cat(paste("---", Sys.time(), "// Done --- \n"))

# --- Put Y_pred back in a raster stack
r0 <- raster(res=res(features), ext=extent(features))

y_hat_m <- apply(y_hat, c(1,2), mean)
y_hat_sd <- apply(y_hat, c(1,2), sd)

# --- Plotting results
plot.new()
par(bg="black", col="white", col.axis = "white", col.lab="white",col.main="white",
    mfrow=c(round(sqrt(ncol(y_hat_m))+1),round(sqrt(ncol(y_hat_m)+1))),
    mar=c(2,2,1,1))

col_matrix <- colmat(pal = brewer.pal(10, "RdYlBu"), saturation = 0,
                     xlab = "Standard deviation", ylab = "Relative Abundance")

for(i in 1:ncol(y_hat_m)){
  r_m <- setValues(r0, y_hat_m[,i])
  r_sd <- setValues(r0, y_hat_sd[,i])
  
  cutx <- seq(0,0.5,0.05)
  cuty <- seq(0,1,0.1)
  
  proj <- bivar_map(rasterx = r_sd, rastery = r_m, colormatrix = col_matrix, cutx = cutx, cuty = cuty)
  
  col <- proj[[2]]
  proj <- synchroniseNA(stack(features[[1]], proj[[1]]))[[2]]

  plot(proj, col = col, legend = FALSE, ylim=c(-90, 90))
  abline(h=c(66, 23, 0, -23, -66), lty = c("dotted","solid","solid","solid","dotted"), lwd = c(1,1,2,1,1), col = "white")
}

par(mar=c(5,8,5,8))
colmat_plot(col_matrix, xlab = "Standard deviation", ylab = "Relative Abundance")
axis(side = 1, at = seq(0,1,0.1), labels = cutx)
axis(side = 2, at = seq(0,1,0.1), labels = cuty)

#Final check on the proportions
print(range(apply(y_hat, c(1,3), sum)))




