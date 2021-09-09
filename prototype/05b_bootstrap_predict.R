#' Replaces the 03, 04 and 05 scripts by re-training the model on all data according
#' to a-priori selected hyperparameters with script 03 and 04.
#' Produces bootstrap-based predictions

source(file = "/home/aschickele/workspace/bluecloud descriptor/00_config.R")
HYPERPARAMETERS <- read_feather(paste0(bluecloud.wd,"/data/HYPERPARAMETERS.feather"))

# --- Load files
features <- stack(paste0(bluecloud.wd,"/data/features"))
Y0 <- as.data.frame(read_feather(paste0(bluecloud.wd,"/data/Y.feather")))
X0 <- as.data.frame(read_feather(paste0(bluecloud.wd,"/data/X.feather")))
N <- nrow(X0)
X <- as.data.frame(getValues(features))

# ================ PART 1 : generating predictions with bootstrap ==============
# --- Create bootstrap X_tr and Y_tr
for(b in 1:NBOOTSTRAP){
  id <- sample(seq(1:nrow(X0)), replace = TRUE)
  X_tr <- X0[id,]
  Y_tr <- Y0[id,]
  
  write_feather(X_tr, paste0(bluecloud.wd,"/data/",b,"_X_tr.feather"))
  write_feather(Y_tr, paste0(bluecloud.wd,"/data/",b,"_Y_tr.feather"))
} # bootstrap loop

# --- Re-fitting models on all data
b <- 1:NBOOTSTRAP
m <- mcmapply(FUN=mbtr_fit, 
              path=paste0(bluecloud.wd, "/data/", b),
              hp_id = as.character(b),
              loss_type='mse',
              n_boosts = as.integer(HYPERPARAMETERS$n_boost),
              min_leaf= HYPERPARAMETERS$MEAN_LEAF,
              learning_rate=HYPERPARAMETERS$LEARNING_RATE,
              lambda_weights=0.0001,
              lambda_leaves=0.0001,
              n_q= as.integer(HYPERPARAMETERS$N_Q),
              early_stopping_rounds = as.integer(HYPERPARAMETERS$n_boost),
              SIMPLIFY = FALSE,
              USE.NAMES = FALSE,
              mc.cores = min(c(MAX_CLUSTER, NBOOTSTRAP)))

# --- Reload best_hp models from python because of "previous session invalidity"
m <- list()
for(b in 1:NBOOTSTRAP){
  m0 <- py_load_object(paste0(bluecloud.wd,"/data/",b,"_",b,"_m"), pickle = "pickle")
  m <- append(m, list(m0))
} #b loop

# --- Predicting values
m <- lapply(m, function(x) {x[[1]]})
y_hat <- NULL

b = 1
zz <- mbtr_predict(model = m[[b]], X_pred = X)
zz <- mapply(FUN = mbtr_predict,
                  model = m[b],
                  X_pred = as.data.frame(X),
                  SIMPLIFY = FALSE,
                  USE.NAMES = FALSE)

b <- 1:NBOOTSTRAP
zz <- lapply(m, FUN = mbtr_predict(m, X))
zz <- mcmapply(FUN = mbtr_predict,
                  model = m[b],
                  X_pred = X,
                  SIMPLIFY = FALSE,
                  USE.NAMES = FALSE,
                  mc.cores = min(c(MAX_CLUSTER, N_FOLD*nrow(HYPERPARAMETERS))))

y_hat <- NULL
cat(paste("---", Sys.time(), "// Projecting maps // This may take some time --- \n"))
for (b in 1:NBOOTSTRAP){
  m0 <- m[[b]][[1]]
  y_hat <- abind(y_hat, mbtr_predict(m0, X), along = 3)
}
cat(paste("---", Sys.time(), "// Done --- \n"))

# ==================  PART 2 : plotting predictions ============================
# --- Initializing raster and projections data
r0 <- raster(res=res(features), ext=extent(features))

y_hat_m <- apply(y_hat, c(1,2), mean)
y_hat_sd <- apply(y_hat, c(1,2), sd)

cutx <- seq(0,0.5,0.05)
cuty <- seq(0,1,0.1)

# --- Rasterizing bivariate projections
col_matrix <- colmat(pal = brewer.pal(10, "RdYlBu"), saturation = 0,
                     xlab = "Standard deviation", ylab = "Relative Abundance")

for(i in 1:ncol(y_hat_m)){
  r_m <- setValues(r0, y_hat_m[,i])
  r_sd <- setValues(r0, y_hat_sd[,i])
  tmp <- bivar_map(rasterx = r_sd, rastery = r_m, colormatrix = col_matrix, cutx = cutx, cuty = cuty)
  
  if(i == 1){
    proj <- tmp[[1]]
    col <- list(tmp[[2]])
  } else {
    proj <- stack(proj, tmp[[1]])
    col[[i]] <- tmp[[2]]
  }
} # i target loop
proj <- synchroniseNA(stack(features[[1]], proj))[[-1]]
names(proj) <- paste(1:ncol(y_hat_m))

# --- Plotting projections
par(bg="black", col="white", col.axis = "white", col.lab="white",col.main="white",
    mfrow=c(round(sqrt(ncol(y_hat_m))+1),round(sqrt(ncol(y_hat_m)+1))),
    mar=c(2,2,2,1))

for(i in 1:ncol(y_hat_m)){
  plot(proj[[i]], col = col[[i]], legend = FALSE, ylim=c(-90, 90), main= paste("target n°",i))
  abline(h=c(66, 23, 0, -23, -66), lty = c("dotted","solid","solid","solid","dotted"), lwd = c(1,1,2,1,1), col = "white")
}

par(mar=c(5,5,5,5))
colmat_plot(col_matrix, xlab = "Standard deviation", ylab = "Relative Abundance")
axis(side = 1, at = seq(0,1,0.1), labels = cutx)
axis(side = 2, at = seq(0,1,0.1), labels = cuty)

# --- Plot pair-wise spatial correlations
par(mar = c(2,2,3,5))
proj_cor <- layerStats(proj, 'pearson', na.rm = TRUE)
proj_cor <- raster(proj_cor[[1]], xmn = 0.5, ymn = 0.5, xmx  = ncol(y_hat_m)+0.5, ymx = ncol(y_hat_m)+0.5)
plot(proj_cor, col = brewer.pal(10, "RdBu"), main = "pair-wise correlation")
abline(h = seq(0:ncol(y_hat_m))+0.5, v = seq(0:ncol(y_hat_m))+0.5)

#Final check on the proportions
print(range(apply(y_hat, c(1,3), sum)))

#Pairwise correlation between targets




