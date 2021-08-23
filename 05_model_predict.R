
#' WARNING :
#' - the model predicts on different environmental variables than the ones
#' used for training. It works, therefore the maps are totally wrong... but useful
#' to test the code !
#' 
#' TO DO LIST:
#' - cleanup
#' - put inside a function

source(file = "/home/aschickele/workspace/bluecloud descriptor/00_config.R")

# --- Load files
m <- py_load_object("m", pickle = "pickle")
features <- stack(paste0(bluecloud.wd,"/data/features"))

# --- Build X_pred from raster
X <- as.data.frame(getValues(features))

# --- Predicting values
y_hat <- NULL

cat(paste("---", Sys.time(), "// Projecting maps // This may take some time --- \n"))
for (cv in 1:N_FOLD){
  m0 <- m[[(hp-1)*N_FOLD+cv]][[1]]
  y_hat <- abind(y_hat, mbtr_predict(m0, X), along = 3)
}
cat(paste("---", Sys.time(), "// Done --- \n"))

# --- Put Y_pred back in a raster stack
r0 <- raster(res=res(features), ext=extent(features))

y_hat_m <- apply(y_hat, c(1,2), mean)
y_hat_sd <- apply(y_hat, c(1,2), sd)

r_m <- setValues(r0, y_hat_m[,1])
r_sd <- setValues(r0, y_hat_sd[,1])

col_matrix <- colmat(pal = brewer.pal(10, "RdYlBu"),
                     saturation = 0,
                     xlab = "Standard deviation",
                     ylab = "Relative Abundance")

proj <- bivar_map(rasterx = r_sd,
                  rastery = r_m,
                  colormatrix = col_matrix,
                  cutx = seq(0,0.1,0.01),
                  cuty = seq(0,1,0.1))

colmat_plot(col_matrix,
            xlab = "Standard deviation",
            ylab = "Relative Abundance")

plot(proj[[1]], col = proj[[2]])



