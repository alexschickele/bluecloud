#' TO DO:

source(file = "/home/aschickele/workspace/bluecloud descriptor/00_config.R")

# --- Load files
setwd(paste0(bluecloud.wd,"/data/"))
m <- py_load_object("m", pickle = "pickle")
features <- stack(paste0(bluecloud.wd,"/data/features"))

# --- Build X_pred from raster
X <- as.data.frame(getValues(features))

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




