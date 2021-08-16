
#' WARNING :
#' - the model predicts on different environmental variables than the ones
#' used for training. It works, therefore the maps are totally wrong... but useful
#' to test the code !
#' 
#' TO DO LIST:
#' - cleanup
#' - put inside a function

ls()
rm(list=ls())

input.wd <- "~/workspace/bluecloud descriptor"
output.wd <- "~/workspace/bluecloud descriptor"

# --- Loading R packages
library(reticulate)
library(feather)
library(raster)
library(rasterVis)
library(virtualspecies)
library(abind)

# --- Parameters
N_FOLD <- 3
hp <- 1

# --- Custom functions
source_python(paste0(input.wd,"/function/mbtr_function.py"))

# --- Load files
m <- py_load_object("m", pickle = "pickle")
features <- stack(paste0(input.wd,"/data/features"))

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

pal <- hsv(h=as.vector(y_hat_m[,1]),
           s=as.vector(1-y_hat_sd[,1]),
           v=rep(1,length(y_hat_m[,1])),
           alpha=rep(1,length(y_hat_m[,1])))

pal <- c("red",rep("black",length(y_hat_m)-1))

r <- setValues(r0, y_hat_m[,1])

proj <- stack(features[[1]], r)
proj <- synchroniseNA(proj)[[-1]]

library(rasterVis)
library(ggplot2)

plot(proj, col=pal)
gplot(proj) +
  geom_tile(aes(fill=factor(value)))



proj <- baseline[[1]]
for (i in 1:ncol(Y_pred)){
  r <- setValues(r0, Y_pred[,i])
  proj <- stack(proj, r)
}
proj <- synchroniseNA(proj)[[-1]]

plot(proj)

writeRaster(proj, paste0(output.wd,"/data/proj"), overwrite = TRUE)

