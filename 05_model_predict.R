
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
library(virtualspecies)

# --- Custom functions
source_python(paste0(input.wd,"/function/mbtr_function.py"))

# --- Load files
m <- py_load_object("m", pickle = "pickle")
X <- as.data.frame(read_feather(paste0(input.wd,"/data/X.feather"))) # TEST
baseline <- stack(paste0(input.wd,"/data/baseline"))

# --- Build X_pred from raster
X_baseline <- as.data.frame(getValues(baseline))

# --- Predicting values
cat(paste("---", Sys.time(), "// Projecting maps // This may take some time --- \n"))
Y_pred <- mbtr_predict(model=m, X_pred=X_baseline)
cat(paste("---", Sys.time(), "// Done --- \n"))

# --- Put Y_pred back in a raster stack
r0 <- raster(res=res(baseline), ext=extent(baseline))
proj <- baseline[[1]]
for (i in 1:ncol(Y_pred)){
  r <- setValues(r0, Y_pred[,i])
  proj <- stack(proj, r)
}
proj <- synchroniseNA(proj)[[-1]]

plot(proj)

writeRaster(proj, paste0(output.wd,"/data/proj"), overwrite = TRUE)

