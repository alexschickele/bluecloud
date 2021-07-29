
#' TO DO LIST
#' - plot with mse function of the number of boosting rounds to see potential
#' overfitting
#' - counting the number of times each feature are used in splits to calculate
#' the feature importance in the model
#' - adapt this script and the 03 to test each fold model against the test set
#' in order to have an uncertainty estimation

ls()
rm(list=ls())

input.wd <- "~/workspace/bluecloud descriptor"
output.wd <- "~/workspace/bluecloud descriptor"

# --- Loading R packages
library(reticulate)
library(feather)
library(abind)

# --- Custom functions
source_python(paste0(input.wd,"/function/mbtr_function.py"))

# --- Loading data
# HYPERPARAMETERS <- read_feather(paste0(input.wd,"/data/HYPERPARAMETERS.feather"))
# best_hp <- HYPERPARAMETERS[which(HYPERPARAMETERS$RMSE_mean==min(HYPERPARAMETERS$RMSE_mean, na.rm = TRUE)),]

m <- py_load_object("m", pickle = "pickle")
Y_te <- as.data.frame(read_feather(paste0(input.wd,"/data/Y_te.feather")))
X_te <- as.data.frame(read_feather(paste0(input.wd,"/data/X_te.feather"))) # TEST

# =========================== TEST FINAL MODEL ==============================

y_hat <- NULL
for (cv in 1:length(m)){
  y_hat <- abind(y_hat,mbtr_predict(m[[cv]], X_te), along = 3)
} # cv loop

y_hat_mean <- apply(y_hat, c(1,2), mean)
y_hat_sd <- apply(y_hat, c(1,2), sd)


# --- Does the lines sum to 1 ?
range(apply(y_hat, c(1,3), sum))

# --- Plot

plot(Y_te[,1], type = 'l', ylim = c(0,max(Y_te)), main = "Y_te (dotted) VS Y_hat (full)", col = "red", lty = "dotted", lwd=2)
lines(Y_te[,2], col = "blue", lty = "dotted", lwd=2)
lines(Y_te[,3], col = "green", lty = "dotted", lwd=2)

lines(y_hat_mean[,1], col = "red")
lines(y_hat_mean[,2], col = "blue")
lines(y_hat_mean[,3], col = "green")

lines(y_hat_mean[,1]-y_hat_sd[,1])
lines(y_hat_mean[,2]-y_hat_sd[,2])
lines(y_hat_mean[,3]-y_hat_sd[,3])

lines(y_hat_mean[,1]+y_hat_sd[,1])
lines(y_hat_mean[,2]+y_hat_sd[,2])
lines(y_hat_mean[,3]+y_hat_sd[,3])



