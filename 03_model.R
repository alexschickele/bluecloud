
#' This script uses X.feather and Y.feather as input for a multivariate
#' boosted regression tree model.
#' The MBTR model is run from python via the reticulate package
#' See MBTR_installation.Rmd for more information on the custom loss function that
#' was implemented

# .rs.restartR()
toto
ls()
rm(list=ls())

input.wd <- "~/workspace/bluecloud descriptor"
output.wd <- "~/workspace/bluecloud descriptor"

# --- Loading R packages
library(reticulate)
library(feather)

# --- Loading Python modules
np <- import("numpy")
mbtr <- import("mbtr")

# --- Load data
X <- as.matrix(read_feather(paste0(input.wd,"/data/X.feather")))
Y <- as.matrix(read_feather(paste0(input.wd,"/data/Y.feather")))
N <- nrow(X)
CHEATCODE <- 1 #testing if multiplying Y is doing something, sadly YES...
Ntarget <- c(2) # which targets do we take into account among 1, 2 and 3

#' In fact, the square function works opposite when the value is <1. That's why 
#' multiplying by 10 or more works. Now need to find a way to correct that.

Y <- Y*CHEATCODE

# --- Plot X data ?
plot(X[,1], type = 'l', ylim = c(-150,50))
lines(X[,2])
lines(X[,3])

# --- Train and test datasets // TO PUT IN A LOOP
ID_tr <- sample(seq(1,N), 0.7*N)
write.table(ID_tr,paste0(output.wd,"/data/ID_tr.txt"))

X_tr <- as.data.frame(X[ID_tr,])
Y_tr <- as.data.frame(Y[ID_tr,Ntarget])

write_feather(X_tr, paste0(output.wd,"/data/X_tr.feather"))
write_feather(Y_tr, paste0(output.wd,"/data/Y_tr.feather"))

ID_te <- seq(1,N)[-ID_tr]
X_te <- as.data.frame(X[ID_te,])
Y_te <- as.data.frame(Y[ID_te,Ntarget])

write_feather(X_te, paste0(output.wd,"/data/X_te.feather"))
write_feather(Y_te, paste0(output.wd,"/data/Y_te.feather"))

# --- Model fit
source_python(paste0(input.wd,"/03b_model_fit.py"))
m0 <- mbtr_fit(path=input.wd, loss_type='mse', 
               learning_rate=0.1, min_leaf=10, lambda_weights=0.1)

head(m0[[2]])
tail(m0[[2]])

head(Y_te)
tail(Y_te)

# --- Does the lines sum to 1 ?
range(apply(m0[[2]], 1, sum))

# --- Model evaluation
rmse_mbt <- sqrt(mean((as.matrix(Y_te/CHEATCODE) - as.matrix(m0[[2]]/CHEATCODE))^2))
print(rmse_mbt)

plot(Y_te[,1], type = 'l', ylim = c(0,max(Y)), main = "Y_te (dotted) VS Y_hat (full)", col = "red", lty = "dotted", lwd=2)
lines(Y_te[,2], col = "blue", lty = "dotted", lwd=2)
lines(Y_te[,3], col = "green", lty = "dotted", lwd=2)

lines(m0[[2]][,1], col = "red")
lines(m0[[2]][,2], col = "blue")
lines(m0[[2]][,3], col = "green")


# --- END

