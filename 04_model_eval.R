#' TO DO LIST
#' - put script in function

source(file = "/home/aschickele/workspace/bluecloud descriptor/00_config.R")

# --- Loading data
setwd(paste0(bluecloud.wd,"/data/"))
m <- py_load_object("m", pickle = "pickle")
Y0 <- as.data.frame(read_feather(paste0(bluecloud.wd,"/data/Y.feather")))
X0 <- as.data.frame(read_feather(paste0(bluecloud.wd,"/data/X.feather")))
HYPERPARAMETERS <- read_feather(paste0(bluecloud.wd,"/data/HYPERPARAMETERS.feather"))

# --- Plotting relative abundance
pal <- brewer.pal(ncol(Y0), "Spectral")

par(bg="black", col="white", col.axis = "white", col.lab="white",col.main="white")
plot(Y0[,1], type='l', ylim = c(0,1), ylab = "relative abundance", xlab = "obs", col="black")
for(cv in 1:N_FOLD){
  m0 <- m[[(hp-1)*N_FOLD+cv]][[1]]
  y_hat <- mbtr_predict(m0, X0)
  
  for(i in 1:ncol(Y0)){
    lines(y_hat[,i], col=pal[i], lwd=1)
    lines(Y0[,i], lty="dotted", col=pal[i], lwd=2)
    
  } # i target
} # k-fold cv  loop

legend(x=nrow(y_hat)-0.1*nrow(y_hat), 1, legend = seq(1:ncol(y_hat)),
       fill = brewer.pal(ncol(y_hat), "Spectral"),
       title = "tar. nb. :", border="white", box.col = "white")

# --- Calculating accuracy metric
accurracy <- NULL

for(cv in 1:N_FOLD){
  m0 <- m[[(hp-1)*N_FOLD+cv]][[1]]
  y_hat <- mbtr_predict(m0, X0)
  
  accurracy <- c(accurracy,calc_rsquared(as.matrix(Y0), y_hat)*100)
}

cat(paste("--- model multidimensional accurracy (R2) is :", round(mean(accurracy),2), "+/-", round(sd(accurracy),2), "% --- \n"))

# --- Calculating variable importance
var_count <- matrix(0, ncol = ncol(X0), nrow=N_FOLD)
colnames(var_count) <- colnames(X0)

for(cv in 1:N_FOLD){
  m0 <- m[[(hp-1)*N_FOLD+cv]][[1]]
  n_tree <- length(m0$trees)
  
  for(t in 1:n_tree){
    n_node <- length(m0$trees[[t]]$g$nodes$`_nodes`)
    
    for(n in 1:n_node){
      var_nb <- m0$trees[[t]]$g$nodes$`_nodes`[[n]]$variable+1 #py index start at 0, R at 1
      
      if(!is.null(var_nb)){
        var_count[cv,var_nb] <- var_count[cv,var_nb]+1
      }
    } # node loop
  } # tree loop
} # fold loop

var_imp <- t(apply(var_count, 1, function(x) (x*100)/sum(x, na.rm = TRUE)))

# --- Plotting variable importance
pal <- rep(brewer.pal(ncol(X0), "Spectral"), each = N_FOLD)

plot(x=rep(seq(1,ncol(X0)), each = N_FOLD), y=var_imp, col = pal,
     pch=18, cex=2, ylim=c(0,100),
     ylab="variable importance (%)", xlab="variable")
abline(h=(seq(0,100,20)), lty="dotted", col="white")
legend(x=ncol(X0)-0.2*ncol(X0), 100, legend = colnames(X0),
       fill = brewer.pal(ncol(X0), "Spectral"),
       title = "variables :", border="white", box.col = "white")

# END
