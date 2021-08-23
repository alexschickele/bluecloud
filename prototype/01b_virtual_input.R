
#' This is a prototype script that generates input data for the MBTR model
#' in the same format as the future environmental (X) and genomic (Y) data
#' By A. Schickele, 24 june 2021 //

library(feather)

# --- working directory
output.wd <- "~/workspace/bluecloud descriptor/data"

# --- Number of sites
N <- 250

# --- Generating virtual environmental features
temperature <- seq(10,25,length.out =  N)
salinity <- seq(32,36,length.out =  N)
bathymetry <- seq(-100,0,length.out =  N)
range <- seq(10,2,length.out =  N)
pp <- exp(seq(1,5,length.out =  N))

X <- data.frame(temperature, bathymetry, salinity, range, pp)
# X <- data.frame(rev(temperature), rev(salinity), rev(bathymetry))
# X <- data.frame(temperature, salinity, bathymetry)
# X <- data.frame(temperature)

# --- Generating virtual abundance targets
y1<- seq(0,0.8,length.out =  N)
y2<- seq(0,0.2,length.out =  N)
y3<- seq(1,0,length.out =  N)

Y0 <- data.frame(y1, y2, y3)

# --- Inducing noise in the abundance values at several station
# Random amplitude change
noise <- rnorm(N,0.2,0.02)
Y0 <- Y0+noise

# Station 1 to 50 have a +2 biais
# Y0[1:50,] <- Y0[1:50,]+2
# Station 1 to 50 have a *2 biais
# Y0[1:50,] <- Y0[1:50,]*2

# --- Standardizing data ? ---
Y <- Y0 # NO
rel_abs <- function(x) {
  if(sum(x)!=0){
    x <- x/sum(x)
  }
  return(x)
}# end function

Y <- data.frame(t(apply(Y0,1, rel_abs))) # YES

plot(Y[,1], type = 'l', ylim = c(0,max(Y)), main = "Y", col = "red")
lines(Y[,2], col = "yellow")
lines(Y[,3], col = "green")

# --- Standardizing feature data ? ---
# X0 <- X
# X <- data.frame(apply(X0, 2, function(x) {x <- x/max(sqrt(x^2), na.rm = TRUE)}))


# --- Saving data in .feather files
write_feather(X, paste0(output.wd,"/X.feather"))
write_feather(Y, paste0(output.wd,"/Y.feather"))
