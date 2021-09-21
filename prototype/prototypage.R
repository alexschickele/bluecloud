
# ============== PART 1 : get world ocean atlas data with range YOUPI ==========

source(file = "/home/aschickele/workspace/bluecloud descriptor/00_config.R")

# --- Initializing some parameters
VAR <- list.files(paste0(data.wd,"/share/WOA/DATA"))
var_names <- VAR # initializing names for the future list actually available
env_raw <- NULL # for storing all env variable values

if(DEPTH == "SUR"){
  DEPTH <- data.frame(top = 0, bottom = 10) # TO IMPROVE
}

# --- Loading environmental variables
for(i in 1: length(VAR)){
  sub_folder <- list.files(path = paste0(data.wd,"/share/WOA/DATA/",VAR[[i]],"/netcdf/"),
                           pattern = c("all|A5B7"))[1]
  setwd(paste0(data.wd,"/share/WOA/DATA/",VAR[[i]],"/netcdf/",sub_folder,"/1.00/"))
  nc_files <- list.files()
  nc_files <- nc_files[grep(paste(paste0(c(rep("0",9), rep("", 3)), seq(1:12), "_01"), collapse = "|"), nc_files)]
  
  month_raw <- NULL # temporary array for monthly data
  
  # Security if the .nc file is missing :
  if (length(nc_files)!=12){
    cat(paste("--- The", VAR[[i]], "all month are not available !---\n",
              "=> Variable has been remove from the list \n"))
    var_names <- var_names[-which(var_names==VAR[i])]
    # .nc file is available, we load it and store it  
  }  else {
    for(n in 1:12){
      nc <- nc_open(nc_files[n])
      var_short <- substr(nc_files[1], nchar(nc_files[1])-8, nchar(nc_files[1])-8)
      month_raw <- abind(month_raw, ncvar_get(nc, paste0(var_short,"_an")), along = 4)
    } # n month loop
    
    # --- Extract salinity for later MLD calculations
    if(i == which(var_names == "salinity")){salinity_raw <- month_raw}
    # --- Extract temperature for later MLD calculations
    if(i == which(var_names == "temperature")){
      temperature_raw <- month_raw
      depth_raw <- ncvar_get(nc, "depth_bnds") %>% apply(2, mean)
      lat_raw <- lat_bnds <- ncvar_get(nc, "lat_bnds") %>% apply(2, mean)
      lon_raw <- lon_bnds <- ncvar_get(nc, "lon_bnds") %>% apply(2, mean)}
    
    # --- Selecting the depth range
    depth_bnds <- ncvar_get(nc, "depth_bnds")
    id_depth <- data.frame(top=head(which(depth_bnds[1,]<=DEPTH$top), n=1),
                           bottom=tail(which(depth_bnds[2,]<=DEPTH$bottom), n=1))
    month_raw <- apply(month_raw[,,c(id_depth$top:id_depth$bottom),],c(1,2,4),
                       function(x){mean(x,na.rm = TRUE)})
  } # end .nc security if
  env_raw <- abind(env_raw, month_raw, along = 4)
} # end i VAR loop

# --- Calculate yearly mean and range
var_names <- paste0(var_names, rep(c("_mean","_range"), each = length(var_names)))
env_data <- abind(apply(env_raw, c(1,2,4), function(x){mean(x, na.rm = TRUE)}),
                  apply(env_raw, c(1,2,4), function(x){max(x, na.rm = TRUE)-min(x, na.rm = TRUE)}),
                  along = 3)

# --- Calculate MLD
tmp <- array(data = NA, dim = dim(temperature_raw))

lon <- apply(tmp, c(2,3,4), function(x){x <- lon_raw})
lat <- apply(tmp, c(1,3,4), function(x){x <- lat_raw}) %>%
  aperm(c(2,1,3,4))
depth <- apply(tmp, c(1,2,4), function(x){x <- depth_raw}) %>%
  aperm(c(2,3,1,4))

pressure_raw <- mcmapply(FUN = swPressure,
                         depth = depth, 
                         latitude = lat,
                         SIMPLIFY = TRUE,
                         mc.cores = 5)

density_raw <- mapply(FUN = swSigma,
                        salinity = salinity_raw,
                        temperature = temperature_raw,
                        pressure = pressure,
                        longitude = lon,
                        latitude = lat,
                        SIMPLIFY = TRUE)


# --- Creating raster stack
lat_bnds <- ncvar_get(nc, "lat_bnds")
lon_bnds <- ncvar_get(nc, "lon_bnds")

r <- raster(xmn = min(lon_bnds[1,]), xmx = max(lon_bnds[2,]),
            ymn = min(lat_bnds[1,]), ymx = max(lat_bnds[2,]),
            resolution = lon_bnds[1,2]-lon_bnds[1,1])

for (i in 1:length(var_names)){
  if (i==1) {
    r_env <- setValues(r, as.vector(env_data[,,i]))
  } else {
    tmp <- setValues(r, as.vector(env_data[,,i]))
    r_env <- stack(r_env, tmp)
  } #  end if
} # end var_names loop

names(r_env) <- var_names
r_env <- flip(r_env, direction = 'y')

# --- Saving resulting raster
plot(r_env)





# ============== Building target_raw dataset from Pavla real data youpi ========

source(file = "/home/aschickele/workspace/bluecloud descriptor/00_config.R")

# --- Loading data
lonlat <- read.csv(paste0(data.wd,"/data/SMAGs_Env.csv"), sep=';', header = TRUE)
data_cluster <- read_feather(paste0(data.wd,"/data/CC_PFAM_Carb_taxo_80SS.feather"))
data_reads <- read.table(paste0(data.wd,"/data/SMAGs-v1.cds.95.mg.matrix_CARB"))

data_reads <- data_reads[1:25,1:6]

# --- Reshape data_reads in single entry dataframe
data_reads <- cbind(rownames(data_reads), data_reads)
data_reads <- reshape(data_reads, varying = list(2:ncol(data_reads)),
                      idvar = 1, ids = rownames(data_reads), times = colnames(data_reads)[-1],
                      direction = "long")

rownames(data_reads) <- NULL
colnames(data_reads) <- c("Genes","code","Reads")

# --- Decompose code names
station <- substr(data_reads$code, start = 2, stop = 4)
depth <- substr(data_reads$code, start = 5, stop = 7)
filter <- substr(data_reads$code, start = 9, stop = 12)

# --- Merge all
data <- merge(x = data_cluster, y = data_reads,
              by = "Genes")
data <- cbind(data, station, depth, filter)

# --- Linking station number with longitude and latitude
get_station <- function(x){substr(x = x, start = 6, stop = nchar(x)-4)}
lonlat$Station<- get_station(lonlat$Station)

data$lon <- NA
data$lat <- NA

for (i in 1:nrow(data)){
  data$lon[i] <- lonlat$Longitude[which(as.numeric(lonlat$Station)==data$station[i])[1]]
  data$lat[i] <- lonlat$Latitude[which(as.numeric(lonlat$Station)==data$station[i])[1]]
}

#' =============================================================================
#' =============================================================================
#' =============================================================================

setMethod("plot", "zzplot", function(x){
  barplot(x)
})


# RASTER TEST
# Proof for JO, to keep !!!

library(raster)
library(RColorBrewer)
source(file = paste0(input.wd,"/function/bivarRasterPlot.R"))

r1 <- raster(matrix(rep(seq(1:10),10), 10, 10))
r2 <- raster(matrix(rep(seq(1:10),each = 10), 10, 10))

col_matrix <- colmat(pal = brewer.pal(5, "RdYlBu"),
                     saturation = 0,
                     xlab = "Standard deviation",
                     ylab = "Relative Abundance")

par(mfrow=c(2,2))
plot(r1, main="SD", col = brewer.pal(5, "Greys"))
plot(r2, main="Mean", col = brewer.pal(5, "RdYlBu"))
colmat_plot(col_matrix,
            xlab = "Standard deviation ->",
            ylab = "Relative Abundance ->")

proj <- bivar_map(rasterx = r1,
                  rastery = r2,
                  colormatrix = col_matrix,
                  cutx = seq(0,5,1),
                  cuty = seq(0,10,2))

plot(proj[[1]], col = proj[[2]])



# PARALLEL PROTOTYPING ======================

ls()
rm(list=ls())

input.wd <- "/home/aschickele/workspace/bluecloud descriptor"
output.wd <- "/home/aschickele/workspace/bluecloud descriptor"

# --- Loading R packages
library(reticulate)
library(feather)
library(RColorBrewer)
library(foreach)
library(doParallel)

# --- Custom functions
kfold <- function(x,n) split(x, cut(seq_along(x), n, labels = FALSE))
source_python(paste0(input.wd,"/function/mbtr_function.py"))

# --- Load data
X0 <- read_feather(paste0(input.wd,"/data/X.feather"))
Y0 <- read_feather(paste0(input.wd,"/data/Y.feather"))
N <- nrow(X0)

# --- Defining PARAMETERS/ model HYPERPARAMETERS
HYPERPARAMETERS <- data.frame(LEARNING_RATE = c(1e-1, 3e-2,1e-2),
                              N_Q = rev(c(10, 20, 50)),
                              MEAN_LEAF = rev(c(10, 30, 50)))
NBOOST <- 10
N_FOLD <- 3

# --- Initialize k-fold cross validation splits
id <- sample(x = seq(1:N), size = N, replace = FALSE)
FOLDS <- kfold(id,N_FOLD)

# --- Apply stuff
for (cv in 1:N_FOLD){
  X_tr <- as.data.frame(X0[sort(unlist(FOLDS[-cv])),])
  Y_tr <- as.data.frame(Y0[sort(unlist(FOLDS[-cv])),])
  
  write_feather(X_tr, paste0(output.wd,"/data/",cv,"_X_tr.feather"))
  write_feather(Y_tr, paste0(output.wd,"/data/",cv,"_Y_tr.feather"))
  
  X_val <- as.data.frame(X0[sort(FOLDS[[cv]]),])
  Y_val <- as.data.frame(Y0[sort(FOLDS[[cv]]),])
  
  write_feather(X_val, paste0(output.wd,"/data/",cv,"_X_val.feather"))
  write_feather(Y_val, paste0(output.wd,"/data/",cv,"_Y_val.feather"))
}

library(parallel)
cl <- makeCluster(min(c(12, N_FOLD*nrow(HYPERPARAMETERS))))

cv <- rep(seq(1:N_FOLD),nrow(HYPERPARAMETERS))
hp <- rep(seq(1:N_FOLD), each = nrow(HYPERPARAMETERS))

m <- mcmapply(FUN=mbtr_fit, 
              path=paste0(output.wd, "/data/", cv),
              hp_id = as.character(hp),
              loss_type='mse',
              n_boosts = as.integer(NBOOST),
              min_leaf= HYPERPARAMETERS$MEAN_LEAF[hp],
              learning_rate=HYPERPARAMETERS$LEARNING_RATE[hp],
              lambda_weights=0,
              lambda_leaves=0,
              n_q= as.integer(HYPERPARAMETERS$N_Q[hp]),
              val_path = paste0(output.wd,"/data/", cv),
              early_stopping_rounds = as.integer(0.5/HYPERPARAMETERS$LEARNING_RATE[hp]),
              SIMPLIFY = FALSE,
              USE.NAMES = FALSE)

stopCluster(cl)

m <- list()
for (hp in 1:nrow(HYPERPARAMETERS)){
  for(cv in 1:N_FOLD){
    m0 <- py_load_object(paste0(output.wd,"/data/",cv,"_",hp,"_m"), pickle = "pickle")
    m <- append(m, list(m0))
  } #cv loop
} #hp loop

for(hp in 1:nrow(HYPERPARAMETERS)){
  cat(paste("---", Sys.time(), "fit :", toString(names(HYPERPARAMETERS)), toString(HYPERPARAMETERS[hp,])), "--- \n")
  
  # foreach(cv = 1:N_FOLD, .packages=c("feather", "reticulate"), .verbose = TRUE) %dopar% {
  for (cv in 1:N_FOLD){
    if(hp==1){
      # --- Preparing training and validation data
      X_tr <- as.data.frame(X0[sort(unlist(FOLDS[-cv])),])
      Y_tr <- as.data.frame(Y0[sort(unlist(FOLDS[-cv])),])
      
      write_feather(X_tr, paste0(output.wd,"/data/",cv,"_X_tr.feather"))
      write_feather(Y_tr, paste0(output.wd,"/data/",cv,"_Y_tr.feather"))
      
      X_val <- as.data.frame(X0[sort(FOLDS[[cv]]),])
      Y_val <- as.data.frame(Y0[sort(FOLDS[[cv]]),])
      
      write_feather(X_val, paste0(output.wd,"/data/",cv,"_X_val.feather"))
      write_feather(Y_val, paste0(output.wd,"/data/",cv,"_Y_val.feather"))
    }
    
    # --- Fitting the model
    m0 <- mbtr_fit(path=paste0(output.wd, "/data/", cv),
                   loss_type='mse',
                   n_boosts = as.integer(NBOOST),
                   min_leaf= HYPERPARAMETERS$MEAN_LEAF[hp],
                   learning_rate=HYPERPARAMETERS$LEARNING_RATE[hp],
                   lambda_weights=0,
                   lambda_leaves=0,
                   n_q= as.integer(HYPERPARAMETERS$N_Q[hp]),
                   val_path = paste0(output.wd,"/data/", cv),
                   early_stopping_rounds = as.integer(0.5/HYPERPARAMETERS$LEARNING_RATE[hp]))
    
    # m <- append(m, list(m0))
    
    if(hp==nrow(HYPERPARAMETERS)){
      # --- Deleting the train and validation files
      file.remove(paste0(output.wd,"/data/",cv,"_X_tr.feather"))
      file.remove(paste0(output.wd,"/data/",cv,"_Y_tr.feather"))
      file.remove(paste0(output.wd,"/data/",cv,"_X_val.feather"))
      file.remove(paste0(output.wd,"/data/",cv,"_Y_val.feather"))
    }
  } # k-fold cross validation loop
} # hyperparameter loop

# stopCluster(cl)




zz <- list(1,2,3)
vars1<-c(1,2,3)
vars2<-c(10,20,30)
mult_one<-function(var1,var2)
{
  list(a=var1, b=var2)
}
zz <- mapply(FUN=mult_one,var1=vars1,var2=vars2, SIMPLIFY = FALSE)






#' Testing some custom loss function on the y and y_hat

Y <- matrix(data = c(0,0.3,0.7,0.3,0,0.7,0,0.7,0.3),
            nrow = 3,
            ncol = 3)

Y_HAT_good <- Y*10

Y_HAT_bad <- matrix(data = rep((1/3),10),
                     nrow = 3,
                     ncol = 3)

Y_HAT_bad <- matrix(data = c(0,0.5,0.5,0.4,0,0.6,0,0.7,0.3),
                    nrow = 3,
                    ncol = 3)

sqrt(mean((Y - Y_HAT_bad)^2))
sqrt(mean((Y - Y_HAT_good)^2))

sqrt(mean((Y - Y_HAT_bad/(sum(Y_HAT_bad)/sum(Y)))^2))
sqrt(mean((Y - Y_HAT_good/(sum(Y_HAT_good)/sum(Y)))^2)) # works 

# NOW how to calculate gradient and hessian of this ? :D
g <- Y - Y_HAT_bad/(sum(Y_HAT_bad)/sum(Y))
h <- Y/Y




# PROTOTYPE DE FONCTION

git clone https://github.com/zzd1992/GBDTMO.git

from gbdtmo import load_lib, GBDTMulti, GBDTSingle
import numpy as np
from pandas import read_feather

LIB = load_lib("/home/aschickele/workspace/custom package/GBDTMO/build/gbdtmo.so")


inp_dim, out_dim = 10, 5
params = {"max_depth": 5, "lr": 0.1, 'loss': b"mse", "num_threads": 1}
booster = GBDTMulti(LIB, out_dim=out_dim, params=params)

inp_dim, out_dim = 5, 10
x_train, y_train = np.random.rand(10000, inp_dim), np.random.rand(10000, out_dim)
x_valid, y_valid = np.random.rand(10000, inp_dim), np.random.rand(10000, out_dim)


x_train, y_train = np.random.rand(10000, inp_dim), np.random.rand(10000, out_dim)
x_valid, y_valid = np.random.rand(10000, inp_dim), np.random.rand(10000, out_dim)
booster.set_data((x_train, y_train), (x_valid, y_valid))






