
#' Here we build the dataset of environmental variables using the World Ocean
#' Atlas data. The variables to be defined are VAR and DEPTH.
#' 
#' TO DO LIST :
#' - add the possibility of having the layer following the bottom
#' - create the X matrix by loading tara ocean depth and coordinates : after
#' the zoom with Pavla D.
#' - add primary production data as well as potential data from Emile & Sakina
#' 
#' @param input.wd path to the World Ocean Atlas clone on the complex server
#' @param output.wd path to the bluecloud descriptor file
#' @param VAR list of environmental variables to extract from the ncdf files
#' @param DEPTH dataframe(top, bottom) with the bathymetrical range to consider
#' 
#' @return a raster stack of environmental variables in the output directory
#' @return a plot of the environmental variables extracted
#' 

# ============== PART 1 : get world ocean atlas data ===========================

source(file = "/home/aschickele/workspace/bluecloud descriptor/00_config.R")

# --- Initializing some parameters
VAR <- list.files(input.wd)
var_names <- VAR # initializing names for the future list actually available
env_raw <- NULL # for storing all env variable values

# --- Loading environmental variables
for(i in 1: length(VAR)){
  sub_folder <- list.files(path = paste0(data.wd,"/share/WOA/DATA/",VAR[[i]],"/netcdf/"),
                           pattern = c("all|A5B7"))[1]
  
  setwd(paste0(data.wd,"/share/WOA/DATA/",VAR[[i]],"/netcdf/",sub_folder,"/1.00/"))
  nc_files <- list.files()
  
  # Security if the .nc file is missing :
  if (length(nc_files)==0){
    cat(paste("--- The", VAR[[i]], "folder is empty, no .nc files !---\n",
              "=> Variable has been remove from the list \n"))
    var_names <- var_names[-which(var_names==VAR[i])]
    # .nc file is available, we load it and store it  
  }  else {
    nc <- nc_open(nc_files[grep("00_01", nc_files)])
    
    var_short <- substr(nc_files[1], nchar(nc_files[1])-8, nchar(nc_files[1])-8)
    env_raw <- abind(env_raw, ncvar_get(nc, paste0(var_short,"_an")), along = 4)
  } # end .nc security if
} # end i VAR loop

# --- Selecting the depth range
depth_bnds <- ncvar_get(nc, "depth_bnds")

id_depth <- data.frame(top=head(which(depth_bnds[1,]<=DEPTH$top), n=1),
                       bottom=tail(which(depth_bnds[2,]<=DEPTH$bottom), n=1))

env_data <- apply(env_raw[,,c(id_depth$top:id_depth$bottom),],c(1,2,4),
                  function(x){mean(x,na.rm = TRUE)})

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
writeRaster(r_env, paste0(bluecloud.wd,"/data/features"), overwrite = TRUE)
plot(r_env)


# ================== Part 2 : add CMEMS Chl data to features ===================

source(file = "/home/aschickele/workspace/bluecloud descriptor/00_config.R")

#--- Initializing
sub_folder <- list.files(path = paste0(data.wd,"/share/cmems/"))
sub_folder <- sub_folder[-grep("ACRI",sub_folder)]

start <- as.numeric(substr(sub_folder[1],1,4))
end <- as.numeric(substr(sub_folder[length(sub_folder)],1,4))

# --- Loading data
for(y in start:end){
  cat(paste(Sys.time(), "--- Loading CHL for year :", y, "--- \n"))
  for(m in 1:12){
    nc <- nc_open(paste0(data.wd,"/share/cmems/",y,"_",m,".nc"))
    r <- raster(ncvar_get(nc, "CHL"))
    r <- aggregate(r, fact = 24, fun=function(x, ...){mean(x, na.rm=TRUE)})
    
    if(m==1 & y==start){
      chl_stack <- r
    } else {
      chl_stack <- addLayer(chl_stack, r)
    }
    
  } # month loop
} # year loop

features <- stack(paste0(bluecloud.wd,"/data/features"))

CHL <- t(mean(chl_stack, na.rm=TRUE))
extent(CHL) <- extent(features)

features <- synchroniseNA(addLayer(features, CHL=CHL))
names(features[[nlayers(features)]]) <- "CHL"

writeRaster(features, paste0(bluecloud.wd,"/data/features"), overwrite = TRUE)

# ============== PART 3 : add supplementary rasters to features ================

source(file = "/home/aschickele/workspace/bluecloud descriptor/00_config.R")

# --- Loading data
features <- stack(paste0(bluecloud.wd,"/data/features"))
bathy <- raster(paste0(data.wd,"/data/environmental_data/bathymetry"))

# --- Creating distance raster
dist <- bathy
dist[is.na(dist)] <- 9999
dist[dist<9999] <- NA

dist <- distance(dist)
dist[dist==9999] <- NA

writeRaster(dist, paste0(data.wd,"/data/environmental_data/distcoast"), overwrite = TRUE)

# --- Adding raster to features
features <- synchroniseNA(stack(features, bathy, dist))
names(features[[10]]) <- "bathymetry"
names(features[[11]]) <- "distcoast"

writeRaster(features, paste0(bluecloud.wd,"/data/features"), overwrite = TRUE)

# --- END