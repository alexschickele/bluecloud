
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
    if(i == which(VAR == "salinity")){salinity_raw <- month_raw}
    # --- Extract temperature for later MLD calculations
    if(i == which(VAR == "temperature")){
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

# --- Calculate pressure and density
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
                         mc.cores = 10)
pressure_raw <- array(pressure_raw, dim = dim(tmp))
density_raw <- mcmapply(FUN = swSigma,
                        salinity = salinity_raw,
                        temperature = temperature_raw,
                        pressure = pressure_raw,
                        longitude = lon,
                        latitude = lat,
                        SIMPLIFY = TRUE,
                        mc.cores = 3)
density_raw <- array(density_raw, dim = dim(tmp))

# --- Interpolating depth to 2m constant resolution
depth_out <- seq(2,300, by = 2)
interp_depth <- function(z){
  if(length(which(!is.na(z))) >= 2){
    z = approx(x = depth_raw[1:29], y = z[1:29], xout = depth_out)$y
  } else {z = rep(NA, length(depth_out))}
}

temperature_raw <- apply(temperature_raw, c(1,2,4), interp_depth) %>% 
  aperm(c(2,3,1,4))
salinity_raw <- apply(salinity_raw, c(1,2,4), interp_depth) %>% 
  aperm(c(2,3,1,4))
density_raw <- apply(density_raw, c(1,2,4), interp_depth) %>% 
  aperm(c(2,3,1,4))

# --- Calculate MLD
MLD <- apply(density_raw, c(1,2,4), function(x) {
  mld(x = x, depth = depth_out, ref.depths = 1:5, default.depth = 50,
      n.smooth = 0, k = 2, criteria = c(0.03, 0.01))})

# --- Add MLD to variables
var_names <- c(var_names, "MLD_mean", "MLD_range")
env_data <- abind(env_data,
                  apply(MLD, c(1,2), function(x){mean(x, na.rm = TRUE)}),
                  apply(MLD, c(1,2), function(x){max(x, na.rm = TRUE)-min(x, na.rm = TRUE)}),
                  along = 3)

# --- Creating raster stack
r <- raster(xmn = min(lon_raw)-0.5, xmx = max(lon_raw)+0.5,
            ymn = min(lat_raw)-0.5, ymx = max(lat_raw)+0.5,
            resolution = lon_raw[2]-lon_raw[1])

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

ZE <- ze_from_surface_chla(CHL)

features <- synchroniseNA(addLayer(features, CHL_mean=CHL, ZE_mean=ZE))
names(features)[19:20] <- c("CHL_mean","ZE_mean")

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
names(features[[21]]) <- "bathymetry"
names(features[[22]]) <- "distcoast"

writeRaster(features, paste0(bluecloud.wd,"/data/features"), overwrite = TRUE)

# --- END