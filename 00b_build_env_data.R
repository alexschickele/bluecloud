#' @concept building the dataset of environmental variables
#' @source World Ocean Atlas data, Copernicus Marine Science, 'castr' package
#' 
#' @param data.wd path to the World Ocean Atlas clone on the complex server
#' @param bluecloud.wd path to the bluecloud descriptor file
#' 
#' @return a raster stack of environmental variables in the bluecloud directory

source(file = "/home/aschickele/workspace/bluecloud descriptor/00_config.R")

# --- Initializing some parameters
VAR <- list.files(paste0(data.wd,"/share/WOA/DATA"))
var_names <- VAR # initializing names for the future list actually available
env_raw <- NULL # for storing all env variable values

if(DEPTH == "SUR"){DEPTH <- data.frame(top = 0, bottom = 10)}

# ============== PART 1 : collect world ocean atlas data =======================
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
      lat_raw <- ncvar_get(nc, "lat_bnds") %>% apply(2, mean)
      lon_raw <- ncvar_get(nc, "lon_bnds") %>% apply(2, mean)}
    
    # --- Selecting the depth range
    depth_bnds <- ncvar_get(nc, "depth_bnds")
    id_depth <- data.frame(top=head(which(depth_bnds[1,]<=DEPTH$top), n=1),
                           bottom=tail(which(depth_bnds[2,]<=DEPTH$bottom), n=1))
    month_raw <- apply(month_raw[,,c(id_depth$top:id_depth$bottom),],c(1,2,4),
                       function(x){mean(x,na.rm = TRUE)})
  } # end .nc security if
  env_raw <- abind(env_raw, month_raw, along = 4)
} # end i VAR loop

# ==================== PART 2 : calculate MLD data =============================
# --- Calculate pressure and density
tmp <- array(data = NA, dim = dim(temperature_raw))

lon <- apply(tmp, c(2,3,4), function(x){x <- lon_raw})
lat <- apply(tmp, c(1,3,4), function(x){x <- lat_raw}) %>% aperm(c(2,1,3,4))
depth <- apply(tmp, c(1,2,4), function(x){x <- depth_raw}) %>% aperm(c(2,3,1,4))

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

density_raw <- apply(density_raw, c(1,2,4), interp_depth) %>% aperm(c(2,3,1,4))

# --- Calculate MLD and add to existing data
MLD <- apply(density_raw, c(1,2,4), function(x) {
  mld(x = x, depth = depth_out, ref.depths = 1:10, default.depth = 80,
      n.smooth = 0, k = 2, criteria = c(0.03, 0.01))})
env_raw <- abind(env_raw, MLD, along = 4)

# ============== PART 3 : collect CMEMS Chl data and calculate ZE ==============
#--- Initializing
sub_folder <- list.files(path = paste0(data.wd,"/share/cmems/"))
sub_folder <- sub_folder[-grep("ACRI",sub_folder)]

start <- as.numeric(substr(sub_folder[1],1,4))
end <- as.numeric(substr(sub_folder[length(sub_folder)],1,4))

chl_raw <- NULL

# --- Loading CMEMS Chl data
for(y in start:end){
  cat(paste(Sys.time(), "--- Loading CHL for year :", y, "--- \n"))
  month_raw <- NULL
  for(m in 1:12){
    nc <- nc_open(paste0(data.wd,"/share/cmems/",y,"_",m,".nc"))
    r <- raster(ncvar_get(nc, "CHL")) %>% 
      aggregate(fact = 24, fun=function(x, ...){mean(x, na.rm=TRUE)}) %>% 
      as.matrix()
    month_raw <- abind(month_raw, r, along = 3)
  } # month loop
  chl_raw <- abind(chl_raw, month_raw, along = 4)
} # year loop

# --- Calculate CHL, ZE and add to existing data
CHL <- apply(chl_raw, c(1,2,3), function(x){mean(x, na.rm = TRUE)})
CHL <- CHL[,180:1,]
ZE <- apply(CHL, c(1,2,3), ze_from_surface_chla)

env_raw <- abind(env_raw, CHL, ZE, along = 4)

# ========= PART 4: calculating mean etc... and store it in a raster stack =====
# --- Creating raster stack
r <- raster(xmn = min(lon_raw)-0.5, xmx = max(lon_raw)+0.5,
            ymn = min(lat_raw)-0.5, ymx = max(lat_raw)+0.5,
            resolution = lon_raw[2]-lon_raw[1])
var_names <- c(var_names, "MLD","CHL","ZE")
var_names <- paste0(var_names, rep(c("mean","sd","med","mad","min","max"), each = length(var_names)))

# --- Calculating the metrics
env_data <- abind(apply(env_raw, c(1,2,4), function(x){mean(x, na.rm = TRUE)}),
                  apply(env_raw, c(1,2,4), function(x){sd(x, na.rm = TRUE)}),
                  apply(env_raw, c(1,2,4), function(x){median(x, na.rm = TRUE)}),
                  apply(env_raw, c(1,2,4), function(x){mad(x, na.rm = TRUE)}),
                  apply(env_raw, c(1,2,4), function(x){min(x, na.rm = TRUE)}),
                  apply(env_raw, c(1,2,4), function(x){max(x, na.rm = TRUE)}),
                  along = 3)

# --- Save in a raster
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
plot(r_env)
writeRaster(r_env, paste0(bluecloud.wd,"/data/features"), overwrite = TRUE)

# =============== PART 5 : add bathymetry and distance to coast ================

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
names(features[[nlayers(features)-1]]) <- "bathymetry"
names(features[[nlayers(features)]]) <- "distcoast"

writeRaster(features, paste0(bluecloud.wd,"/data/features"), overwrite = TRUE)

# --- END