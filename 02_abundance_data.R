
#' Here we build the relative abundance dataset with the TARA OCEAN data
#' The output are
#' 1. A relative abundance datasset : Ytargets * Nobs
#' 2. The corresponding environmental features : Xfeatures * Nobs
#' 
#' TO DO LIST:
#' - do the script as soon as i have the data from Pavla D.

input.wd <- "~/complex/data"
output.wd <- "~/workspace/bluecloud descriptor"

# --- Loading R packages
library(raster)
library(feather)

# --- Input parameters


# --- Loading data
lonlat <- read.csv(paste0(input.wd,"/SMAGs_Env.csv"), sep=';', header = TRUE)
data <- read_feather(paste0(input.wd,"/FF_metaG_Unknown_annot_subset.feather"))

get_station <- function(x){substr(x = x, start = 6, stop = nchar(x)-4)}
lonlat$Station<- get_station(lonlat$Station)

data$lon <- NA
data$lat <- NA

for (i in 1:nrow(data)){
  data$lon[i] <- lonlat$Longitude[which(as.numeric(lonlat$Station)==data$station[i])[1]]
  data$lat[i] <- lonlat$Latitude[which(as.numeric(lonlat$Station)==data$station[i])[1]]
}

data <- data[which(!is.na(data$PfamName)),]

fam <- levels(as.factor(data$PfamName))

# --- Protein-specific dataset ?
# Test on how many sites it is present and how many genes

