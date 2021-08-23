#' Here we build the relative abundance dataset with the TARA OCEAN data
#' The output are
#' 1. A relative abundance datasset : Ytargets * Nobs
#' 2. The corresponding environmental features : Xfeatures * Nobs
#' 
#' TO DO LIST:
#' - to update for unknowns as soon as I have the data
#' - cleanup useless objects for speed and RAM usage
#' - find a way to make a user choice or a loop after checking numbers of genes
#' and location for each protein family, probably in the master script by separating
#' in two functions with a default loop mode

# ================================== PART 1 ====================================
# Building raw dataset from MATOU and the Tara Ocean locations
# - Should not be necessary with the final data from Pavla

source(file = "/home/aschickele/workspace/bluecloud descriptor/00_config.R")

# --- Loading data
lonlat <- read.csv(paste0(data.wd,"/data/SMAGs_Env.csv"), sep=';', header = TRUE)
data <- read_feather(paste0(data.wd,"/data/FF_metaG_Unknown_annot_subset.feather"))

# --- Linking station number with longitude and latitude
get_station <- function(x){substr(x = x, start = 6, stop = nchar(x)-4)}
lonlat$Station<- get_station(lonlat$Station)

data$lon <- NA
data$lat <- NA

for (i in 1:nrow(data)){
  data$lon[i] <- lonlat$Longitude[which(as.numeric(lonlat$Station)==data$station[i])[1]]
  data$lat[i] <- lonlat$Latitude[which(as.numeric(lonlat$Station)==data$station[i])[1]]
}

# --- Save
write_feather(data, path = paste0(bluecloud.wd,"/data/target_raw.feather"))

# ================================== PART 2 ====================================
# Summary of the data

source(file = "/home/aschickele/workspace/bluecloud descriptor/00_config.R")

# --- Load data
target0 <- read_feather(path = paste0(bluecloud.wd,"/data/target_raw.feather"))
target0$PfamName[is.na(target0$PfamName)] <- "Unknown"

# --- Building data summary
Y_summary <- data.frame(PfamName=character(),
                        filter=character(),
                        depth=character(),
                        nb.geneID=integer(),
                        nb.station=integer())

for(p in levels(as.factor(target0$PfamName))){
  for(f in levels(as.factor(target0$filter))){
    for(d in levels(as.factor(target0$depth))){
      tmp <- target0[which(target0$PfamName==p & target0$filter==f & target0$depth==d),]
      nb.gene <- length(unique(tmp$geneID))
      nb.station <- length(unique(tmp$station))
      
      if(nb.gene>MIN.GENE & nb.station>MIN.STATION){
        Y_summary[nrow(Y_summary)+1,] <- list(p,f,d,nb.gene,nb.station)
      }
    } # depth
  } # filter
} # prot family name

print(Y_summary)

# ================================== PART 3 ====================================
# Selecting which protein family to model and build the X and Y dataset from

source(file = "/home/aschickele/workspace/bluecloud descriptor/00_config.R")

# --- Load data
target0 <- read_feather(path = paste0(bluecloud.wd,"/data/target_raw.feather"))
target0$PfamName[is.na(target0$PfamName)] <- "Unknown"
feature0 <- stack(paste0(bluecloud.wd,"/data/features"))

# --- Selecting subset of the target0 data
# TO DO with the final dataset from Pavla, for now i test on unknown surface and ssuu
target1 <- target0[which(target0$PfamName==PFAM & target0$depth==DEPTH & target0$filter==FILTER),]

# --- Building Y
obs <- unique(target1$station)
tar <- unique(target1$geneID)

Y0 <- matrix(0, ncol = length(tar), nrow = length(obs))
dimnames(Y0) <- list(obs,tar)

for (i in 1:length(obs)){
  for (j in 1:length(tar)){
    tmp <- target1[which(target1$geneID==tar[j]  & target1$station==obs[i]),]
    Y0[i,j] <- sum(tmp$readCount)
  }
}

Y <- as.data.frame(Y0/apply(Y0, 1, sum))

# --- Building X
obs_xy <- unique(data.frame(obs=target1$station, x=target1$lon, y=target1$lat))
X <- as.data.frame(extract(feature0, obs_xy[,-1]))

# --- Saving data
out <- which(is.na(X[,1]))
Y <- Y[-out,]
X <- X[-out,]

write_feather(X, path = paste0(bluecloud.wd,"/data/X.feather"))
write_feather(Y, path = paste0(bluecloud.wd,"/data/Y.feather"))

# --- END