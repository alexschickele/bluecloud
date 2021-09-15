#' Here we build the relative abundance dataset with the TARA OCEAN data
#' The output are
#' 1. A relative abundance datasset : Ytargets * Nobs
#' 2. The corresponding environmental features : Xfeatures * Nobs
#' 
#' TO DO LIST:



# ================================== PART 1a ====================================
# Building raw dataset from MATOU and the Tara Ocean locations
# Adapted for the final data

source(file = "/home/aschickele/workspace/bluecloud descriptor/00_config.R")

# --- Loading data
lonlat <- read.csv(paste0(data.wd,"/data/SMAGs_Env.csv"), sep=';', header = TRUE)
data_cluster <- read_feather(paste0(data.wd,"/data/CC_PFAM_Carb_taxo_80SS.feather"))
# data_depth <- read.table(paste0(data.wd,"/data/SMAGS_insitu_WOA.txt"), sep = "\t", header = TRUE)
# data_reads <- read.table(paste0(data.wd,"/data/SMAGs-v1.cds.95.mg.matrix_CARB"))
data_reads <- read.table(paste0(data.wd,"/data/SMAGs-v1.cds.95.mt.matrix_CARB"))

# --- Reshape data_reads in single entry dataframe
data_reads <- cbind(rownames(data_reads), data_reads)
data_reads <- reshape(data_reads, varying = list(2:ncol(data_reads)),
                      idvar = 1, ids = rownames(data_reads), times = colnames(data_reads)[-1],
                      direction = "long")

rownames(data_reads) <- NULL
colnames(data_reads) <- c("Genes","code","readCount")

data_reads <- data_reads[grep("DCM|SUR", data_reads$code),] # keeping only used ones
data_reads <- data_reads[grep("SSUU|MMQQ|QQSS|GGMM", data_reads$code),] # keeping only used ones

# --- Decompose code names
Station <- substr(data_reads$code, start = 2, stop = regexpr("DCM|SUR", data_reads$code)-1)
Station <- sprintf("%03d",as.numeric(Station))
depth <- substr(data_reads$code, start = regexpr("DCM|SUR", data_reads$code), stop = regexpr("DCM|SUR", data_reads$code)+2)
filter <- substr(data_reads$code, start = nchar(data_reads$code)-5, stop = nchar(data_reads$code)-2)

# --- Re-defining station names
get_station <- function(x){substr(x = x, start = 6, stop = nchar(x)-4)}
lonlat$Station<- get_station(lonlat$Station)
lonlat <- unique(lonlat[,c(1,3,4)])

# --- Merge all
data <- cbind(data_reads, Station, depth, filter)
data <- merge(x = data, y = lonlat,
              by = "Station")
data <- merge(x = data, y = data_cluster,
              by = "Genes")

# --- Save
write_feather(data, path = paste0(bluecloud.wd,"/data/target_raw.feather"))

# ================================== PART 1b ====================================
# Building raw dataset from MATOU and the Tara Ocean locations
# - Currently not in use

source(file = "/home/aschickele/workspace/bluecloud descriptor/00_config.R")

# --- Loading data
lonlat <- read.csv(paste0(data.wd,"/data/SMAGs_Env.csv"), sep=';', header = TRUE)
data <- read_feather(paste0(data.wd,"/data/FF_metaT_Unknown_annot_subset.feather"))

# --- Linking station number with longitude and latitude
get_station <- function(x){substr(x = x, start = 6, stop = nchar(x)-4)}
lonlat$Station<- get_station(lonlat$Station)

data$lon <- NA
data$lat <- NA

for (i in 1:nrow(data)){
  data$lon[i] <- lonlat$Longitude[which(as.numeric(lonlat$Station)==data$station[i])[1]]
  data$lat[i] <- lonlat$Latitude[which(as.numeric(lonlat$Station)==data$station[i])[1]]
}
names(data) <- c("Genes", "geneLength.x", "readCount", "couvCum", "couveBase",
                 "sampleName", "Station", "depth", "iteration", "filter", "CC_ID",
                 "PfamDesc", "taxId", "taxName", "taxRank", "taxLineage", "Longitude", "Latitude")

# --- Save
write_feather(data, path = paste0(bluecloud.wd,"/data/target_raw.feather"))

# ================================== PART 2 ====================================
# Summary of the data

source(file = "/home/aschickele/workspace/bluecloud descriptor/00_config.R")

# --- Load data
target0 <- read_feather(path = paste0(bluecloud.wd,"/data/target_raw.feather"))
target0$PFAMs[is.na(target0$PFAMs)] <- "Unknown"

# --- Building data summary
Y_summary <- data.frame(CC_ID=character(),
                        filter=character(),
                        depth=character(),
                        nb.geneID=integer(),
                        nb.station=integer(),
                        sum.reads=integer())

for(c in levels(as.factor(target0$CC_ID))){
  for(f in levels(as.factor(target0$filter))){
    for(d in levels(as.factor(target0$depth))){
      tmp <- target0[which(target0$CC_ID==c & target0$filter==f & target0$depth==d),]
      nb.gene <- length(unique(tmp$Genes))
      nb.station <- length(unique(tmp$Station))
      sum.reads <- sum(tmp$readCount)
      
      if(nb.gene>MIN.GENE & nb.station>MIN.STATION){
        Y_summary[nrow(Y_summary)+1,] <- list(c,f,d,nb.gene,nb.station, sum.reads)
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
target0$PFAMs[is.na(target0$PFAMs)] <- "Unknown"
feature0 <- stack(paste0(bluecloud.wd,"/data/features"))

# --- Selecting subset of the target0 data
# TO DO with the final dataset from Pavla, for now i test on unknown surface and ssuu
target1 <- target0[which(target0$CC_ID==CLUSTER),]
target1 <- target1[grep(pattern = paste0(DEPTH, collapse = "|"), target1$depth),]
target1 <- target1[grep(pattern = paste0(FILTER, collapse = "|"), target1$filter),]

# --- Building Y
# obs <- unique(target1$Station)
obs <- as.data.frame(unique(target1[,c("Station", "filter", "depth")]))
tar <- unique(target1$Genes)

Y0 <- matrix(0, ncol = length(tar), nrow = nrow(obs))
dimnames(Y0) <- list(paste0(obs$Station,obs$filter, obs$depth),tar)

for (i in 1:nrow(obs)){
  for (j in 1:length(tar)){
    tmp <- target1[which(target1$Genes==tar[j]  & target1$Station==obs[i,1] & target1$filter==obs[i,2] & target1$depth==obs[i,3]),]
    Y0[i,j] <- sum(tmp$readCount)
  }
}

# Y0 <- Y0+rep(seq(1,ncol(Y0)), each = nrow(Y0)) # to generate more data
# Y <- as.data.frame(Y0) # absolute abundance
Y <- as.data.frame(Y0/apply(Y0, 1, sum))

# --- Building X
obs_xy <- unique(data.frame(target1[,c("Station", "filter", "depth","Longitude", "Latitude")]))
if(length(which(duplicated(obs_xy[,c("Station","filter", "depth")])==TRUE))>0){
  obs_xy <- obs_xy[-which(duplicated(obs_xy[,c("Station","filter", "depth")])==TRUE),] # some station have two different locations... take the first one
}

X <- as.data.frame(raster::extract(feature0, obs_xy[,c("Longitude","Latitude")]))

# --- Saving data
out <- which(is.na(X[,1]) | is.na(Y[,1])) #remove point on land or no relative abundance
Y <- Y[-out,]
X <- X[-out,]

write_feather(X, path = paste0(bluecloud.wd,"/data/X.feather"))
write_feather(Y, path = paste0(bluecloud.wd,"/data/Y.feather"))

# ==================================== PART 4 ==================================
# Visual check on the target data relation to the environment

source(file = "/home/aschickele/workspace/bluecloud descriptor/00_config.R")

# --- Load data
X0 <- as.data.frame(read_feather(paste0(bluecloud.wd,"/data/X.feather")))
Y0 <- as.data.frame(read_feather(paste0(bluecloud.wd,"/data/Y.feather")))

# --- Plotting data
pal <- rep(brewer.pal(nrow(HYPERPARAMETERS), "Spectral"), each = 11)

pdf(paste0(bluecloud.wd,"/graphic/raw_data.pdf"))
for(yi in 1:ncol(Y0)){
  par(mfrow=c(3,4), bg="black", col="white", col.axis = "white", col.lab="white",col.main="white",
      mar=c(5,3,3,2))
  for(xi in 1:ncol(X0)){
    plot(as.vector(X0[,xi]), as.vector(Y0[,yi]), 
         xlab = names(X0)[xi], ylab = "relative abundance", main = paste("target n°", yi),
         pch = 16, col = pal[xi])
    grid()
  } # yi target loop
} # xi feature loop

D <- mean(apply(Y0, 1, function(x){length(which(x>0))}))
print(paste("average number of gene present by obs :", round(D, 2), "/", ncol(Y0)))

D2 <- sum(apply(Y0, 1, function(x){length(which(x>0.5))}))
print(paste("number of obs where a gene represent more than 50% of relative abundance:", D2, "/", nrow(Y0)))

while (dev.cur() > 1) dev.off()


# --- END