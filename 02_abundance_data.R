#' Here we build the relative abundance dataset with the TARA OCEAN data
#' The output are
#' 1. A relative abundance datasset : Ytargets * Nobs
#' 2. The corresponding environmental features : Xfeatures * Nobs
#' 
#' TO DO LIST:
#' - check if its needed to extract feature data from nearest neighbor pixel
#' - call the query parameters from the config file
#' - create one database per filter size and name it dynamically ?

# ================== PART 1 : CREATING THE DATABASE ============================
# To perform fast queries on the data and extract the necessary for the target
# Adapted to the final data layout

source(file = "/home/aschickele/workspace/bluecloud descriptor/00_config.R")

# --- 1. Create and open RSQLite database
unlink(paste0(bluecloud.wd, "/omic_data/omic_DB.sqlite"))
db <- dbConnect(RSQLite::SQLite(), paste0(bluecloud.wd, "/omic_data/omic_DB.sqlite"))

# --- 2. Open "clusters" and sort by n_genes
# We also create the "cluster_sort" table for future filtering of the DB
clusters <- read_feather(paste0(bluecloud.wd, "/omic_data/CC_PFAM_taxo_80cutoff.feather")) %>% 
  dplyr::select(c("Genes","PFAMs","CC_ID", "Class", "Order", "Family", "Genus"))
copy_to(db, clusters, temporary = FALSE, overwrite = TRUE)

cluster_sort <- tbl(db, "clusters") %>%
  group_by(CC_ID) %>% 
  summarise(n_genes = n_distinct(Genes), 
            unknown_rate = sum(is.na(PFAMs))/n(),
            .groups = "drop")

# --- 3. Open "reads" and filter according to n_genes
# We filter now to reduce the DB and upcoming calculation size
reads <- vroom(file = paste0(bluecloud.wd, "/omic_data/SMAGs-v1.cds.95.mg.matrix_CC_corr_80cutoff")) %>% 
  rename(Genes = 1) %>% 
  dplyr::select(contains(c("Genes", DEPTH))) %>%
  dplyr::select(contains(c("Genes", FILTER))) %>% 
  pivot_longer(!Genes, names_to = "code", values_to = "readCount") %>% 
  mutate(code = paste0("00", code),
         Station = str_sub(code, -13, -11)) %>% 
  select(-code) %>% 
  inner_join(select(clusters, "Genes"))
copy_to(db, reads, temporary = FALSE, overwrite = TRUE)

# --- 4. Open "locs" and calculate sum_reads by station
# Used for normalisation of the reads
locs <- vroom(file = paste0(bluecloud.wd, "/omic_data/SMAGs_Env_lonlat.csv")) %>% 
  dplyr::select("Station","Latitude","Longitude") %>% 
  mutate(Station = str_sub(Station, 6, 8)) %>% 
  distinct()
copy_to(db, locs, temporary = FALSE, overwrite = TRUE)

sum_station <- tbl(db, "reads") %>% 
  group_by(Station) %>% 
  summarise(sum_reads = sum(readCount, na.rm = TRUE), .groups = "drop") %>% 
  left_join(tbl(db, "locs"), by = "Station") %>% 
  select("Station","sum_reads")
copy_to(db, sum_station, temporary = FALSE, overwrite = TRUE)

# --- 5. Join "reads", "clusters" and "locs" into "data"
# To calculate supplementary filtering metrics
data <- tbl(db, "reads") %>% 
  inner_join(tbl(db, "clusters"), by = "Genes") %>% 
  left_join(tbl(db, "locs"), by = "Station")
copy_to(db, data, temporary = FALSE, overwrite = TRUE)
dbSendQuery(db, "create index by_cluster on data (CC_ID )")

# --- 6. Add "n_station", "sum_reads", "unknown_rate" to cluster_sort
tmp <- tbl(db, "data") %>% 
  group_by(CC_ID, Station) %>% 
  summarise(n_reads = sum(readCount, na.rm = TRUE), .groups = "drop_last") %>%
  filter(n_reads > 0) %>% 
  summarise(n_station = n_distinct(Station), 
            sum_reads = sum(n_reads, na.rm = TRUE))

cluster_sort <- cluster_sort %>% 
  inner_join(tmp)
copy_to(db, cluster_sort, temporary = FALSE, na.rm = TRUE)

# --- 7. Close connection
dbDisconnect(db)

# ===================== PART 2 : querying database =============================
# Extract the necessary data for the model according to the pre-defined filters
# i.e. number of genes and station per clusters

source(file = "/home/aschickele/workspace/bluecloud descriptor/00_config.R")
db <- dbConnect(RSQLite::SQLite(), paste0(bluecloud.wd, "/omic_data/omic_DB.sqlite"))

# --- 1. Filter "data" by "cluster_sort"
query <- tbl(db, "cluster_sort") %>% 
  filter(n_station >= 80 & n_genes >= 5 & n_genes <= 25)

target <- query %>% 
  select(CC_ID) %>% 
  inner_join(tbl(db, "data")) %>% 
  select(c("Genes", "CC_ID", "readCount", "Station", "Longitude", "Latitude"))

# --- 2. Reshape and normalize "data" by "sum_station"
target <- target %>% 
  group_by(CC_ID, Station, Latitude, Longitude) %>% 
  summarise(reads = sum(readCount), .groups = "drop") %>% 
  pivot_wider(names_from = CC_ID, values_from = reads) %>% 
  left_join(tbl(db, "sum_station")) %>% 
  mutate_at(.vars = vars(c(-Station, -Latitude, -Longitude, -sum_reads)), .funs = ~ . / sum_reads) %>%
  collect()

# --- 3. Building the final feature table "X"
# TO DO : Due to the coarse resolution of the features, environmental data are taken from nearest neighbor in case of NA
feature0 <- stack(paste0(bluecloud.wd,"/data/features"))
obs_xy <- target %>% 
  select(c(Longitude, Latitude))

X <- as.data.frame(raster::extract(feature0, obs_xy[,c("Longitude","Latitude")]))
write_feather(X, path = paste0(bluecloud.wd,"/data/X.feather"))

# --- 4. Building the final target table "Y"
Y <- target %>% 
  select(contains("CC"))
write_feather(Y, path = paste0(bluecloud.wd,"/data/Y.feather"))

# --- 5. Close connection
dbDisconnect(db)

# ==================================== PART 3 ==================================
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