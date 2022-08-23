#' @concept building a RSQLite database from which we extract feature and
#' target data
#' 
#' @source TARA Ocean MetaG reads
#' @source TARA Ocean Cluster of sequence similarity table, anotated by KEGG and Pathways
#' @source TARA Ocean station locations
#' @source Environmental data calculated from script 01_
#' 
#' @param bluecloud.wd path to the bluecloud descriptor file
#' @param FILTER class size from which the database is build
#' @param CLUSTER_SELEC list of filters to select clusters of appropriate size :
#' i.e. minimum number of stations, minimum number of genes, maximum number of genes
#' 
#' @details the extraction of environmental variable corresponding to each station
#' location is done from nearest neighbor in case of NA, due to the coarse 1° resolution
#' 
#' @return an RSQLite database containing all necessary data for the models and queries
#' @return X : n_obs x n_env_variable .feather of features
#' @return Y : n_obs x n_clusters .feather of targets

# ================== PART 1 : CREATING THE DATABASE ============================
# To perform fast queries on the data and extract the necessary for the target
# Adapted to the final data layout

source(file = "./code/00a_config.R")

# --- 1. Create and open RSQLite database on Marie
# unlink(paste0(bluecloud.wd, "/omic_data/",FILTER,"_DB.sqlite"))
db <- dbConnect(RSQLite::SQLite(), paste0(bluecloud.wd, "/omic_data/",FILTER,"_DB.sqlite"))

# --- 1B. Create and open RPostgreSQL database on BlueCloud
# db <- dbConnect(
#   drv=PostgreSQL(),
#   host="postgresql-srv.d4science.org",
#   dbname="bluecloud_demo2",
#   user="bluecloud_demo2_writer",
#   password="toto",
#   port=5432
# )

# --- 2. Open "clusters" and sort by n_genes
taxo <- read_feather(paste0(bluecloud.wd, "/omic_data/CC_PFAM_taxo_80cutoff.feather")) %>% 
  dplyr::select(c("Genes","Phylum","Class","Order","Family","Genus"))
clusters <- read_feather(paste0(bluecloud.wd, "/omic_data/CC_80_withallannot.feather")) %>% 
  dplyr::select(c("Genes","CC", "COG_category","GOs", "KEGG_ko", "KEGG_Pathway", "KEGG_Module", "PFAMs", "Description")) %>% 
  left_join(taxo)
copy_to(db, clusters, temporary = FALSE, overwrite = TRUE)

# --- 3. Open "reads" and filter according to n_genes
# We filter now to reduce the DB and upcoming calculation size
reads <- vroom(file = paste0(bluecloud.wd, "/omic_data/SMAGs-v1.cds.95.mg.matrix_CC_corr_80cutoff")) %>% 
  rename(Genes = 1) %>% 
  dplyr::select(contains(c("Genes", DEPTH))) %>%
  dplyr::select(contains(c("Genes", FILTER))) %>% 
  pivot_longer(!Genes, names_to = "code", values_to = "readCount") %>% 
  mutate(code = paste0("00", code),
         Station = str_sub(code, -13, -11)) %>% 
  dplyr::select(-code) %>% 
  inner_join(dplyr::select(clusters, "Genes"))

gene_length <- vroom(file = paste0(bluecloud.wd, "/omic_data/EukMAGS_SMAGs_nucl_concat_length"), col_names = c("Genes", "kb")) %>% 
  inner_join(reads) %>% 
  mutate(readperkb = readCount / kb)

reads <- gene_length %>% 
  dplyr::select(contains(c("Genes","readperkb","Station"))) %>% 
  rename(readCount = readperkb)

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
  dplyr::select("Station","sum_reads")
copy_to(db, sum_station, temporary = FALSE, overwrite = TRUE)

# --- 5. Join "reads", "clusters" and "locs" into "data"
# To calculate supplementary filtering metrics
data <- tbl(db, "reads") %>% 
  inner_join(tbl(db, "clusters"), by = "Genes") %>% 
  left_join(tbl(db, "locs"), by = "Station")
copy_to(db, data, temporary = FALSE, overwrite = TRUE)
dbSendQuery(db, "create index by_cluster on data (CC)")
# dbSendQuery(db, 'CREATE INDEX by_cluster ON public.data ("CC")')

# --- 6. Create cluster_sort and add "n_station", "sum_reads", "unknown_rate" to cluster_sort
# /!\  created from "data" instead of "clusters" because data does not account for 0 reads genes
tmp <- tbl(db, "data") %>% 
  group_by(CC, Station) %>% 
  summarise(n_reads = sum(readCount, na.rm = TRUE), .groups = "drop_last") %>%
  filter(n_reads > 0) %>% 
  summarise(n_station = n_distinct(Station), 
            sum_reads = sum(n_reads, na.rm = TRUE))

cluster_sort <- tbl(db, "data") %>% 
  inner_join(tmp)  %>% 
  dplyr::group_by(CC, n_station, sum_reads) %>% 
  dplyr::summarise(n_genes = n_distinct(Genes), 
                   unknown_rate = sum(is.na(PFAMs))*100/n(), #NOT WORKING WITH PostgreSQL !!
                   .groups = "drop")
copy_to(db, cluster_sort, temporary = FALSE, na.rm = TRUE)
dbSendQuery(db, "create index by_cluster_sort on cluster_sort (CC)")
# dbSendQuery(db, 'CREATE INDEX by_cluster_sort ON public.cluster_sort ("CC")')

# --- 7. Add correspondence Cluster - KEGG_Pathway
kegg_sort <- tbl(db, "data") %>% 
  dplyr::select("Genes", "CC", "KEGG_Pathway", "KEGG_Module", "KEGG_ko") %>% 
  rename(kegg_pathway = KEGG_Pathway, kegg_module = KEGG_Module, kegg_ko = KEGG_ko) %>% #Otherwise I cannot do the query() in the next script !! ... on BlueCLoud postgresql
  collect() %>% 
  mutate(n_kegg = str_count(kegg_pathway, "ko"), n_mod = str_count(kegg_module, "M"), n_ko = str_count(kegg_ko, "K")) #Calculated out of SQL here to have read-only query later
copy_to(db, kegg_sort, temporary = FALSE, na.rm = FALSE)

# --- 8. Pre-calculate feature table from nearest non-NA values
# To only query the station later and not extract during queries section
X0 <- tbl(db, "locs") %>% collect()
xy <- X0 %>% dplyr::select(x = Longitude, y = Latitude)

features <- stack(paste0(bluecloud.wd,"/data/features")) %>% 
  readAll()

sample_raster_NA <- function(r, xy){
  apply(X = xy, MARGIN = 1, 
        FUN = function(xy) r@data@values[which.min(replace(distanceFromPoints(r, xy), is.na(r), NA))])
}
sampled <- mclapply(features@layers, function(a_layer) sample_raster_NA(a_layer, xy), mc.cores = 10) %>% 
  as.data.frame() %>% 
  mutate(Station = X0$Station, .before = 1)
X0 <- sampled
names(X0) <- c("Station", names(features))
copy_to(db, X0, temporary = FALSE, na.rm = TRUE)

# --- 9. Close connection
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