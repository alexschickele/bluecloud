
# ===================== PART 1 : querying database =============================
# Extract the necessary data for the model according to the pre-defined filters
# i.e. number of genes and station per clusters

query_data <- function(config_file = "/home/aschickele/workspace/bluecloud descriptor/code/00a_config.R",
                       KEGG_p = "00190",
                       CLUSTER_SELEC = list(MIN_STATIONS = 80, MIN_GENES = 5, MAX_GENES = 25),
                       ENV_METRIC = c("mean","sd","med","mad","dist","bathy")){
  
  source(config_file)
  db <- dbConnect(RSQLite::SQLite(), paste0(bluecloud.wd, "/omic_data/",FILTER,"_DB.sqlite"))
  
  # --- 1. Filter "data" by "cluster_sort"
  # the "!!" are necessary for some unknown reasons
  # the KEGG_p query needs to be updated in order to avoid the copy = TRUE
  query <- dbGetQuery(db, paste('SELECT * FROM KEGG_sort WHERE KEGG_Pathway LIKE "%', KEGG_p, '%"', sep = "")) %>% 
    group_by(CC) %>% 
    summarise(KEGG_p = KEGG_p, n_KEGG = str_count(KEGG_Pathway, "ko")) %>% 
    filter(n_KEGG == 1) %>% 
    inner_join(tbl(db, "cluster_sort"), copy = TRUE) %>% 
    filter(n_station >= !!CLUSTER_SELEC$MIN_STATIONS & n_genes >= !!CLUSTER_SELEC$MIN_GENES & n_genes <= !!CLUSTER_SELEC$MAX_GENES)
  copy_to(db, query, temporary = TRUE)
  
  target <- tbl(db, "query") %>% 
    select(CC) %>% 
    inner_join(tbl(db, "data")) %>% 
    select(c("Genes", "CC", "readCount", "Station", "Longitude", "Latitude"))
  
  # --- 2. Reshape and normalize "data" by "sum_station"
  target <- target %>% 
    group_by(CC, Station, Latitude, Longitude) %>% 
    summarise(reads = sum(readCount), .groups = "drop") %>% 
    pivot_wider(names_from = CC, values_from = reads) %>% 
    left_join(tbl(db, "sum_station")) %>% 
    mutate_at(.vars = vars(c(-Station, -Latitude, -Longitude, -sum_reads)), .funs = ~ . / sum_reads)
  
  # --- 3. Building the final feature table "X"
  # TO DO : Due to the coarse resolution of the features, environmental data are taken from nearest neighbor in case of NA
  X <- target %>% 
    select(Station) %>% 
    inner_join(tbl(db, "X0")) %>% 
    select(-Station) %>% 
    select(contains(ENV_METRIC)) %>% 
    collect()
  
  write_feather(X, path = paste0(bluecloud.wd,"/data/X.feather"))
  
  # --- 4. Building the final target table "Y"
  Y <- target %>% 
    select(contains("CC")) %>% 
    collect()
  Y <- Y/max(Y, na.rm = TRUE)
  write_feather(Y, path = paste0(bluecloud.wd,"/data/Y.feather"))
  
  # --- 5. Close connection
  dbDisconnect(db)
  print(paste("Number of stations :", nrow(X)))
  print(paste("Number of environmental features :", ncol(X)))
  print(paste("Number of gene cluster targets :", ncol(Y)))
  
  return(list(X=as.data.frame(X), Y=as.data.frame(Y)))
} # end function
