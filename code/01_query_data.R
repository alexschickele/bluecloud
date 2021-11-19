
# ===================== PART 1 : querying database =============================
# Extract the necessary data for the model according to the pre-defined filters
# i.e. number of genes and station per clusters

query_data <- function(bluecloud.wd = "/home/jovyan/bluecloud",
                       KEGG_p = "00190",
                       CLUSTER_SELEC = list(MIN_STATIONS = 50, MIN_GENES = 5, MAX_GENES = 25),
                       ENV_METRIC = c("mean","sd","med","mad","dist","bathy")){
  
  # --- For local database
  # db <- dbConnect(RSQLite::SQLite(), paste0(bluecloud.wd, "/omic_data/",FILTER,"_DB.sqlite"))
  
  # --- For bluecloud
  db <- dbConnect(
    drv=PostgreSQL(),
    host="postgresql-srv.d4science.org",
    dbname="bluecloud_demo2",
    user="bluecloud_demo2_u",
    password="6a26c54a05ec5dede958a370ca744a",
    port=5432
  )
  
  # --- 1. Filter "data" by "cluster_sort"
  # the "!!" are necessary for some unknown reasons
  # the KEGG_p query needs to be updated in order to avoid the copy = TRUE
  query <- dbGetQuery(db, paste0("SELECT * FROM public.kegg_sort WHERE kegg_pathway LIKE '%", KEGG_p,"%'")) %>% 
    dplyr::group_by(CC) %>% 
    dplyr::summarise(kegg_p = KEGG_p, n_kegg = max(str_count(kegg_pathway, "ko"))) %>% 
    filter(n_kegg == 1) %>% 
    inner_join(tbl(db, "cluster_sort"), copy = TRUE) %>% 
    filter(n_station >= !!CLUSTER_SELEC$MIN_STATIONS & n_genes >= !!CLUSTER_SELEC$MIN_GENES & n_genes <= !!CLUSTER_SELEC$MAX_GENES)
  copy_to(db, query, temporary = FALSE)
  
  if(nrow(query)!=0){
    target <- tbl(db, "query") %>% 
      select(CC) %>% 
      inner_join(tbl(db, "data")) %>% 
      select(c("Genes", "CC", "readCount", "Station", "Longitude", "Latitude", "Description"))
    
    # --- 2. Get cluster functional description
    CC_desc <- target %>% 
      collect() %>% 
      group_by(CC, Description) %>% 
      summarise() %>% 
      summarise(desc = paste(Description, collapse = ", "))
    write_feather(CC_desc, path = paste0(bluecloud.wd,"/data/CC_desc.feather"))
    
    # --- 3. Reshape and normalize "data" by "sum_station"
    target <- target %>% 
      group_by(CC, Station, Latitude, Longitude) %>% 
      summarise(reads = sum(readCount), .groups = "drop") %>% 
      pivot_wider(names_from = CC, values_from = reads) %>% 
      left_join(tbl(db, "sum_station")) %>% 
      mutate_at(.vars = vars(c(-Station, -Latitude, -Longitude, -sum_reads)), .funs = ~ . / sum_reads)
    
    # --- 4. Building the final feature table "X"
    X <- target %>% 
      select(Station) %>% 
      inner_join(tbl(db, "X0")) %>% 
      select(-Station) %>% 
      select(contains(ENV_METRIC)) %>% 
      collect()
  
    write_feather(X, path = paste0(bluecloud.wd,"/data/X.feather"))
    
    # --- 5. Building the final target table "Y"
    Y <- target %>% 
      select(contains("CC")) %>% 
      collect()
    Y <- Y/max(Y, na.rm = TRUE)
    write_feather(Y, path = paste0(bluecloud.wd,"/data/Y.feather"))
    
    # --- 6. Close connection
    
    print(paste("Number of stations :", nrow(X)))
    print(paste("Number of environmental features :", ncol(X)))
    print(paste("Number of gene cluster targets :", ncol(Y)))
    
    return(list(X=as.data.frame(X), Y=as.data.frame(Y), CC_desc = as.data.frame(CC_desc)))
  } else {message("STOP: the query has no result. Try to increase the gene range or minimum number of stations")} # if continue
  
  # --- Close connection
  dbRemoveTable(db, "query")
  dbDisconnect(db)
} # end function
