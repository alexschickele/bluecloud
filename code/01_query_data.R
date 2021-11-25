
# ===================== PART 1 : querying database =============================
# Extract the necessary data for the model according to the pre-defined filters
# i.e. number of genes and station per clusters

query_data <- function(bluecloud.wd = bluecloud_dir,
                       KEGG_p = "00190",
                       CLUSTER_SELEC = list(MIN_STATIONS = 80, MIN_GENES = 5, MAX_GENES = 25),
                       ENV_METRIC = c("mean","sd","dist","bathy"),
                       relative = TRUE){
  
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
  query <- tbl(db, "kegg_sort") %>% 
    filter(str_detect(kegg_pathway, KEGG_p)) %>% 
    dplyr::group_by(CC) %>% 
    dplyr::summarise(max_kegg = max(n_kegg, na.rm = TRUE)) %>% 
    filter(max_kegg == 1) %>% 
    inner_join(tbl(db, "cluster_sort"), copy = TRUE) %>% 
    filter(n_station >= !!CLUSTER_SELEC$MIN_STATIONS & n_genes >= !!CLUSTER_SELEC$MIN_GENES & n_genes <= !!CLUSTER_SELEC$MAX_GENES)
  
  check_query <- query %>% collect()

  if(nrow(check_query)!=0){
    target <- query %>% 
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
      select(contains(c(ENV_METRIC))) %>% 
      collect()
  
    write_feather(X, path = paste0(bluecloud.wd,"/data/X.feather"))
    
    # --- 5. Building the final target table "Y"
    Y <- target %>% 
      select(contains(c("Station","CC"))) %>%
      collect() %>% 
      arrange(Station) %>% 
      select(-Station)
    Y <- Y/max(Y, na.rm = TRUE)
    if(relative == TRUE){
      Y <- apply(as.matrix(Y), 1, function(x){if(sum(x)>0){x = x/sum(x, na.rm = TRUE)} else {x = x}}) %>%
        aperm(c(2,1)) %>%
        as.data.frame()
    }
    write_feather(Y, path = paste0(bluecloud.wd,"/data/Y.feather"))
    
    # --- 6. Extract the vector of station names
    ID <- target %>% 
      select(Station) %>% 
      collect()
    ID <- sort(ID$Station)
    write_feather(data.frame(Station = ID), path = paste0(bluecloud.wd,"/data/Station_ID.feather"))
    
    # --- 7. Close connection
    print(paste("Number of stations :", nrow(X)))
    print(paste("Number of environmental features :", ncol(X)))
    print(paste("Number of gene cluster targets :", ncol(Y)))
    
    return(list(X=as.data.frame(X), Y=as.data.frame(Y), CC_desc = as.data.frame(CC_desc)))
  } else {message("STOP: the query has no result. Try to increase the gene range or minimum number of stations")} # if continue
  
  # --- Close connection
  dbDisconnect(db)
} # end function
