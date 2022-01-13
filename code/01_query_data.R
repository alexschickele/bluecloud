
# ===================== PART 1 : querying database =============================
# Extract the necessary data for the model according to the pre-defined filters
# i.e. number of genes and station per clusters

query_data <- function(bluecloud.wd = bluecloud_dir,
                       CC_id = NULL,
                       KEGG_m = 165:172,
                       CLUSTER_SELEC = list(N_CLUSTERS = 25, MIN_GENES = 2, MAX_GENES = 25),
                       ENV_METRIC = c("mean","sd","dist","bathy"),
                       relative = TRUE){
  
  # --- For local database
  db <- dbConnect(RSQLite::SQLite(), paste0(bluecloud.wd, "/omic_data/",FILTER,"_DB.sqlite"))
  
  # --- For bluecloud
  # db <- dbConnect(drv=PostgreSQL(), host="postgresql-srv.d4science.org", dbname="bluecloud_demo2",
  #   user="bluecloud_demo2_u", password="6a26c54a05ec5dede958a370ca744a", port=5432)

  # --- 1. Filter "data" by "cluster_sort"
  if(is.null(CC_id)){
    query <- dbGetQuery(db, paste0("SELECT * FROM kegg_sort WHERE kegg_module LIKE '%", 
                                   paste(KEGG_m, collapse = "%' OR kegg_module LIKE '%"), "%'")) %>% 
      # query <- tbl(db, "kegg_sort") %>% 
      #   filter(str_detect(kegg_pathway, KEGG_p)) %>% 
      mutate(exclusivity = str_count(kegg_module, paste(c("-","NA",KEGG_m), collapse = "|"))/n_mod) %>% 
      dplyr::group_by(CC) %>% 
      dplyr::summarise(max_kegg = max(n_kegg, na.rm = TRUE), max_mod = max(n_mod, na.rm = TRUE), min_exl = min(exclusivity)) %>% 
      filter(max_kegg > 0 & max_mod > 0 & min_exl == 1) %>% 
      inner_join(tbl(db, "cluster_sort"), copy = TRUE) %>% 
      arrange(desc(n_station)) %>% 
      filter(n_genes >= !!CLUSTER_SELEC$MIN_GENES & n_genes <= !!CLUSTER_SELEC$MAX_GENES) %>% 
      slice(1:CLUSTER_SELEC$N_CLUSTERS)
    
    copy_to(db, query, overwrite = TRUE)
    query <- tbl(db, "query")

    target <- query %>% 
      select(CC) %>% 
      inner_join(tbl(db, "data")) %>% 
      select(c("Genes", "CC", "readCount", "Station", "Longitude", "Latitude", "KEGG_ko","KEGG_Module","Description", "Class", "Genus"))
  } else {
    target <- tbl(db, "data") %>% 
      filter(CC == CC_id)
  }
  
  # --- 2. Get cluster functional description
  if(is.null(CC_id)){
    CC_desc <- target %>% 
      inner_join(tbl(db, "cluster_sort")) %>% 
      collect() %>% 
      group_by(CC, unknown_rate) %>% 
      summarise(kegg_ko = paste(unique(KEGG_ko), collapse = ", "),
                kegg_module = paste(unique(KEGG_Module), collapse = ", "),
                desc = paste(unique(Description), collapse = ", "),
                class = paste(unique(Class), collapse = ", "),
                genus = paste(unique(Genus), collapse = ", "))
  } else {
    CC_desc <- target %>% 
      select("Genes","KEGG_ko","KEGG_Module", "Description", "Class", "Genus") %>% 
      distinct() %>% 
      collect()
  }
  write_feather(CC_desc, path = paste0(bluecloud.wd,"/data/CC_desc.feather"))

  # --- 3. Reshape and normalize "data" by "sum_station"
  if(is.null(CC_id)){
    target <- target %>% 
      group_by(CC, Station, Latitude, Longitude) %>% 
      summarise(reads = sum(readCount), .groups = "drop") %>% 
      pivot_wider(names_from = CC, values_from = reads) %>% 
      left_join(tbl(db, "sum_station")) %>% 
      mutate_at(.vars = vars(c(-Station, -Latitude, -Longitude, -sum_reads)), .funs = ~ . / sum_reads)
  } else {
    target <- target %>% 
      group_by(Genes, Station, Latitude, Longitude) %>% 
      summarise(reads = sum(readCount), .groups = "drop") %>% 
      pivot_wider(names_from = Genes, values_from = reads) %>% 
      left_join(tbl(db, "sum_station")) %>% 
      mutate_at(.vars = vars(c(-Station, -Latitude, -Longitude, -sum_reads)), .funs = ~ . / sum_reads)
  }

  # --- 4. Building the final feature table "X"
  X <- target %>% 
    select(Station) %>% 
    inner_join(tbl(db, "X0")) %>% 
    select(contains(c(ENV_METRIC))) %>% 
    collect()
  
  write_feather(X, path = paste0(bluecloud.wd,"/data/X.feather"))
  
  # --- 5. Building the final target table "Y"
  Y <- target %>% 
    select(-Latitude, -Longitude, -sum_reads) %>%
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
  
  # --- 7. Plot CC vs kegg_modules
  if(is.null(CC_id)){
      CC_module <- matrix(NA, ncol = length(KEGG_m), nrow = CLUSTER_SELEC$N_CLUSTERS, 
                          dimnames = list(CC_desc$CC, KEGG_m))
      for(j in 1:nrow(CC_module)){
        for(k in 1:ncol(CC_module)){
          if(str_detect(CC_desc$kegg_module[j], as.character(KEGG_m[k])) == TRUE){CC_module[j,k] <- 1}
        } #k module
      } # j CC
      CC_module <- CC_module[do.call(order, as.data.frame(CC_module)),]
  } else {CC_module <- NULL}


  # --- 8. Close connection
  print(paste("Number of stations :", nrow(X)))
  print(paste("Number of environmental features :", ncol(X)))
  print(paste("Number of gene/cluster targets :", ncol(Y)))
  
  return(list(X=as.data.frame(X), Y=as.data.frame(Y), CC_desc = as.data.frame(CC_desc), CC_module = CC_module))
  
  # --- Close connection
  dbDisconnect(db)
} # end function
