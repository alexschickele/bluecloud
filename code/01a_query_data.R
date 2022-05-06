
# ===================== PART 1 : querying database =============================
# Extract the necessary data for the model according to the pre-defined filters
# i.e. number of genes and station per clusters

query_data <- function(bluecloud.wd = bluecloud_dir,
                       CC_id = NULL,
                       KEGG_m = 165:172,
                       CLUSTER_SELEC = list(N_CLUSTERS = 30, MIN_GENES = 5, EXCLUSIVITY_R = 1),
                       ENV_METRIC = c("mean","sd","dist","bathy"),
                       relative = TRUE){
  
  # --- For local database
  db <- dbConnect(RSQLite::SQLite(), paste0(bluecloud.wd, "/omic_data/",FILTER,"_DB.sqlite"))
  
  # --- For bluecloud
  # db <- dbConnect(drv=PostgreSQL(), host="postgresql-srv.d4science.org", dbname="bluecloud_demo2",
  #   user="bluecloud_demo2_u", password="6a26c54a05ec5dede958a370ca744a", port=5432)

  # --- 1. Filter "data" by "cluster_sort"
  if(is.null(CC_id)){
    # query <- dbGetQuery(db, paste0("SELECT CC FROM kegg_sort WHERE kegg_module LIKE '%",
    #                                paste(KEGG_m, collapse = "%' OR kegg_module LIKE '%"), "%'")) %>%
    query <- dbGetQuery(db, paste0("SELECT CC FROM kegg_sort WHERE kegg_ko LIKE '%",
                                  paste(KEGG_m, collapse = "%' OR kegg_ko LIKE '%"), "%'")) %>%
      unique() %>% 
      inner_join(tbl(db, "kegg_sort"), copy = TRUE) %>%
      # query <- tbl(db, "kegg_sort") %>% 
      #   filter(str_detect(kegg_pathway, KEGG_p)) %>% 
      # mutate(exclusivity = str_count(kegg_module, paste(c("-","NA",KEGG_m), collapse = "|"))/n_mod) %>%
      mutate(exclusivity = str_count(kegg_ko, paste(c("-","NA",KEGG_m), collapse = "|"))/n_ko) %>%
      dplyr::group_by(CC) %>% 
      dplyr::summarise(max_kegg = max(n_kegg, na.rm = TRUE), max_mod = max(n_mod, na.rm = TRUE), min_exl = min(exclusivity, na.rm = TRUE)) %>% 
      filter(min_exl >= CLUSTER_SELEC$EXCLUSIVITY_R) %>%
      inner_join(tbl(db, "cluster_sort"), copy = TRUE) %>% 
      filter(n_genes >= !!CLUSTER_SELEC$MIN_GENES & n_station >= 10)
    
    query_check <- query
    copy_to(db, query, overwrite = TRUE)
    query <- tbl(db, "query")

    target <- query %>% 
      dplyr::select("CC") %>% 
      inner_join(tbl(db, "data")) %>% 
      dplyr::select(c("Genes", "CC", "readCount", "Station", "Longitude", "Latitude", "KEGG_ko","KEGG_Module","Description", "Phylum","Class", "Genus")) %>% 
      collect()
    target <- target %>% mutate(MAG = gsub("(.*_){2}(\\d+)_.+", "\\2", Genes))
    copy_to(db, target, overwrite = TRUE)
    target <- tbl(db, "target")
    
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
                phylum = paste(unique(Phylum), collapse = ", "),
                class = paste(unique(Class), collapse = ", "),
                genus = paste(unique(Genus), collapse = ", "),
                mag = paste(unique(MAG), collapse = ", ")) %>% 
      inner_join(query_check)
  } else {
    CC_desc <- target %>% 
      select("Genes","KEGG_ko","KEGG_Module", "Description", "Class", "Genus") %>% 
      distinct() %>% 
      collect()
    names(CC_desc) <- c("Genes","kegg_ko","kegg_module","desc","class","genus")
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
  
  # --- 4. Selecting the clusters most representative of the variability
    e <- target %>% 
      dplyr::select(c(-Station, -Latitude, -Longitude, -sum_reads)) %>% 
      collect() %>% 
      escouf()
    
    target0 <- target %>% 
      dplyr::select(c(Station, Latitude, Longitude, e$vr[1:length(e$vr)]+3)) # for Y0
    target <- target %>% 
      dplyr::select(c(Station, Latitude, Longitude, e$vr[1:min(length(e$vr),CLUSTER_SELEC$N_CLUSTERS)]+3))

  # --- 5. Building the final feature table "X"
  X <- target %>% 
    dplyr::select(Station) %>% 
    inner_join(tbl(db, "X0")) %>% 
    dplyr::select(contains(c(ENV_METRIC))) %>% 
    collect()
  
  write_feather(X, path = paste0(bluecloud.wd,"/data/X.feather"))
  
  # --- 6. Building the final target table "Y0" for interpretation
  Y0 <- target0 %>% 
    dplyr::select(-Latitude, -Longitude) %>%
    collect() %>% 
    arrange(Station) %>% 
    dplyr::select(-Station)
  Y0 <- Y0/max(Y0, na.rm = TRUE)
  if(relative == TRUE){
    Y0 <- apply(as.matrix(Y0), 1, function(x){if(sum(x)>0){x = x/sum(x, na.rm = TRUE)} else {x = x}}) %>%
      aperm(c(2,1)) %>%
      as.data.frame()
  }
  
  # --- 7. Building the final target table "Y"
  Y <- Y0[1:50]
  Y <- apply(as.matrix(Y), 1, function(x){if(sum(x)>0){x = x/sum(x, na.rm = TRUE)} else {x = x}}) %>%
        aperm(c(2,1)) %>%
        as.data.frame()
  
  write_feather(Y, path = paste0(bluecloud.wd,"/data/Y.feather"))
  
  # --- 8. Extract the vector of station names
  ID <- target %>% 
    dplyr::select(Station) %>% 
    collect()
  ID <- sort(ID$Station)
  write_feather(data.frame(Station = ID), path = paste0(bluecloud.wd,"/data/Station_ID.feather"))
  
  # --- 9. Close connection
  print(paste("Number of stations :", nrow(X)))
  print(paste("Number of environmental features :", ncol(X)))
  print(paste("Number of gene/cluster targets :", ncol(Y)))
  
  return(list(X=as.data.frame(X), Y=as.data.frame(Y), Y0=as.data.frame(Y0), CC_desc = as.data.frame(CC_desc), e = e))
  
  # --- Close connection
  dbDisconnect(db)
} # end function
