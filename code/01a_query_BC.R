
# ===================== PART 1 : querying database =============================
# Extract the necessary data for the model according to the pre-defined filters
# i.e. number of genes and station per clusters

query_data <- function(bluecloud.wd = bluecloud_dir,
                       EXCLUDE = exclude,
                       KEGG_m = paste0("K",c("01595","00051","00028","00029","00814","14272","01006","14454",
                                             "14455","00024","00025","00026","01610")),
                       CLUSTER_SELEC = list(N_CLUSTERS = 50, MIN_GENES = 5, EXCLUSIVITY_R = 1),
                       ENV_METRIC = c("mean","sd","dist","bathy"),
                       relative = TRUE){
  
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
  tmp <- paste(KEGG_m, collapse = "|")
  query <- tbl(db, "kegg_sort") %>%
    dplyr::filter(grepl(pattern = tmp, kegg_ko)) %>% 
    select(CC) %>% 
    distinct() %>% 
    inner_join(tbl(db, "kegg_sort")) %>%
    dplyr::group_by(CC) %>% 
    dplyr::summarise(max_kegg = max(n_kegg, na.rm = TRUE), max_mod = max(n_mod, na.rm = TRUE)) %>% 
    inner_join(tbl(db, "cluster_sort"), copy = TRUE) %>% 
    filter(n_genes >= !!CLUSTER_SELEC$MIN_GENES & n_station >= 10)
  
  query_check <- query %>% collect()
  
  target <- query %>% 
    dplyr::select("CC") %>% 
    inner_join(tbl(db, "data")) %>% 
    dplyr::select(c("Genes", "CC", "readCount", "Station", "Longitude", "Latitude", "KEGG_ko","KEGG_Module","Description", "Phylum","Class", "Genus")) %>% 
    collect()
  
  # --- 1b. Doing exclusivity filtering in memory
  query <- target %>% 
    mutate(exclusivity = str_count(KEGG_ko, paste(c("-","NA",KEGG_m), collapse = "|"))/str_count(KEGG_ko, paste(c("-","NA","K"), collapse = "|"))) %>% 
    dplyr::group_by(CC) %>% 
    dplyr::summarise(min_exl = min(exclusivity, na.rm = TRUE)) %>% 
    filter(min_exl >= CLUSTER_SELEC$EXCLUSIVITY_R) %>% 
    inner_join(query_check)
  
  tmp <- query %>% select("CC")
  
  target <- target %>% 
    inner_join(tmp)
  
  # --- END bluecloud specificity, back to normal logic  
  
  target <- target %>% mutate(MAG = gsub("(.*_){2}(\\d+)_.+", "\\2", Genes)) %>% # getting MAGs
    dplyr::filter(!grepl(EXCLUDE, Class)) %>% # Exclude non sense classes
    group_by(CC, Station) %>%
    mutate(Unknown_rate = sum(is.na(KEGG_ko))*100/n()) %>% collect() # defining unknown rate here instead of 00c_build_omic_data.R

  # --- 2. Get cluster functional description
  CC_desc <- target %>% 
    inner_join(tbl(db, "cluster_sort"), copy = TRUE) %>% 
    group_by(CC) %>% 
    summarise(unknown_rate = paste(unique(Unknown_rate), collapse = ", "),
              kegg_ko = paste(unique(KEGG_ko), collapse = ", "),
              kegg_module = paste(unique(KEGG_Module), collapse = ", "),
              desc = paste(unique(Description), collapse = ", "),
              phylum = paste(unique(Phylum), collapse = ", "),
              class = paste(unique(Class), collapse = ", "),
              genus = paste(unique(Genus), collapse = ", "),
              mag = paste(unique(MAG), collapse = ", ")) %>% 
    inner_join(query_check)

  write_feather(CC_desc, path = paste0(bluecloud.wd,"/data/CC_desc.feather"))

  # --- 3. Reshape and normalize "data" by "sum_station"
  target <- target %>% 
    group_by(CC, Station, Latitude, Longitude) %>% 
    summarise(reads = sum(readCount), .groups = "drop") %>% 
    pivot_wider(names_from = CC, values_from = reads) %>% 
    left_join(tbl(db, "sum_station"), copy = TRUE) %>% 
    mutate_at(.vars = vars(c(-Station, -Latitude, -Longitude, -sum_reads)), .funs = ~ . / sum_reads)
  
  # --- 3b. Security to put target back in the right order (needed for BlueCloud)
  id <- order(colnames(target[,4:ncol(target)]))
  target <- target[, c(1,2,3, id+3)]
  
  # --- 4. Selecting the clusters most representative of the variability
  e <- target %>% 
    dplyr::select(c(-Station, -Latitude, -Longitude, -sum_reads)) %>% 
    escouf()
  
  target0 <- target %>% 
    dplyr::select(c(Station, Latitude, Longitude, e$vr[1:length(e$vr)]+3)) # for Y0
  target <- target %>% 
    dplyr::select(c(Station, Latitude, Longitude, e$vr[1:min(length(e$vr),CLUSTER_SELEC$N_CLUSTERS)]+3))

  # --- 5. Building the final feature table "X"
  X <- target %>% 
    dplyr::select(Station) %>% 
    inner_join(tbl(db, "X0"), copy = TRUE) %>% 
    arrange(Station) %>% # security too re-arrange station names alphabetically
    dplyr::select(contains(c(ENV_METRIC)))
  
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
  Y <- Y0[1:min(ncol(Y0), CLUSTER_SELEC$N_CLUSTERS)]
  Y <- apply(as.matrix(Y), 1, function(x){if(sum(x)>0){x = x/sum(x, na.rm = TRUE)} else {x = x}}) %>%
        aperm(c(2,1)) %>%
        as.data.frame()
  
  write_feather(Y, path = paste0(bluecloud.wd,"/data/Y.feather"))
  
  # --- 8. Extract the vector of station names
  ID <- target %>% 
    dplyr::select(Station)
  ID <- sort(ID$Station)
  write_feather(data.frame(Station = ID), path = paste0(bluecloud.wd,"/data/Station_ID.feather"))
  
  # --- 9. Correspondence analysis
  # Assign points to the nearest escoufier selected point (correspondence analysis)
  res_ca <- CA(Y0, graph = FALSE, ncp = 50)
  ndim_ca <- which(res_ca$eig[,3]>80)[1]
  dist_ca <- dist(res_ca$col$coord[,1:ndim_ca]) %>% as.matrix()
  dist_ca <- dist_ca[,1:ncol(Y)]
  
  # Quantify dominance of CC
  tmp <- apply(Y0,2,sum)
  
  # Concatenate in dataframe
  nn_ca <- data.frame(CC = dimnames(dist_ca)[[1]],
                      pos_CC = 1:nrow(dist_ca),
                      nn_CC = apply(dist_ca, 1, function(x){x = names(which(x == min(x)))}),
                      pos_nn_CC = apply(dist_ca, 1, function(x){x = which(x == min(x))}),
                      sum_CC = tmp,
                      scale_CC = tmp/tmp[apply(dist_ca, 1, function(x){x = which(x == min(x))})])
  
  # --- 10. Close connection
  print(paste("Number of stations :", nrow(X)))
  print(paste("Number of environmental features :", ncol(X)))
  print(paste("Number of gene/cluster targets :", ncol(Y)))
  
  return(list(X=as.data.frame(X), Y=as.data.frame(Y), Y0=as.data.frame(Y0), CC_desc = as.data.frame(CC_desc), e = e, nn_ca = nn_ca))
  
  # --- Close connection
  dbDisconnect(db)
} # end function
