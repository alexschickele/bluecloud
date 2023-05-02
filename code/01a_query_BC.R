
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
  
  # --- 1. Enter database credentials for bluecloud
  db <- dbConnect(
    drv=PostgreSQL(),
    host="postgresql-srv.d4science.org",
    dbname="bluecloud_demo2",
    user="bluecloud_demo2_u",
    password="6a26c54a05ec5dede958a370ca744a",
    port=5432
  )

  # --- 2. Initial database query
  # 2.1. List KEGG modules of interest
  tmp <- paste(KEGG_m, collapse = "|")
  # 2.2. Extract related Connected Components and filter by minimum genes and number of stations
  query <- tbl(db, "kegg_sort") %>%
    dplyr::filter(grepl(pattern = tmp, kegg_ko)) %>% 
    select(CC) %>% 
    distinct() %>% 
    inner_join(tbl(db, "kegg_sort")) %>%
    dplyr::group_by(CC) %>% 
    dplyr::summarise(max_kegg = max(n_kegg, na.rm = TRUE), max_mod = max(n_mod, na.rm = TRUE)) %>% 
    inner_join(tbl(db, "cluster_sort"), copy = TRUE) %>% 
    filter(n_genes >= !!CLUSTER_SELEC$MIN_GENES & n_station >= 10)
  
  # 2.3. Save initial query in another object for later
  query_check <- query %>% collect()
  
  # 2.4. Extract reads per station and corresponding functional annotation
  target <- query %>% 
    dplyr::select("CC") %>% 
    inner_join(tbl(db, "data")) %>% 
    dplyr::select(c("Genes", "CC", "readCount", "Station", "Filter","Longitude", "Latitude", "KEGG_ko","KEGG_Module","Description", "Phylum","Class", "Genus")) %>% 
    collect()
  
  # --- 3. Filter the connected component by the exclusivity criteria
  # 3.1. Calculate the exclusivity and filter the previous table
  query <- target %>% 
    mutate(exclusivity = str_count(KEGG_ko, paste(c("-","NA",KEGG_m), collapse = "|"))/str_count(KEGG_ko, paste(c("-","NA","K"), collapse = "|"))) %>% 
    dplyr::group_by(CC) %>% 
    dplyr::summarise(min_exl = min(exclusivity, na.rm = TRUE)) %>% 
    filter(min_exl >= CLUSTER_SELEC$EXCLUSIVITY_R) %>% 
    inner_join(query_check)
  
  # 3.2. Extract the list of corresponding Connected Component
  tmp <- query %>% select("CC")
  
  # 3.3. Filter the station - level table accordingly
  target <- target %>% 
    inner_join(tmp)
  
  # --- 4. Supplementary columns
  # Extract the MAG identifier and the unknown rate per Connected Component in dedicated columns
  target <- target %>% mutate(MAG = gsub("(.*_){2}(\\d+)_.+", "\\2", Genes)) %>% # getting MAGs
    dplyr::filter(!grepl(EXCLUDE, Class)) %>% # Exclude non sense classes
    group_by(CC, Station) %>%
    mutate(Unknown_rate = sum(is.na(KEGG_ko))*100/n()) %>% collect() # defining unknown rate here instead of 00c_build_omic_data.R

  # --- 5. Get cluster functional description
  # Concatenate the functional annotations at the Connected Component level
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

  # --- 6. Normalize the read counts by the sum of reads per stations
  # Necessary to alleviate any variation in the sequencing depth between stations
  target <- target %>% 
    group_by(CC, Station, Filter, Latitude, Longitude) %>% 
    summarise(reads = sum(readCount), .groups = "drop") %>% 
    pivot_wider(names_from = CC, values_from = reads) %>% 
    left_join(tbl(db, "sum_station"), copy = TRUE) %>% 
    mutate_at(.vars = vars(c(-Station, -Filter, -Latitude, -Longitude, -sum_reads)), .funs = ~ . / sum_reads)
  
  # --- 7. Average GGMM and GGZZ
  # Average the values of the two filters in their common station to avoid creating a sampling bias
  # The composition in terms of reads per station and connected components between filters has been tested
  target <- target %>% 
    dplyr::select(-Filter, -sum_reads) %>% 
    group_by(Station, Latitude, Longitude) %>% 
    mutate(across(everything(), mean)) %>% 
    ungroup() %>% 
    distinct() %>% 
    collect() %>% 
    mutate(Station = str_pad(Station, 3, "0", side = "left")) %>% 
    dplyr::filter(!is.na(Latitude) & !is.na(Longitude))
  
  # --- 8. Security to put target back in the right order (needed for BlueCloud)
  id <- order(colnames(target[,4:ncol(target)]))
  target <- target[, c(1,2,3, id+3)]
  
  # --- 9. Selecting the clusters most representative of the total variance
  # 9.1. Performing an escoufier dimensional reduction
  e <- target %>% 
    dplyr::select(c(-Station, -Latitude, -Longitude)) %>% 
    escouf()
  
  # 9.2. Updating the target table accordingly for later
  target0 <- target %>% 
    dplyr::select(c(Station, Latitude, Longitude, e$vr[1:length(e$vr)]+3)) # for Y0
  target <- target %>% 
    dplyr::select(c(Station, Latitude, Longitude, e$vr[1:min(length(e$vr),CLUSTER_SELEC$N_CLUSTERS)]+3))

  # --- 10. Building the final feature table "X"
  # Extracting the environmental values at the corresponding station coordinates
  X <- target %>% 
    dplyr::select(Station) %>% 
    inner_join(tbl(db, "X0"), copy = TRUE) %>% 
    arrange(Station) %>% # security too re-arrange station names alphabetically
    dplyr::select(contains(c(ENV_METRIC)))
  
  write_feather(X, path = paste0(bluecloud.wd,"/data/X.feather"))
  
  # --- 11. Building the final target table "Y0"
  # This table contains all targets, selected or not by escoufier, for later use when reconstructing the projections
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
  
  # --- 12. Building the final target table "Y"
  # Only containing the escoufier selected connected components, for direct model fitting purpose
  Y <- Y0[1:min(ncol(Y0), CLUSTER_SELEC$N_CLUSTERS)]
  Y <- apply(as.matrix(Y), 1, function(x){if(sum(x)>0){x = x/sum(x, na.rm = TRUE)} else {x = x}}) %>%
        aperm(c(2,1)) %>%
        as.data.frame()
  
  write_feather(Y, path = paste0(bluecloud.wd,"/data/Y.feather"))
  
  # --- 13. Extract the vector of station names to keep track of the cross validation splits
  ID <- target %>% 
    dplyr::select(Station)
  ID <- sort(ID$Station)
  write_feather(data.frame(Station = ID), path = paste0(bluecloud.wd,"/data/Station_ID.feather"))
  
  # --- 14. Correspondence analysis
  # This analysis is done to assign each non-escoufier-selected target to the escoufier-selected one
  # that has the nearest observed abundance pattern
  # 14.1. Perform the correspondance analysis
  res_ca <- CA(Y0, graph = FALSE, ncp = 50)
  ndim_ca <- which(res_ca$eig[,3]>80)[1]
  dist_ca <- dist(res_ca$col$coord[,1:ndim_ca]) %>% as.matrix()
  dist_ca <- dist_ca[,1:ncol(Y)]
  
  # 14.2. Quantify dominance of each CC in terms of observed relative abundance
  tmp <- apply(Y0,2,sum)
  
  # 14.3. Concatenate the objects in a dataframe
  nn_ca <- data.frame(CC = dimnames(dist_ca)[[1]],
                      pos_CC = 1:nrow(dist_ca),
                      nn_CC = apply(dist_ca, 1, function(x){x = names(which(x == min(x)))}),
                      pos_nn_CC = apply(dist_ca, 1, function(x){x = which(x == min(x))}),
                      sum_CC = tmp,
                      scale_CC = tmp/tmp[apply(dist_ca, 1, function(x){x = which(x == min(x))})])
  
  # --- 15. Close connection to database and inform the user
  print(paste("Number of stations :", nrow(X)))
  print(paste("Number of environmental features :", ncol(X)))
  print(paste("Number of gene/cluster targets :", ncol(Y)))
  
  return(list(X=as.data.frame(X), Y=as.data.frame(Y), Y0=as.data.frame(Y0), CC_desc = as.data.frame(CC_desc), e = e, nn_ca = nn_ca))
  
  # --- Close connection
  dbDisconnect(db)
} # end function
