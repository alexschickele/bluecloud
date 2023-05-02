
# ===================== PART 1 : querying database =============================
# Extract the necessary data for the model according to the pre-defined filters
# i.e. number of genes and station per clusters

query_analysis <- function(bluecloud.wd = bluecloud_dir,
                           EXCLUDE = exclude,
                           KEGG_m = 165:172,
                           CLUSTER_SELEC = list(N_CLUSTERS = 30, MIN_GENES = 5, EXCLUSIVITY_R = 1),
                           ENV_METRIC = c("mean","sd","dist","bathy"),
                           relative = TRUE){
  
  # --- For local database
  # db <- dbConnect(RSQLite::SQLite(), paste0(bluecloud.wd, "/omic_data/",FILTER,"_DB_clean.sqlite"))
  db <- dbConnect(RSQLite::SQLite(), paste0(bluecloud.wd, "/omic_data/Picoeuk_DB_clean.sqlite"))
  
  # --- 1. Filter "data" by "cluster_sort"
  query <- dbGetQuery(db, paste0("SELECT CC FROM kegg_sort WHERE kegg_ko LIKE '%",
                                 paste(KEGG_m, collapse = "%' OR kegg_ko LIKE '%"), "%'")) %>%
    unique() %>% 
    inner_join(tbl(db, "kegg_sort"), copy = TRUE) %>%
    mutate(exclusivity = str_count(kegg_ko, paste(c("-","NA",KEGG_m), collapse = "|"))/n_ko) %>%
    dplyr::group_by(CC) %>% 
    dplyr::summarise(max_kegg = max(n_kegg, na.rm = TRUE), max_mod = max(n_mod, na.rm = TRUE), min_exl = min(exclusivity, na.rm = TRUE)) %>% 
    inner_join(tbl(db, "cluster_sort"), copy = TRUE) %>% 
    filter(n_genes >= !!CLUSTER_SELEC$MIN_GENES)
  query_check <- query
  
  target <- query %>% 
    dplyr::select("CC", "min_exl") %>% 
    inner_join(tbl(db, "data"), copy = TRUE) %>% 
    dplyr::select(c("Genes", "CC", "readCount", "Station", "Filter","Longitude", "Latitude", "KEGG_ko","KEGG_Module","Description", "Phylum","Class", "Genus", "min_exl"))
  
  target <- target %>% mutate(MAG = gsub("(.*_){2}(\\d+)_.+", "\\2", Genes))  %>% 
    dplyr::filter(!grepl(EXCLUDE, Class))
  
  # --- 1b. Excluding taxa by expert knowledge if needed
  target <- tbl(db, "target")
  
  # --- 2. Get cluster functional description
  
  CC_desc <- target %>% 
    inner_join(tbl(db, "cluster_sort"), copy = TRUE) %>% 
    group_by(CC, unknown_rate) %>% 
    summarise(kegg_ko = paste(unique(KEGG_ko), collapse = ", "),
              kegg_module = paste(unique(KEGG_Module), collapse = ", "),
              desc = paste(unique(Description), collapse = ", "),
              phylum = paste(unique(Phylum), collapse = ", "),
              class = paste(unique(Class), collapse = ", "),
              genus = paste(unique(Genus), collapse = ", "),
              mag = paste(unique(MAG), collapse = ", ")) %>% 
    inner_join(query_check)
    
  # ---  3.  Initialization
  # Color palette for superposed plots
  pal <- alpha(c("white","red","blue","black"), c(0.5,0.3,0.3,1))
  dens <- c(NA,NA,NA,20)
  
  # CC_desc tables with different filters
  list_all <- list(one = list(CC_desc = CC_desc),
                   two = list(CC_desc = CC_desc %>% 
                                filter(min_exl >= 0.5)),
                   three = list(CC_desc = CC_desc %>% 
                                  filter(min_exl >= 1)),
                   four = list(CC_desc = CC_desc %>% 
                                 filter(min_exl >= 1 & n_station >= 10)))
  
  # --- 4. Building distribution tables in list_all
  for(i in 1:length(list_all)){
    # Extracting Taxonomic, Functional and MAG annotations
    CC_desc_tmp <- list_all[[i]]$CC_desc
    factor_raw <- list(CC_desc_tmp$kegg_ko, CC_desc_tmp$class, CC_desc_tmp$mag)
    if(i==1){
      factor_names <- lapply(factor_raw, function(x){x <- gsub(x, pattern = " ", replacement = "") %>% 
        strsplit(split = ",") %>% 
        unlist() %>% unique()
      if(length(which(x == "-" | x == "NA")) > 0){
        x <- x[-which(x == "-" | x == "NA")]
      } else {x <-  x}
      })
      factor_names[[1]] <- kegg_m[paste0("ko:",kegg_m)%in%factor_names[[1]]]
    } # if: we want the baseline list of names only

    
    #  Building matrices
    df_list <- list()
    for(k in 1:length(factor_names)){
      df <- matrix(0, nrow = nrow(CC_desc_tmp), ncol = length(factor_names[[k]]), dimnames = list(CC_desc_tmp$CC, factor_names[[k]]))
      for(m in 1:dim(df)[[1]]){
        for(n in 1:dim(df)[[2]]){
          if(str_detect(factor_raw[[k]][m], factor_names[[k]][n])==TRUE){df[m,n] <- df[m,n]+1}
        }  # m
      } # n
      df <- apply(df,2,sum)
      # df <- df/sum(df)
      df_list[[k]] <- df
    } # k
    list_all[[i]][["df_list"]] <- df_list
  } # i
  
  # --- 5. Plotting
  par(mfrow = c(2,2), mar = c(8,5,3,1))
  
  # Sum of reads plot
  list_all[[1]][["hist"]][["sum_reads"]] <- hist(log1p(list_all[[1]]$CC_desc$sum_reads), breaks = seq(0,10,0.5), col = pal[1], main = "sum of Reads (log1p)", xlab = "", xlim = c(0,7))
  abline(h = seq(0,250,50), col = "black")
  for(i in 2:length(list_all)){
    list_all[[i]][["hist"]][["sum_reads"]] <- hist(log1p(list_all[[i]]$CC_desc$sum_reads), breaks = seq(0,10,0.5), col = pal[i], xlab = "",
                                                   density = dens[i], border = NA, add = TRUE)
    print("Sum of reads significancy test :")
    print(chisq.test(list_all[[1]][["hist"]][["sum_reads"]]$counts,list_all[[i]][["hist"]][["sum_reads"]]$counts))
  }
  
  # Nb. of genes plot
  list_all[[1]][["hist"]][["n_genes"]] <- hist(log1p(list_all[[1]]$CC_desc$n_genes), breaks = seq(0,10,0.5), col = pal[1], main = "Nb. of genes (log1p)", xlab = "", xlim = c(0,9))
  abline(h = seq(0,250,50), col = "black")
  for(i in 2:length(list_all)){
    list_all[[i]][["hist"]][["n_genes"]] <- hist(log1p(list_all[[i]]$CC_desc$n_genes), breaks = seq(0,10,0.5), col = pal[i], xlab = "", xlim = c(0,9),
                                                 density = dens[i], border = NA, add = TRUE)
    print("Nb. of genes significancy test :")
    print(chisq.test(list_all[[1]][["hist"]][["n_genes"]]$counts,list_all[[i]][["hist"]][["n_genes"]]$counts))
  }
  
  # Functional KO plot
  barplot(list_all[[1]]$df_list[[1]], las = 2, cex.names = 0.6, col = pal[1], main = "KO composition", ylab = "Frequency")
  abline(h = seq(0,100,20), col = "black")
  for(i in 2:length(list_all)){
    barplot(list_all[[i]]$df_list[[1]], las = 2, cex.names = 0.6, col = pal[i], density = dens[i], border = NA,
            add = TRUE, axisnames = FALSE, axes = FALSE)
    print("Functional KO significancy test :")
    print(chisq.test(list_all[[1]]$df_list[[1]],list_all[[i]]$df_list[[1]]))
  }
  
  # Taxonomic class plot
  barplot(list_all[[1]]$df_list[[2]], las = 2, cex.names = 0.6, col = pal[1], main = "Taxo. class composition", ylab = "Frequency")
  abline(h = seq(0,100,20), col = "black")
  for(i in 2:length(list_all)){
    barplot(list_all[[i]]$df_list[[2]], las = 2, cex.names = 0.6, col = pal[i], density = dens[i], border = NA,
            add = TRUE, axisnames = FALSE, axes = FALSE)
    print("Taxonomic class significancy test :")
    print(chisq.test(list_all[[1]]$df_list[[2]],list_all[[i]]$df_list[[2]]))
  }
  
  # --- Close connection
  dbDisconnect(db)
} # end function
