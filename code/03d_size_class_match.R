
# ==============================================================================
# GGMM & GGZZ correlation - misc. analysis
# ==============================================================================

# This code computes a set of metric to ensure that the GGMM and GGZZ size classes
# do not present major differences in RUBISCO and C4-enzymes' read composition
# across common stations.

# NECESSARY: Run 01_query_BC until step 6 (included)
# Now we get target just before the average across
toto <- target %>% collect()

# --- 1. Extract stations
# Get stations where we have GGMM
GGMM_locs <- toto %>% 
  dplyr::select(Station, Filter) %>% 
  dplyr::filter(Filter == 'GGMM')

# Get stations where we have GGZZ
GGZZ_locs <- toto %>% 
  dplyr::select(Station, Filter) %>% 
  dplyr::filter(Filter == 'GGZZ')

# Get common stations
common_locs <- GGZZ_locs %>% 
  inner_join(GGMM_locs, by = "Station")

# --- 2. Extract read data
# Read matrix from GGMM
GGMM_reads <- toto %>% 
  dplyr::filter(Filter == "GGMM") %>% 
  dplyr::filter(grepl(paste(common_locs$Station, collapse = "|"), Station)) %>% 
  dplyr::select(-Station, -Filter, -Longitude, -Latitude, -sum_reads) %>% 
  apply(1, function(x){if(sum(x)>0){x = x/sum(x, na.rm = TRUE)} else {x = x}}) %>% 
  aperm(c(2,1))

# Read matrix from GGZZ
GGZZ_reads <- toto %>% 
  dplyr::filter(Filter == "GGZZ") %>% 
  dplyr::filter(grepl(paste(common_locs$Station, collapse = "|"), Station)) %>% 
  dplyr::select(-Station, -Filter, -Longitude, -Latitude, -sum_reads) %>% 
  apply(1, function(x){if(sum(x)>0){x = x/sum(x, na.rm = TRUE)} else {x = x}}) %>% 
  aperm(c(2,1))

# --- 3. Are the read composition significantly different ?
# Mantel test
library(vegan)
mantel(dist(GGMM_reads), dist(GGZZ_reads))
# mantel(dist(GGMM_reads), dist(GGMM_reads)) # Control

# Simple GGMM vs GGZZ plot
plot(GGMM_reads,GGZZ_reads)
grid()
abline(coef = c(0,1), col = "red")

# --- 4. Supplementary metrics
# Protein functional cluster detection rate
tmp1 <- which(GGMM_reads != 0) %>% length()
message(paste("Number of Protein Functional Clusters detected across stations for the 0.8 - 5 μm class (GGMM):", tmp1))
tmp2 <- which(GGZZ_reads != 0) %>% length() 
message(paste("Number of Protein Functional Clusters detected across stations for the 0.8 - 2000 μm class (GGZZ):", tmp2))
tmp3 <- which(GGZZ_reads != 0 & GGMM_reads != 0) %>% length()
message(paste("Number of Protein Functional Clusters detected across stations, common for both size classes:", tmp3))
message(paste("Corresponding detection rate:", tmp3/tmp2))

# Save the union id's
id <- which(GGZZ_reads != 0 & GGMM_reads != 0)

# Protein functional cluster read rate
tmp1 <- which(GGMM_reads != 0) %>% sum()
message(paste("Number of Protein Functional Clusters reads across stations for the 0.8 - 5 μm class (GGMM):", tmp1))
tmp2 <- which(GGZZ_reads != 0) %>% sum()
message(paste("Number of Protein Functional Clusters reads across stations for the 0.8 - 2000 μm class (GGZZ):", tmp2))
tmp3 <- which(GGZZ_reads != 0 & GGMM_reads != 0) %>% sum()
message(paste("Number of Protein Functional Clusters reads across stations, common for both size classes:", tmp3))
message(paste("Corresponding read rate:", tmp3/tmp2))

# Correlation
tmp <- cor(GGMM_reads[id], GGZZ_reads[id], method = "pearson")
message(paste("Pearson correlation on the reads of the enzymes detected in both size classes:", tmp))
tmp <- cor(GGMM_reads[id], GGZZ_reads[id], method = "spearman")
message(paste("Spearman correlation on the reads of the enzymes detected in both size classes:", tmp))

