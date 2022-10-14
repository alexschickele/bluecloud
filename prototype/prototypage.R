
# ==============================================================================
# TARA Stations map
# ==============================================================================

# --- For local database
db <- dbConnect(RSQLite::SQLite(), paste0(bluecloud.wd, "/omic_data/",FILTER,"_DB.sqlite"))

# --- Collect data
tara_station <- tbl(db, "sum_station") %>% 
  inner_join(tbl(db, "locs")) %>% 
  select(-sum_reads) %>% 
  collect()

features <- stack(paste0(data_dir, "/features"))

# --- Plot
library(maps)
pal = colorRampPalette(brewer.pal(9,"Spectral"))(100) %>% rev()
plot(features[["temperaturemean"]], col = pal)
map("world", add = TRUE)
points(x = tara_station$Longitude, y = tara_station$Latitude, pch = 19, cex = 0.8)

# --- Barplots
par(mfrow = c(2,1))
hist(query$X$temperaturemean, breaks = 10, freq = TRUE, col = "gray50", xlab = "temperature")
hist(query$X$oxygenmean, breaks = 10, freq = TRUE, col = "gray50", xlab = "oxygen")

# ==============================================================================
# Investigating the HEXANAUPLIA & co case
# ==============================================================================

# all classes tag
unique(CC_desc_e$class)

# look for one suspect
suspect <- "Oomycota"
id <- grep(pattern = suspect, CC_desc_e$class)
CC_suspect <- CC_desc_e[id,]
CC_suspect[c("CC", "class", "kegg_ko", "mag")]

# look for the corresponding genes
CC_id <- CC_suspect$CC
db <- dbConnect(RSQLite::SQLite(), paste0(bluecloud.wd, "/omic_data/",FILTER,"_DB.sqlite"))

query <- dbGetQuery(db, paste0("SELECT * FROM data WHERE CC LIKE '%",
                               paste(CC_id, collapse = "%' OR CC LIKE '%"), "%'")) 

unique(query$Class)
unique(query$KEGG_ko)




# ==============================================================================
# Which enzyme is detected in which taxa ?
# ==============================================================================
# Adapted to Nitrogen case for now

# --- creating storage table
enz_detect <- matrix(data = NA,
                     nrow = length(kegg_m),
                     ncol = length(factor_names[[3]]),
                     dimnames = list(kegg_m, factor_names[[3]]))

# --- iteratively filling up the table
for(i in 1:nrow(enz_detect)){
  for(j in 1:ncol(enz_detect)){
    id <- CC_desc_e$pos_CC[which(str_detect(CC_desc_e$kegg_ko, kegg_m[i])==TRUE & str_detect(CC_desc_e$class, factor_names[[3]][j])==TRUE)]
    scale <- CC_desc_e$sum_CC[id] %>% sum()
    
    if(length(id) != 0){enz_detect[i,j] <- scale}
  
  } # j col
} # i row

enz_detect <- enz_detect/max(enz_detect, na.rm = TRUE)

# --- Plot
par(mar = c(10,5,1,6))
image(t(log(enz_detect)), col = pal, axes = FALSE)
box()

axis(side = 1, at = seq(0,1, length.out = ncol(enz_detect)), labels = colnames(enz_detect), las = 2, cex.axis = 0.8)
abline(v = seq(0.023,1.023, length.out = ncol(enz_detect)))

axis(side = 2, at = seq(0,1, length.out = nrow(enz_detect)), labels = rownames(enz_detect), las = 2, cex.axis = 0.8)
abline(h = seq(0.05,1.05, length.out = nrow(enz_detect)))

axis(side = 4, at = c(0.07, 0.24, 0.5, 0.69, 0.88), labels = c("N. Transferases", "Ass. N. Red.", "Diss. N. Red.", "NH3 to Urea", "NH3 to GLU"), 
     las = 2, tick = FALSE, cex.axis = 0.8)
abline(h = c(0.14, 0.34, 0.63, 0.76), lwd = 4)

# --- Legend plot
image(1:100, 1, as.matrix(1:100), col = pal)

# ==============================================================================
# ==============================================================================

library(leaflet)
library(raster)
library(tidyverse)

# marche po
range_r <- getValues(proj$proj[[1]]) %>% range(na.rm = TRUE)

value_r <- proj$proj[[1]] %>% 
  getValues()
value_r <- value_r[!is.na(value_r)] %>% 
  sort()

pal <- proj$col_matrix
pal <- pal[value_r] %>% sort() %>% rev()

leaflet() %>%
  addRasterImage(proj$proj[[1]], colors = pal, opacity = 1)

# marche
pal <- proj$col[[1]]
plot(proj$proj[[1]], col = pal)

# autre manière de définir le vecteur de couleur à partir de la matrice de base:
# i.e. une cellule = une valeur x incertitude; 
# les valeurs du raster sont l'index de la matrice correspondant à la combinaison valeur x incertitude
pal <- proj$col_matrix[range_r[1]:range_r[2]]




# ==============================================================================
# Uncertainty propagation in the aggregated models
one <- runif(100, 0, 20)
two <- runif(100, 80, 100)

mean(one)
sd(one)/mean(one)
mean(two)
sd(two)/mean(two)

three <- c(one, two)
mean(three)
sd(three)/mean(three)




# ==============================================================================
# Functional vs env plot with raw data
par(mfrow = c(3,3))
for(j in 1:length(plot_list)){
  # Extract functional data
  id <- which(str_detect(CC_desc_e$kegg_ko, plot_list[[j]])==TRUE)
  x <- NULL
  y <- NULL
  for(k in 1:length(id)){
    x <- cbind(x, query$X$CHLsd)
    y <- cbind(y, query$Y0[,id[k]])
  }
  
  plot(x, y, ylim = c(0, 0.5), main = names(plot_list)[j], pch = 20, col = alpha("black", 0.2))
  grid(col = "gray30")
  lines(x = x[order(x[,1]),1], y = apply(y, 1, mean)[order(x[,1])], col = "red")

} # for j



# ==============================================================================
# --- Testing different between sum of columns in Y and y_hat_m
sum_Y <- apply(query$Y, 2, max)
sum_Y <- sum_Y/max(sum_Y)

sum_y_hat <- apply(eval$y_hat, 2, function(x) (x = max(x, na.rm = TRUE)))
sum_y_hat <- sum_y_hat/max(sum_y_hat, na.rm = TRUE)

# column plot
barplot(sum_Y, las = 2, cex.names = 0.5, col = alpha("red", 0.3))
barplot(sum_y_hat, las = 2, cex.names = 0.5, col = alpha("blue",0.3), add = TRUE)

# --- Testing different between sum of columns in Y and y_hat_m
sum_Y <- apply(query$Y, 1, sum)
sum_Y <- sum_Y/max(sum_Y)

sum_y_hat <- apply(eval$y_hat, 1, function(x) (x = sum(x, na.rm = TRUE)))
sum_y_hat <- sum_y_hat/max(sum_y_hat, na.rm = TRUE)

# line plot
barplot(sum_Y, las = 2, cex.names = 0.5, col = alpha("red", 0.3))
barplot(sum_y_hat, las = 2, cex.names = 0.5, col = alpha("blue",0.3), add = TRUE)


# ==============================================================================
# --- Testing special clustering
# ==============================================================================

# Clustering func vs mag correlation with taxonomic information on mag cluster
library(dendextend)
library(vegan)

# --- Removing uncertainty
r0 <- stack(paste0(data.wd,"/features"))[[1]]
mag_data0 <- apply(mag_data, c(1,3), mean)
func_data0 <- apply(func_data, c(1,3), mean)

# --- Calculate correlation matrix and clusterings
cor_mat <- cor(mag_data0, func_data0, use = 'pairwise.complete.obs')
func_clust <- hclust(dist(t(cor_mat)), method = "ward.D2")
mag_clust <- hclust(dist(cor_mat), method = "ward.D2")
mag_group <- cutree(mag_clust, k = 6)

# --- Compute list of taxo labels
taxo_lab <- NULL
for(j in 1:length(mag_clust$labels)){
  tmp <- CC_desc_e$class[which(str_detect(CC_desc_e$mag, mag_clust$labels[j])==TRUE)] %>% 
    strsplit(split = ",") %>% unlist() %>% unique()
  if(length(tmp) > 1){tmp <- "Multiple"}
  taxo_lab <- c(taxo_lab, tmp)
}

# --- Compute scale corresponding to each mag
mag_scale <- apply(mag_data0, 2, function(x) (x = sum(x, na.rm = TRUE)))

# --- Color palette for labels
taxo_pal <- colorRampPalette(col = rev(brewer.pal(10,"Spectral")))(length(unique(taxo_lab)))
names(taxo_pal) <- unique(taxo_lab)
taxo_pal[c(2)] <- "white"

if(scaled == TRUE){pdf(paste0(bluecloud_dir,"/output/", output_dir, "/MAG_correlations_scaled.pdf"))
}else{pdf(paste0(bluecloud_dir,"/output/", output_dir, "/MAG_correlations.pdf"))}


# --- Creating pie charts
par(mfrow = c(3,2), mar = c(3,3,3,3))
for(j in 1:max(mag_group)){
  tmp <- data.frame(lab = taxo_lab[which(mag_group[mag_clust$order] == j)],
                    scale = mag_scale[which(mag_group[mag_clust$order] == j)]) %>% 
    group_by(lab) %>% 
    summarize(sum_scale = sum(scale)) %>% 
    mutate(sum_scale = sum_scale/sum(sum_scale))
  
  pie(tmp$sum_scale, labels = tmp$lab, col = taxo_pal[tmp$lab],
      main = paste("group", j, "\n", "H:", round(diversity(tmp[which(tmp$lab != "Multiple"),]$sum_scale), 2)))
}

# --- Creating barplot version with top 5
par(mfrow = c(3,2), mar = c(3,10,3,3))
for(j in 1:max(mag_group)){
  tmp <- data.frame(lab = taxo_lab[which(mag_group[mag_clust$order] == j)],
                    scale = mag_scale[which(mag_group[mag_clust$order] == j)]) %>% 
    group_by(lab) %>% 
    summarize(sum_scale = sum(scale)) %>% 
    mutate(sum_scale = sum_scale/sum(sum_scale))
  tmp <- tail(tmp[order(tmp$sum_scale),], 5)
  barplot(tmp$sum_scale, names = tmp$lab, las = 2, col = "black", horiz = TRUE,  xlim = c(0,1),
          main = paste("group", j))
  abline(v = seq(0,1,0.1), lty = "dotted")
}


# --- Create mag group maps ?
par(mfrow = c(3,2), mar = c(3,3,3,3))
for(j in 1:max(mag_group)){
  tmp <- mag_data0[,which(mag_group[mag_clust$order] == j)] %>% 
    apply(1, sum)
  mag_r <- setValues(r0, tmp)
  mag_r <- synchroniseNA(stack(r0, mag_r))[[-1]]
  plot(mag_r, col = pal, main = paste("group", j))
}

# --- Plot heatmap
heatmap(cor_mat, Rowv = as.dendrogram(mag_clust), Colv = as.dendrogram(func_clust),
        col = rev(brewer.pal(9, "RdBu")), cexCol = 0.8)

dev.off()


