



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


