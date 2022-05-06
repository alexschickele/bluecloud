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


# --- Testing special clustering
# Clustering func vs mag correlation with taxonomic information on mag cluster
library(dendextend)

# --- Calculate correlation matrix and clusterings
cor_mat <- cor(mag_data, func_data, use = 'pairwise.complete.obs')
func_clust <- hclust(dist(t(cor_mat)), method = "ward.D2")
mag_clust <- hclust(dist(cor_mat), method = "ward.D2")
mag_group <- cutree(mag_clust, k = 8)

# --- Compute list of taxo labels
taxo_lab <- NULL
for(j in 1:length(mag_clust$labels)){
  tmp <- CC_desc_e$phylum[which(str_detect(CC_desc_e$mag, mag_clust$labels[j])==TRUE)] %>% 
    strsplit(split = ",") %>% unlist() %>% unique()
  if(length(tmp) > 1){tmp <- "Multiple"}
  taxo_lab <- c(taxo_lab, tmp)
}

# --- Color palette for labels
taxo_pal <- colorRampPalette(col = rev(brewer.pal(10,"Spectral")))(length(unique(taxo_lab)))
names(taxo_pal) <- unique(taxo_lab)
taxo_pal[c(2)] <- "white"

# --- Creating pie charts
par(mfrow = c(3,3), mar = c(3,3,3,3))
for(j in 1:max(mag_group)){
  tmp <- taxo_lab[which(mag_group[mag_clust$order] == j)] %>% as.factor() %>% summary()
  pie(tmp, labels = names(tmp), main = paste("group", j),
      col = taxo_pal[names(tmp)])
}

# --- Plot heatmap
heatmap(cor_mat, Rowv = as.dendrogram(mag_clust), Colv = as.dendrogram(func_clust),
        col = brewer.pal(9, "RdBu"), cexCol = 0.8)


toto <- color_branches(dend = as.dendrogram(mag_clust), col = taxo_pal[as.factor(taxo_lab)])
labels_cex(toto) <- 0.5
labels_colors(toto) <- taxo_pal[as.factor(taxo_lab)]
labels(toto) <- taxo_lab
plot(toto)
points(x = seq(1:length(mag_clust$labels)), y = rep(0, length(mag_clust$labels)), pch = 21, bg = taxo_pal[taxo_lab])


toto <- heatmap(cor_mat)
