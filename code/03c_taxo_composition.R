

# ==============================================================================
# --- Testing special clustering
# Objective: get the taxonomic annotation composition for each enzyme map
# based on the corresponding MAGs
# ==============================================================================

# --- Load additional libraries
library(dendextend)
library(vegan)

# --- Initial parameters
bluecloud.wd = bluecloud_dir

# --- Build table and query to get the MAG - taxa correspondence
db <- dbConnect(RSQLite::SQLite(), paste0(bluecloud.wd, "/omic_data/",FILTER,"_DB.sqlite"))

tmp <- tbl(db, "data") %>% 
  dplyr::select(c("Genes", "Class")) %>% 
  collect()

taxomag <- tmp %>% mutate(MAG = gsub("(.*_){2}(\\d+)_.+", "\\2", Genes)) %>% 
  dplyr::select(c("MAG", "Class")) %>% 
  dplyr::distinct()

save(list = c("taxomag"), file = paste0(bluecloud.wd, "/data/taxomag.RData"))

# --- Load files
load(file = paste0(bluecloud.wd, "/data/taxomag.RData"))

# --- Removing uncertainty
r0 <- stack(paste0(data.wd,"/features"))[[1]]
mag_data0 <- mag_data
func_data0 <- apply(func_data, c(1,3), mean)

# --- Calculate correlation matrix and clusterings
cor_mat <- cor(mag_data0, func_data0, use = 'pairwise.complete.obs')
func_clust <- hclust(dist(t(cor_mat)), method = "ward.D2")
mag_clust <- hclust(dist(cor_mat), method = "ward.D2")
mag_group <- cutree(mag_clust, k = 3)

# --- Compute list of taxo labels
taxo_lab <- matrix(0, 
                   nrow = length(dimnames(taxo_data)[[3]]),
                   ncol = length(dimnames(mag_data)[[2]]))
rownames(taxo_lab) <- dimnames(taxo_data)[[3]]
colnames(taxo_lab) <- dimnames(mag_data)[[2]]

for(j in 1:dim(taxo_lab)[[2]]){
  tmp <- taxomag[which(taxomag$MAG==colnames(taxo_lab)[j]),]
  taxo_lab[grep(pattern = paste(tmp$Class, collapse = "|"), rownames(taxo_lab)),j] <- taxo_lab[grep(pattern = paste(tmp$Class, collapse = "|"), rownames(taxo_lab)),j]+1
}

# --- Compute scale corresponding to each mag
mag_scale <- apply(mag_data0, 2, function(x) (x = sum(x, na.rm = TRUE)))
taxo_lab <- apply(taxo_lab, 1, function(x) (x = x*mag_scale)) %>% t()

# --- Color palette for labels
taxo_pal <- colorRampPalette(col = rev(brewer.pal(10,"Spectral")))(length(rownames(taxo_lab)))
names(taxo_pal) <- rownames(taxo_lab)

# --- Plot all in a PDF
if(scaled == TRUE){pdf(paste0(bluecloud_dir,"/output/", output_dir, "/MAG_correlations_scaled.pdf"))
}else{pdf(paste0(bluecloud_dir,"/output/", output_dir, "/MAG_correlations.pdf"))}

# --- Creating pie charts
par(mfrow = c(3,2), mar = c(3,3,3,3))
for(j in 1:max(mag_group)){
  tmp <- taxo_lab[,which(mag_group == j)] %>% 
    apply(1, sum)
  tmp <- tmp/sum(tmp)
  
  pie(tmp, labels = names(tmp), col = taxo_pal[names(tmp)],
      main = paste("group", j, "\n", "H:", round(diversity(tmp), 2)))
}

# --- Creating barplot version with top 5
par(mfrow = c(3,2), mar = c(3,10,3,3))
for(j in 1:max(mag_group)){
  tmp <- taxo_lab[,which(mag_group == j)] %>% 
    apply(1, sum)
  tmp <- tmp/sum(tmp)
  tmp <- tail(tmp[order(tmp)], 5)
  barplot(tmp, names = names(tmp), las = 2, col = "black", horiz = TRUE,  xlim = c(0,1),
          main = paste("group", j))
  abline(v = seq(0,1,0.1), lty = "dotted")
}

# --- Create mag group maps ?
par(mfrow = c(3,2), mar = c(3,3,3,3))
for(j in 1:max(mag_group)){
  tmp <- mag_data0[,which(mag_group == j)] %>% 
    apply(1, sum)
  mag_r <- setValues(r0, tmp)
  mag_r <- synchroniseNA(stack(r0, mag_r))[[-1]]
  plot(mag_r, col = pal, main = paste("group", j))
}

# --- Plot heatmap
heatmap(cor_mat, Rowv = as.dendrogram(mag_clust), Colv = as.dendrogram(func_clust),
        col = rev(brewer.pal(9, "RdBu")), cexCol = 0.8)

dev.off()
