# This script is to draw heatmap for all TFs in different tissues.
# Tissue expression was from Stelpflug et al., 2015 Plant Genome paper.
# Do mean calculation for each tissue, only keep one value each tissue.
# So one column each tissue. Plot 1409 TFs that expressed in all four tissues.
# Files in: 1) genes in each tissue: cpm_tissue.RData
#           2) Atlas expression data: PlantGenomeS1_fourTissue_only.txt
###############################################################################


# load library and setwd
library(tidyverse)
library(pheatmap)
library(RColorBrewer)
load("~/projects/NTWK/tissue_ntwk/data/right_seed/cpm_tissue.RData")

# Find TFs in each tissue type.
tf_leaf <- intersect(rownames(cpm_leaf), tf_all)
tf_root <- intersect(rownames(cpm_root), tf_all)
tf_sam <- intersect(rownames(cpm_sam), tf_all)
tf_seed <- intersect(rownames(cpm_seed), tf_all)

# Find TFs in four tissues.
tf_four <- Reduce(intersect, list(tf_leaf, tf_root, tf_sam, tf_seed))

# clean
rm(cpm_leaf, cpm_root, cpm_sam, cpm_seed)
rm(tf_leaf, tf_root, tf_seed, tf_sam)

###############################################################################
# PLOT 
# Plot mean value for each tissue. One value for each tissue.
# Cut rows to 15 clusters. 15 = C(4,1) + C(4,2) + C(4,3) +C(4,4)
################################################################################

# read expression table
pg_expr <- read_tsv("data/PlantGenomeS1_fourTissue_only.txt", col_names = T) %>% 
  filter(geneid %in% tf_four)

# calculate mean value and prepare matrix.
mydata <- cbind(MeanI = rowMeans(pg_expr[2:19],na.rm = T))
mydata <- cbind(mydata, MeanI = rowMeans(pg_expr[20:38],na.rm = T))
mydata <- cbind(mydata, MeanI = rowMeans(pg_expr[39:40],na.rm = T))
mydata <- cbind(mydata, MeanI = rowMeans(pg_expr[41:64],na.rm = T))
colnames(mydata) <- c("root", "leaf", "SAM", "seed")

# remove TFs with sd = 0.
expr_sd <- apply(mydata,1, sd)
pg_expr <- mydata[(expr_sd != 0),]

# Plot the heatmap. save as PNG.
# PDF's color is different, png is ok.
png(filename = "results/TF_1409_expr_four_tissue_cutrow_15.png", 
    width = 7, height = 7, units = "in", res = 300)
pheatmap(pg_expr, scale = "row", cluster_cols = T, treeheight_row = 0,
         treeheight_col = 20, fontsize_col = 20, cutree_rows = 15,
         color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100))
dev.off()




