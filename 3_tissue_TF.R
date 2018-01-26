# This script is to find TFs and Genes expressed in each tissue.
# The criteria is defined in "2_expr_norm.R". Genes need to have CPM >1 in more
# than 10% of total libraries of that tissue.

# File in: 1. Four matix of normalized log2(CPM+1) values for each tissue.;
#              (cpm_leaf, cpm_root, cpm_sam and cpm_seed)
#          2. All maize TFs defined by GRASSIUS.
#
# File out: 1. Four txt files with TFs in each tissue:
#              (tf_leaf_specifc_10%.txt, tf_root_specifc_10%.txt,
#               tf_sam_specifc_10%.txt, tf_seed_specifc_10%.txt)
#            2.  Four txt files with Genes in each tissue:
#              (genes_leaf_specifc_10%.txt, genes_root_specifc_10%.txt,
#               genes_sam_specifc_10%.txt, genes_seed_specifc_10%.txt)
################################################################################


library(tidyverse)
setwd("~/projects/NTWK/tissue_ntwk/")
load("data/right_seed/cpm_tissue.RData")

# read all TF genes and keep only genes have version 3 IDs.
tf_all <- read_tsv("data/maize_all_TF_Grassius.txt",col_names = T)

tf_all <- intersect(tf_all$Geneid,gene_name_all) # 

# Find TFs in each tissue type.
tf_leaf <- intersect(rownames(cpm_leaf), tf_all)
tf_root <- intersect(rownames(cpm_root), tf_all)
tf_sam <- intersect(rownames(cpm_sam), tf_all)
tf_seed <- intersect(rownames(cpm_seed), tf_all)

tf_four <- Reduce(intersect, list(tf_leaf, tf_root, tf_sam, tf_seed))


# write TFs in each tissue.
write.table(tf_leaf, file = "data/right_seed/tf_leaf_specifc_10%.txt",quote = F, 
            sep = "\t", row.names = F, col.names = "tf_leaf")
write.table(tf_root, file = "data/right_seed/tf_root_specifc_10%.txt",quote = F, 
            sep = "\t", row.names = F, col.names = "tf_root")
write.table(tf_sam, file = "data/right_seed/tf_sam_specifc_10%.txt",quote = F, 
            sep = "\t", row.names = F, col.names = "tf_sam")
write.table(tf_seed, file = "data/right_seed/tf_seed_specifc_10%.txt",quote = F, 
            sep = "\t", row.names = F, col.names = "tf_seed")

save.image(file = "data/right_seed/cpm_tissue.RData")

# write all genes in each tissue.
write.table(row.names(cpm_leaf), file = "data/right_seed/genes_leaf_10%.txt", 
            quote = F, sep = "\t", row.names = F, col.names = "genes_leaf")

write.table(row.names(cpm_root), file = "data/right_seed/genes_root_10%.txt", 
            quote = F, sep = "\t", row.names = F, col.names = "genes_root")

write.table(row.names(cpm_sam), file = "data/right_seed/genes_sam_10%.txt", 
            quote = F, sep = "\t", row.names = F, col.names = "genes_sam")

write.table(row.names(cpm_seed), file = "data/right_seed/genes_seed_10%.txt", 
            quote = F, sep = "\t", row.names = F, col.names = "genes_seed")
