# This script calculates weigh matrix from GENIE3 on four tissue expression matrix.
# Then save the top 1 million edges from the GRN,
#
# File in: 1. Four matix of normalized log2(CPM+1) values for each tissue.
#              (cpm_leaf, cpm_root, cpm_sam and cpm_seed)
#          2. GENIE3 source code (R/C version)
#
# File out: 1. Four weigh matrix as RData: (wm_leaf_tf.RData, wm_root_tf.RData, 
#                wm_sam_tf.RData, wm_seed_tf.RData)
#           2. Four top one million edges as link lists:
#                 link_list_four_tissue_1million.RData
#           3. The same link lists file as tab-delimited files:
#                 (ll_leaf.txt, ll_root.txt, ll_sam.txt, ll_seed.txt)
################################################################################


# load library and tissue specific expression matrix.
library(tidyverse)
source("GENIE3.R")
load("data/right_seed/cpm_tissue.RData")
setwd("~/projects/NTWK/tissue_ntwk/")

# Calculate weigh matrix.

# weigh matrix (wm) for root
ptm <- proc.time()
wm_root <- GENIE3(cpm_root, ncores = 22, seed = 123,
                  regulators = tf_root, K = "sqrt")
proc.time() - ptm
save(wm_root,file = "wm_root_tf.RData")

# weigh matrix (wm) for leaf
ptm <- proc.time()
wm_leaf <- GENIE3(cpm_leaf, ncores = 22, seed = 123,
                  regulators = tf_leaf, K = "sqrt")
proc.time() - ptm # 28694s
save(wm_leaf,file = "wm_leaf_tf.RData")

# weigh matrix (wm) for sam
ptm <- proc.time()
wm_sam <- GENIE3(cpm_sam, ncores = 22, seed = 123,
                  regulators = tf_sam, K = "sqrt")
proc.time() - ptm # 40789s
save(wm_sam,file = "wm_sam_tf.RData")

# weigh matrix (wm) for seed
ptm <- proc.time()
wm_seed <- GENIE3(cpm_seed, ncores = 22, seed = 123,
                 regulators = tf_seed, K = "sqrt")
proc.time() - ptm # 5099s
save(wm_seed,file = "data/right_seed/wm_seed_tf.RData")

# Get links list (ll) for four tissue types.
# Did for 100,000, 1,000,000 and 10,000,000.
top_links <- 10000000

ll_leaf <- get.link.list(wm_leaf, report.max = top_links)
ll_root <- get.link.list(wm_root, report.max = top_links) 
ll_sam <- get.link.list(wm_sam, report.max = top_links)
ll_seed <- get.link.list(wm_seed, report.max = top_links)


# save data. As binary and as txt to import to Cytoscape.
save(ll_leaf, ll_root, ll_sam, ll_seed, 
     file = "data/right_seed/link_list_four_tissue_10million.RData")

save(ll_leaf, ll_root, ll_sam, ll_seed, 
     file = "data/right_seed/link_list_four_tissue_1million.RData")


# save(ll_leaf, ll_root, ll_sam, ll_seed, 
#     file = "link_list_four_tissue_100thousand.RData")

# save the tab-delimited txt.
write.table(ll_root, file = "ll_root.txt", quote = F, sep = "\t", row.names = F)
write.table(ll_leaf, file = "ll_leaf.txt", quote = F, sep = "\t", row.names = F)
write.table(ll_seed, file = "data/right_seed/ll_seed.txt", quote = F, 
            sep = "\t", row.names = F)
write.table(ll_sam, file = "ll_sam.txt", quote = F, sep = "\t", row.names = F)


# Remove unecessary data
rm(cpm_leaf, cpm_sam, cpm_seed, cpm_root)
rm(ptm, top_links, GENIE3, get.link.list, read.expr.matrix)

