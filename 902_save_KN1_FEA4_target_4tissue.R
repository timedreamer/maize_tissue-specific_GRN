# This script is to save KN1 and FEA4 targets in four tissues
# top 1 million edges. The genes were used by GO analysis later, as well as 
# finding overlap and distinct genes.
#
# File in: 1) four tissue GRN (top 1million):
#             ll_leaf, ll_root, ll_sam and ll_seed. (Or use ll_four)
# 
# File out: 1) KN1 prediced targets in sam and seed:
#              kn1_sam_targets_758.txt, kn1_seed_targets_1044.txt
#           2) FEA4 prediced targets in leaf, root, sam and seed:
#              fea4_leaf_targets_283.txt, fea4_root_targets_55.txt
#              fea4_sam_targets_501.txt, fea4_seed_targets_28.txt
################################################################################


# load libraries and data.
setwd("~/projects/NTWK/tissue_ntwk/")
load("data/right_seed/link_list_four_tissue_1million.RData")
library(tidyverse)


# combine four link_list into one.
ll_leaf <- mutate(ll_leaf, tissue = "leaf")
ll_root <- mutate(ll_root, tissue = "root")
ll_sam <- mutate(ll_sam, tissue = "sam")
ll_seed <- mutate(ll_seed, tissue = "seed")

ll_four <- bind_rows(ll_leaf, ll_root, ll_sam, ll_seed)
rm(ll_leaf,ll_root, ll_sam, ll_seed)

# 2.1 write KN1 targets into files
gene <- "GRMZM2G017087"

kn1_sam <- ll_four %>% 
  filter(regulatory.gene == gene, tissue == "sam")

kn1_seed <- ll_four %>% 
  filter(regulatory.gene == gene, tissue == "seed")

write_tsv(x = kn1_sam, path = "results/kn1_sam_targets_758.txt",col_names = T)
write_tsv(x = kn1_seed, path = "results/kn1_seed_targets_1044.txt",col_names = T)

# 2.2 write FEA4 targets into files
gene <- "GRMZM2G133331"

fea4_leaf <- ll_four %>% 
  filter(regulatory.gene == gene, tissue == "leaf")

fea4_root <- ll_four %>% 
  filter(regulatory.gene == gene, tissue == "root")

fea4_sam <- ll_four %>% 
  filter(regulatory.gene == gene, tissue == "sam")

fea4_seed <- ll_four %>% 
  filter(regulatory.gene == gene, tissue == "seed")

write_tsv(x = fea4_leaf, path = "results/fea4_leaf_targets_283.txt",col_names = T)
write_tsv(x = fea4_root, path = "results/fea4_root_targets_55.txt",col_names = T)
write_tsv(x = fea4_sam, path = "results/fea4_sam_targets_501.txt",col_names = T)
write_tsv(x = fea4_seed, path = "results/fea4_seed_targets_28.txt",col_names = T)




