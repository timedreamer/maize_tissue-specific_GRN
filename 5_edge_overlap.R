# This script is to calculate edges overlap between
# four tissu GRN and the PlantRegMap. These numbers are used for Venn diagram.
# Save top 1 million edges name from four tissues,
# and all edges name from platRegMap.
# Also, in Part II, it compared tissue GRN with Briggs atlas RNA/protein GRNs. 

# File in: 1. Four top 1 million edges:
#             (ll_leaf, ll_root, ll_sam, ll_seed)
#          2. Network file from PlantRegMap. plantReg
#          3. Three Briggs GRNs (part II)

# File out: 1. Top 1 million edges, each line as string.
#             edge_top1M_five_ntwk.RData: edge_leaf, edge_root, edge_sam,
#             edge_seed, edge_plantReg.
#           2. A table for edges shared by four tissue GRNs: 
#              overlap_four_2679.txt
#           3. Other numbers were recorded after run script.
################################################################################

# Load packages and data. 
setwd("~/projects/NTWK/tissue_ntwk/")
load("data/right_seed/link_list_four_tissue_1million.RData")
library(tidyverse)
library(stringr)


# Read network file from PlantRegMap

plantReg <- read_tsv("regulation_merged_Zma.txt",
                     col_names = c("regulator", "direction", 
                                   "target", "species", "type"))

head(plantReg)

################################################################################
# PART 1. CALCULATE EDGES OVERLAPS
################################################################################
#  Numbers were used to draw Venn diagram.

# combine regulator and target as edge
edge_plantReg <- plantReg %>% 
  mutate(name = paste(regulator, target, sep = "+")) %>% 
  arrange(name) %>% 
  select(name)


# concatenate regulator with target, so each $name is a edge
edge_leaf <- ll_leaf %>% 
  mutate(name = paste(regulatory.gene, target.gene, sep = "+")) %>% 
  arrange(name) %>% 
  select(name) 


edge_root <- ll_root %>% 
  mutate(name = paste(regulatory.gene, target.gene, sep = "+")) %>% 
  arrange(name) %>% 
  select(name)

edge_seed <- ll_seed %>% 
  mutate(name = paste(regulatory.gene, target.gene, sep = "+")) %>% 
  arrange(name) %>% 
  select(name)


edge_sam <- ll_sam %>% 
  mutate(name = paste(regulatory.gene, target.gene, sep = "+")) %>% 
  arrange(name) %>% 
  select(name)

## save these five edges

save(edge_plantReg, edge_leaf, edge_root, edge_seed, edge_sam,
     file = "data/right_seed/edge_top1M_five_ntwk.RData")


# Calculate each categories value to draw venn diagram.
load("data/right_seed/edge_top1M_five_ntwk.RData")

# Two overlap
leaf_root <- intersect(edge_leaf$name, edge_root$name)
leaf_sam <- intersect(edge_leaf$name, edge_sam$name)
leaf_seed <- intersect(edge_leaf$name, edge_seed$name)
root_seed <- intersect(edge_root$name, edge_seed$name)
root_sam <- intersect(edge_root$name, edge_sam$name)
sam_seed <- intersect(edge_sam$name, edge_seed$name)

# Three overlap
leaf_root_sam <- Reduce(intersect, list(edge_leaf$name, 
                                        edge_root$name, edge_sam$name))
leaf_root_seed <- Reduce(intersect, list(edge_leaf$name, 
                                         edge_root$name, edge_seed$name))
leaf_seed_sam <- Reduce(intersect, list(edge_leaf$name, 
                                        edge_seed$name, edge_sam$name))
root_seed_sam <- Reduce(intersect, list(edge_root$name, 
                                        edge_seed$name, edge_sam$name))


# Four overlap
four_overlap <- Reduce(intersect, list(edge_leaf$name, edge_root$name, 
                                       edge_seed$name, edge_sam$name))

write.table(four_overlap, file = "results/overlap_four_2679.txt", 
            quote = F, sep = "\t", row.names = F)


###############################################################################
# PART 2. COMPARE WITH BRIGGS TOP 1MILLION NETWORKS
###############################################################################

# Briggs networks were downloaded from https://goo.gl/9YYgBX.

# 1. Read briggs networks. 1 million networks
brig_protein <- read_tsv(file = "data/briggs/briggs_protein_only.txt")
brig_rna <- read_tsv(file = "data/briggs/briggs_rna_only.txt")
brig_three <- read_tsv(file = "data/briggs/briggs_protein+rna+phospho.txt")

# 2. Get all edges 
edge_brig_protein <- brig_protein %>% 
  mutate(name = paste(regulator, target, sep = "+")) %>% 
  arrange(name) %>% 
  select(name)

edge_brig_rna <- brig_rna %>% 
  mutate(name = paste(regulator, target, sep = "+")) %>% 
  arrange(name) %>% 
  select(name)


edge_brig_three <- brig_three %>% 
  mutate(name = paste(regulator, target, sep = "+")) %>% 
  arrange(name) %>% 
  select(name)


# Overlap Briggs protein only with RNA only.
length(intersect(edge_brig_protein$name, edge_brig_rna$name)) # 48574

# Overlap Briggs with my leaf GRN.
length(intersect(edge_leaf$name, edge_brig_protein$name)) # 18314
length(intersect(edge_leaf$name, edge_brig_rna$name)) # 29040
length(intersect(edge_leaf$name, edge_brig_three$name)) # 19004

# Calculate how many TFs in Brigs are also in my leaf GRN.
brig_protein_leaf <- intersect(ll_leaf$regulatory.gene, brig_protein$regulator)
brig_rna_leaf <- intersect(ll_leaf$regulatory.gene, brig_rna$regulator)
brig_three_leaf <- intersect(ll_leaf$regulatory.gene, brig_three$regulator)

# Clean everything
rm(list = ls())
gc()
