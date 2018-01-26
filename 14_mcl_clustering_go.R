# This script is to caclualte MCL clustering. Then send each cluster to
# goprofile to test GO enrichment.
#
# File in: 1) txt format of top 1 million of four tissue GRNs:
#             ll_leaf.txt, ll_root.txt, ll_sam.txt, ll_seed.txt
# 
# File out: 1) MCL cluster results for each tissue: mcl_1M_leaf_I2.5.txt,
#              mcl_1M_root_I2.5.txt, mcl_1M_sam_I2.5.txt and mcl_1M_seed_I2.5.txt.
#           2) GO enrichment result for each cluster:
#              all saved in "results/cluster_go" directory.
#
################################################################################



# load required packages
library(tidyverse)
library(stringr)
library(gProfileR)

################################################################################
# 1. Calcuate cluters using MCL. Run one time!!
################################################################################

# To call MCL for clustering from R
mcl_path <- "/home5/jhuang/local/bin/mcl"
system("mkdir -p results/cluster_go")

# 604 clusters in leaf.
system(paste(mcl_path, "data/right_seed/ll_leaf.txt --abc",
             "-I 2.5 -te 20 -o results/cluster_go/mcl_1M_leaf_I2.5.txt")) 

# 737 clusters in root.
system(paste(mcl_path, "data/right_seed/ll_root.txt --abc",
             "-I 2.5 -te 20 -o results/cluster_go/mcl_1M_root_I2.5.txt")) 

# 844 clusters in sam.
system(paste(mcl_path, "data/right_seed/ll_sam.txt --abc",
             "-I 2.5 -te 20 -o results/cluster_go/mcl_1M_sam_I2.5.txt")) 

# 399 clusters in seed.
system(paste(mcl_path, "data/right_seed/ll_seed.txt --abc",
             "-I 2.5 -te 20 -o results/cluster_go/mcl_1M_seed_I2.5.txt")) 

rm(mcl_path)

################################################################################
# 2. Read in cluster results
################################################################################

# 2.1 make new directories ready. Run one-time!
system("mkdir results/cluster_go/leaf results/cluster_go/root")
system("mkdir results/cluster_go/sam results/cluster_go/seed")

# 2.2 Set the gprofile server to Ensembl Genome 31 which used AGPv3 geneID.
set_base_url("https://biit.cs.ut.ee/gprofiler_archive2/r1622_e84_eg31/web/")

# 2.3 Read clusters from each as lists. Choose clusters with more than 10 genes.

# 2.3.1 leaf. 232 clusters > 10 genes.
clust_leaf <- read_lines(file = "results/cluster_go/mcl_1M_leaf_I2.5.txt") %>% 
  str_split(., "\t")
clust_leaf <- clust_leaf[lengths(clust_leaf) > 10]

# 2.3.2 root. 278 clusters > 10 genes.
clust_root <- read_lines(file = "results/cluster_go/mcl_1M_root_I2.5.txt") %>% 
  str_split(., "\t")
clust_root <- clust_root[lengths(clust_root) > 10]

# 2.3.3 sam. 268 clusters > 10 genes.
clust_sam <- read_lines(file = "results/cluster_go/mcl_1M_sam_I2.5.txt") %>% 
  str_split(., "\t")
clust_sam <- clust_sam[lengths(clust_sam) > 10]

# 2.3.4 seed. 166 clusters > 10 genes.
clust_seed <- read_lines(file = "results/cluster_go/mcl_1M_seed_I2.5.txt") %>% 
  str_split(., "\t")
clust_seed <- clust_seed[lengths(clust_seed) > 10]

################################################################################
# 3. Query gprofile with those clusters for Gene ontology
################################################################################
# Only kept BP and gene lists less than 1000 genes. 
# Cluster name was based on the lowest p-value.
# Each time query about 30 clusters, otherwise had `Error: Proxy Error`.

# 3.1 Query for leaf
for (i in c(152:length(clust_leaf))) {
  go_enrich <- gprofiler(clust_leaf[[i]], organism = "zmays",
                         max_p_value = 0.05, correction_method = "gSCS", 
                         hier_filtering = "moderate") %>% 
    filter(domain == "BP" & term.size < 1000) %>% 
    arrange(., p.value)
  
  first_term <- go_enrich[1,12]
  write_tsv(go_enrich, path = paste0("results/cluster_go/leaf/", 
                                     "module_", i, "_", first_term, "_BP_GO.txt"))
}

# 3.2 Query for root
for (i in c(265:length(clust_root))) {
  go_enrich <- gprofiler(clust_root[[i]], organism = "zmays",
                         max_p_value = 0.05, correction_method = "gSCS", 
                         hier_filtering = "moderate") %>% 
    filter(domain == "BP" & term.size < 1000) %>% 
    arrange(., p.value)
  
  first_term <- go_enrich[1,12]
  write_tsv(go_enrich, path = paste0("results/cluster_go/root/", 
                                     "module_", i, "_", first_term, "_BP_GO.txt"))
}

# 3.3 Query for sam
for (i in c(151:length(clust_sam))) {
  go_enrich <- gprofiler(clust_sam[[i]], organism = "zmays",
                         max_p_value = 0.05, correction_method = "gSCS", 
                         hier_filtering = "moderate") %>% 
    filter(domain == "BP" & term.size < 1000) %>% 
    arrange(., p.value)
  
  first_term <- go_enrich[1,12]
  write_tsv(go_enrich, path = paste0("results/cluster_go/sam/", 
                                     "module_", i, "_", first_term, "_BP_GO.txt"))
}

# 3.4 Query for seed
for (i in c(160:length(clust_seed))) {
  go_enrich <- gprofiler(clust_seed[[i]], organism = "zmays",
                         max_p_value = 0.05, correction_method = "gSCS", 
                         hier_filtering = "moderate") %>% 
    filter(domain == "BP" & term.size < 1000) %>% 
    arrange(., p.value)
  
  first_term <- go_enrich[1,12]
  write_tsv(go_enrich, path = paste0("results/cluster_go/seed/", 
                                     "module_", i, "_", first_term, "_BP_GO.txt"))
}


# clean
rm(list = ls())
gc()
