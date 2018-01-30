# This script is to use g:profile to computer GO enrichment for KN1, FEA4 and O2
# targets at each tissue.
# Have to use an archive version of g:profile (eg31), because I use AGPv3 geneID.
#
# File in: 1) Total four million edges from four tissues:
#             link_list_allFour_1million.RDS
# File out: 1) combined GO enrichment result: gene_ontology_three_tfs_all.txt
#
################################################################################

# load libraries and data.
## setwd("~/projects/NTWK/tissue_ntwk/")
ll_four <- readRDS(file = "data/right_seed/link_list_allFour_1million.RDS")
library(tidyverse)
library(gProfileR)


################################################################################
# 1. PREPARE GRN LIST FROM TF AND TISSUE
################################################################################
# KN1 
gene <- "GRMZM2G017087"
kn1_sam <- ll_four %>% 
  filter(tissue == "sam" & regulatory.gene == gene)

kn1_seed <- ll_four %>% 
  filter(tissue == "seed" & regulatory.gene == gene)

# FEA4
gene <- "GRMZM2G133331"
fea4_leaf <- ll_four %>% 
  filter(tissue == "leaf" & regulatory.gene == gene)

fea4_root <- ll_four %>% 
  filter(tissue == "root" & regulatory.gene == gene)

fea4_sam <- ll_four %>% 
  filter(tissue == "sam" & regulatory.gene == gene)

fea4_seed <- ll_four %>% 
  filter(tissue == "seed" & regulatory.gene == gene)

# O2 
gene <- "GRMZM2G015534"
o2_seed <- ll_four %>% 
  filter(tissue == "seed" & regulatory.gene == gene)


################################################################################
## 2. USE gPROFILER TO COMPUTE GENE ONTOLOGY
################################################################################

# Need to use archive for AGPv3
set_base_url("https://biit.cs.ut.ee/gprofiler_archive2/r1622_e84_eg31/web/")
  
# KN1
go_kn1_sam <- gprofiler(kn1_sam$target.gene, organism = "zmays",
                        max_p_value = 0.05, correction_method = "gSCS", 
                        hier_filtering = "none") %>% 
  mutate(gene = "KN1", tissue = "SAM")

go_kn1_seed <- gprofiler(kn1_seed$target.gene, organism = "zmays",
                        max_p_value = 0.05, correction_method = "gSCS", 
                        hier_filtering = "none") %>% 
  mutate(gene = "KN1", tissue = "seed")

# FEA4
go_fea4_leaf <- gprofiler(fea4_leaf$target.gene, organism = "zmays",
                        max_p_value = 0.05, correction_method = "gSCS", 
                        hier_filtering = "none") %>% 
  mutate(gene = "FEA4", tissue = "leaf")

go_fea4_root <- gprofiler(fea4_root$target.gene, organism = "zmays",
                         max_p_value = 0.05, correction_method = "gSCS", 
                         hier_filtering = "none") %>% 
  mutate(gene = "FEA4", tissue = "root")

go_fea4_sam <- gprofiler(fea4_sam$target.gene, organism = "zmays",
                        max_p_value = 0.05, correction_method = "gSCS", 
                        hier_filtering = "none") %>% 
  mutate(gene = "FEA4", tissue = "SAM")

go_fea4_seed <- gprofiler(fea4_seed$target.gene, organism = "zmays",
                         max_p_value = 0.05, correction_method = "gSCS", 
                         hier_filtering = "none") %>% 
  mutate(gene = "FEA4", tissue = "seed")

# O2
go_o2_seed <- gprofiler(o2_seed$target.gene, organism = "zmays",
                          max_p_value = 0.05, correction_method = "gSCS", 
                          hier_filtering = "none") %>% 
  mutate(gene = "O2", tissue = "seed")

# Combine all results into one tibble
go_three_tf <- bind_rows(go_kn1_sam, go_kn1_seed,
               go_fea4_leaf, go_fea4_root, go_fea4_sam, go_fea4_seed,
               go_o2_seed)

# Save result into table.
write.table(go_three_tf, file = "results/gene_ontology_three_tfs_all.txt",
            quote = F, row.names = F, sep = "\t")

