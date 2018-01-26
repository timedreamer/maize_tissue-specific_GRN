################################################################################
# DISCARD THIS SCRIPT
# USE TARGET LISTS PLUS http://bioinformatics.psb.ugent.be/webtools/Venn/!!
################################################################################

# This is to caculate edges differences between different tissue GRNs.
# Numbers reported were used to draw Venn diagrams.

setwd("~/projects/NTWK/tissue_ntwk/")
load("data/right_seed/link_list_four_tissue_1million.RData")
library(tidyverse)

# combine four link_list into one.
ll_leaf <- mutate(ll_leaf, tissue = "leaf")
ll_root <- mutate(ll_root, tissue = "root")
ll_sam <- mutate(ll_sam, tissue = "sam")
ll_seed <- mutate(ll_seed, tissue = "seed")

ll_four <- bind_rows(ll_leaf, ll_root, ll_sam, ll_seed)

saveRDS(object = ll_four, 
        file = "data/right_seed/link_list_allFour_1million.RDS")

# KN1 GRMZM2G017087 (SAM and Seed)
gene <- "GRMZM2G017087"

kn1_sam <- ll_four %>% 
  filter(regulatory.gene == gene, tissue == "sam") %>% 
  mutate(name = paste(regulatory.gene, target.gene, sep = "+")) %>%
  arrange(name) %>% 
  select(name)

kn1_seed <- ll_four %>% 
  filter(regulatory.gene == gene, tissue == "seed") %>% 
  mutate(name = paste(regulatory.gene, target.gene, sep = "+")) %>%
  arrange(name) %>% 
  select(name)

length(intersect(kn1_sam$name, kn1_seed$name))

# FEA4 GRMZM2G133331
gene <- "GRMZM2G133331"

fea4_leaf <- ll_four %>% 
  filter(regulatory.gene == gene, tissue == "leaf") %>% 
  mutate(name = paste(regulatory.gene, target.gene, sep = "+")) %>%
  arrange(name) %>% 
  select(name)

fea4_root <- ll_four %>% 
  filter(regulatory.gene == gene, tissue == "root") %>% 
  mutate(name = paste(regulatory.gene, target.gene, sep = "+")) %>%
  arrange(name) %>% 
  select(name)

fea4_sam <- ll_four %>% 
  filter(regulatory.gene == gene, tissue == "sam") %>% 
  mutate(name = paste(regulatory.gene, target.gene, sep = "+")) %>%
  arrange(name) %>% 
  select(name)

fea4_seed <- ll_four %>% 
  filter(regulatory.gene == gene, tissue == "seed") %>% 
  mutate(name = paste(regulatory.gene, target.gene, sep = "+")) %>%
  arrange(name) %>% 
  select(name)

# Two overlap
length(intersect(fea4_leaf$name, fea4_root$name))
length(intersect(fea4_leaf$name, fea4_sam$name))
length(intersect(fea4_leaf$name, fea4_seed$name))
length(intersect(fea4_root$name, fea4_sam$name))
length(intersect(fea4_root$name, fea4_seed$name))
length(intersect(fea4_sam$name, fea4_seed$name))

# Three overlap
length(Reduce(intersect, list(fea4_leaf$name, 
                              fea4_root$name, fea4_sam$name)))
length(Reduce(intersect, list(fea4_leaf$name, 
                              fea4_root$name, fea4_seed$name)))
length(Reduce(intersect, list(fea4_leaf$name, 
                              fea4_seed$name, fea4_sam$name)))
length(Reduce(intersect, list(fea4_root$name, 
                              fea4_seed$name, fea4_sam$name)))

# Four overlap zero

## Clean everything
rm(list = ls())
gc()