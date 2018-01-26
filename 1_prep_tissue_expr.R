# This script is to prepare four expressio tables for four tissues
# including: sam, leaf, root and seed.

# File in: 1. Gene expression matrix (raw count) for 1266 libraries.
#          2. A table indicating libraris from which tissue.

# File out: 1. Four matrix for each tissues gene expression.
#             (expr_leaf, expr_root, expr_sam, expr_seed)
###############################################################################


# load libraries
setwd("~/projects/NTWK/tissue_ntwk/data/")
library(tidyverse)

# read the originial expression table
expr_all <- read_tsv(file = "ALL_FC_noDuplicateLib_biggerThan5Million_70allignmentRate_1266.txt",
                     col_names = T)

tissue_lib <- read_tsv(file = "Four_tissue_specific_libraries_correct_seed.txt",col_names = T)

head(expr_all)
head(tissue_lib)

# get gene names. remove unecessary columns.
gene_name_all <- unname(expr_all$Geneid)
expr_all <- expr_all[,-c(2:6)]

# get the library name for each tissue
tissue_sam <- tissue_lib$sam %>% 
  na.omit %>% 
  as.vector()

tissue_leaf <- tissue_lib$leaf %>% 
  na.omit %>% 
  as.vector()

tissue_root <- tissue_lib$root %>% 
  na.omit %>% 
  as.vector()

tissue_seed <- tissue_lib$seed %>% 
  na.omit %>% 
  as.vector()

# generate tissue specific expression table. These will be used by next script.
expr_sam <- select(expr_all, one_of(tissue_sam))
expr_leaf <- select(expr_all, one_of(tissue_leaf))
expr_root <- select(expr_all, one_of(tissue_root))
expr_seed <- select(expr_all, one_of(tissue_seed))


# remove all files except tissue expression matrix and gene name.
# expr_leaf, expr_root, expr_sam, expr_seed and
# gene_name_all were kept.
rm(expr_all,tissue_lib,tissue_leaf, tissue_root,
   tissue_sam, tissue_seed)
