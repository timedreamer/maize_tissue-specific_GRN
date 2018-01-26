# This script is to 1) compare the edges in top 1 million GRNs
# with interactions found in ChIP-Seq data (KN1, FEA4 and O2).
# 2) compare the top 100,000 edges in each tissue GRNs
# with interactions found in ChIP-Seq data (KN1, FEA4 and O2).
# 3) compare the top 4000 targets in each tissue GRNs
# with interactions found in ChIP-Seq data (KN1, FEA4 and O2).
# The result can be used in chi-square test.
# 4) Do the same thing for Briggs RNA and protein GRNs, top 1 million edges
# and top 100,000 edges. 
# 5) All comparisons can be tested by one-tail Fisher's exact test and/or 
# chi-square test.

# File in: 1) top 1 million edges for each tissue GRN:
#             ll_leaf, ll_root, ll_sam and ll_seed
#          2) top 100,000 edges for each tissue GRN, same name.
#          3) ChIP target genes for KN1, FEA4 and O2

# File out: 1) Number of the overlap (record while running script)
#           2) Top 1 million edges from four tissues were combined into one
#              long table (link_list_allFour_1million.RData)
###############################################################################

setwd("~/projects/NTWK/tissue_ntwk/")
load("data/right_seed/link_list_four_tissue_1million.RData")
load("data/right_seed/link_list_four_tissue_100thousand.RData")
library(tidyverse)


# Read chip-seq target list. ChIP only,
kn1 <- read_tsv("./data/KN1_chip_binding_10k.txt", col_names = "geneid")
kn1 <- as.vector(kn1$geneid)
fea4 <- read_tsv("data/FEA4_binding_10k.txt", col_names = "geneid")
fea4 <- as.vector(fea4$geneid)
o2 <- read_tsv("./data/o2_bind_gene.txt", col_names = "geneid")
o2 <- as.vector(o2$geneid)
# all geneids are in `length(kn1 %in% gene_name_all)`

# Try list with expression changes (true binding).
# O2 is from 
kn1 <- read_tsv("./data/kn1_true.txt", col_names = "geneid")
kn1 <- as.vector(kn1$geneid)

kn1 <- read_tsv("data/kn1_true_withleaf.txt", col_names = "geneid")
kn1 <- as.vector(kn1$geneid)

fea4 <- read_tsv("data/FEA4_true_binding.txt", col_names = "geneid")
fea4 <- as.vector(fea4$geneid)
o2 <- read_tsv("data/o2_true_binding.txt", col_names = "geneid")
o2 <- as.vector(o2$geneid)

# combine four link_list into one.
ll_leaf <- mutate(ll_leaf, tissue = "leaf")
ll_root <- mutate(ll_root, tissue = "root")
ll_sam <- mutate(ll_sam, tissue = "sam")
ll_seed <- mutate(ll_seed, tissue = "seed")

ll_four <- bind_rows(ll_leaf, ll_root, ll_sam, ll_seed)

saveRDS(object = ll_four, 
        file = "data/right_seed/link_list_allFour_10million.RData")

# make sure the number is correct
ll_four %>% 
  group_by(tissue) %>% 
  summarise(n = n())


# Query gene (KN1, P1, O2 or FEA4)
# KN1: GRMZM2G017087; O2: GRMZM2G015534; P1: GRMZM2G084799; FEA4: GRMZM2G133331

# write a function to calcule edges numbers
edges_in_chip <- function(gene, chip){
  edg <- ll_four %>% 
    filter(regulatory.gene == gene) %>% 
    group_by(tissue) %>% 
    summarise(n = n())
  
  edg_chip <- ll_four %>% 
    filter(regulatory.gene == gene) %>% 
    filter(target.gene %in% chip) %>% 
    group_by(tissue) %>% 
    summarise(n = n())
  
  left_join(edg, edg_chip, by = "tissue")
}

# the first number is how many edges in our top 1 million GRN.
# the second number is how many are overlap with ChIP-Seq targets.

edges_in_chip("GRMZM2G017087", kn1) # KN1: GRMZM2G017087 
edges_in_chip("GRMZM2G084799", p1)  # P1: GRMZM2G084799
edges_in_chip("GRMZM2G015534", o2) # O2: GRMZM2G015534
# FEA4: GRMZM2G133331 number is similar using KN1
edges_in_chip("GRMZM2G133331", fea4) 

# I can use chi-square test to test enrichment. 
# columns are "Regulated", "not-Regulated"; rows are "binding", "not binding"
#######    Regulated  not-Regulated
# Bind        5           4060
# Not Bind    67         35347

test <- matrix(c(51,743,2357,36328), nrow = 2)
test
set.seed(123)
fisher.test(test, simulate.p.value = T, alternative = "greater")
chisq.test(test,simulate.p.value = T)
#phyper(191,758,38721,4274, lower.tail = T)

plantReg <- read_tsv("regulation_merged_Zma.txt",
                     col_names = c("regulator", "direction", 
                                   "target", "species", "type"))

head(plantReg)

edg_chip <- ll_four %>% 
  filter(regulatory.gene == gene) %>% 
  filter(target.gene %in% chip) %>% 
  group_by(tissue) %>% 
  summarise(n = n())


plantReg %>% 
  filter(regulator == "GRMZM2G017087") %>% 
  filter(target %in% kn1) %>%
  summarise(n = n())

################################################################################
## GET TOP ~4000 GENES AND FIND OVERLAP WITH CHIP TARGET
################################################################################

load("data/right_seed/wm_tissue_4tf.RData")

# KN1 top 4274
kn1_sam_4274 <- names(sort(wm_sam_4tf["GRMZM2G017087",], 
                           decreasing = T)[1:4274])

kn1_seed_4274 <- names(sort(wm_seed_4tf["GRMZM2G017087",], 
                           decreasing = T)[1:4274])

print(paste("Top 4274 predicted KN1 target from SAM GRN has", 
            length(intersect(kn1_sam_4274, kn1)), 
      "overlap with ChIP binding targets") )

print(paste("Top 4274 predicted KN1 target from seed GRN has", 
            length(intersect(kn1_seed_4274, kn1)), 
            "overlap with ChIP binding targets") )

# FEA4 top 4065
fea4_leaf_4065 <- names(sort(wm_leaf_4tf["GRMZM2G133331",], 
                           decreasing = T)[1:4065])

fea4_root_4065 <- names(sort(wm_root_4tf["GRMZM2G133331",], 
                            decreasing = T)[1:4065])

fea4_sam_4065 <- names(sort(wm_sam_4tf["GRMZM2G133331",], 
                           decreasing = T)[1:4065])

fea4_seed_4065 <- names(sort(wm_seed_4tf["GRMZM2G133331",], 
                            decreasing = T)[1:4065])

print(paste("Top 4065 predicted FEA4 target from leaf GRN has", 
            length(intersect(fea4_leaf_4065, fea4)), 
            "overlap with ChIP binding targets") )

print(paste("Top 4065 predicted FEA4 target from root GRN has", 
            length(intersect(fea4_root_4065, fea4)), 
            "overlap with ChIP binding targets") )

print(paste("Top 4065 predicted FEA4 target from SAM GRN has", 
            length(intersect(fea4_sam_4065, fea4)), 
            "overlap with ChIP binding targets") )

print(paste("Top 4065 predicted FEA4 target from seed GRN has", 
            length(intersect(fea4_seed_4065, fea4)), 
            "overlap with ChIP binding targets") )


# O2 top 2408
o2_seed_2408 <- names(sort(wm_seed_4tf["GRMZM2G015534",], 
                             decreasing = T)[1:2408])

print(paste("Top 2408 predicted O2 target from seed GRN has", 
            length(intersect(o2_seed_2408, o2)), 
            "overlap with ChIP binding targets") )

################################################################################
# random overlap
################################################################################

# KN1
kn1_random_4274 <- c()

for (i in c(1:10000)) {
  set.seed(i)
  rd_sample <- sample(x = names(wm_seed_4tf[1,]), replace = F, size = 4274)
  rd_number <- length(intersect(rd_sample, kn1))
  kn1_random_4274 <- append(kn1_random_4274, rd_number)
} 

print(mean(kn1_random_4274))
hist(kn1_random_4274)

# FEA4
fea4_random_4065 <- c()

for (i in c(1:10000)) {
  set.seed(i)
  rd_sample <- sample(x = names(wm_seed_4tf[1,]), replace = F, size = 4065)
  rd_number <- length(intersect(rd_sample, fea4))
  fea4_random_4065 <- append(fea4_random_4065, rd_number)
} 

print(mean(fea4_random_4065))
hist(fea4_random_4065)


################################################################################
# BRIGGS NETWORK OVERLAP CHIPSEQ
################################################################################
# Briggs networks were downloaded from https://goo.gl/9YYgBX.

# 1. Read briggs networks. 1 million networks
brig_protein <- read_tsv(file = "data/briggs/briggs_protein_only.txt")
brig_rna <- read_tsv(file = "data/briggs/briggs_rna_only.txt")
brig_three <- read_tsv(file = "data/briggs/briggs_protein+rna+phospho.txt")

# 2. Calculate overlap of Briggs network with ChIP bindings.

# 2.1 KN1
gene <- "GRMZM2G017087" #KN1

# How many edges in each network?
tf_brig_protein <- brig_protein %>% 
  filter(regulator == gene)

tf_brig_rna <- brig_rna %>% 
  filter(regulator == gene)

tf_brig_three <- brig_three %>% 
  filter(regulator == gene)

# How many genes have ChIP binding?
tf_brig_protein_true <- brig_protein %>% 
  filter(regulator == gene) %>%
  filter(target %in% kn1)

tf_brig_rna_true <- brig_rna %>% 
  filter(regulator == gene) %>% 
  filter(target %in% kn1)

tf_brig_three_true <- brig_three %>% 
  filter(regulator == gene) %>% 
  filter(target %in% kn1)

# Print the result.
print("Below is the result for KN1.")
print(paste("In Briggs protein only network,",
            dim(tf_brig_protein_true)[1], "out of",
            dim(tf_brig_protein)[1], "confirmed by ChIP"))

print(paste("In Briggs RNA only network,",
            dim(tf_brig_rna_true)[1], "out of",
            dim(tf_brig_rna)[1], "confirmed by ChIP"))

print(paste("In Briggs three network,",
            dim(tf_brig_three_true)[1], "out of",
            dim(tf_brig_three)[1], "confirmed by ChIP"))

# 2.2 FEA4
gene <- "GRMZM2G133331" #FEA4

# How many edges in each network?
tf_brig_protein <- brig_protein %>% 
  filter(regulator == gene)

tf_brig_rna <- brig_rna %>% 
  filter(regulator == gene)

tf_brig_three <- brig_three %>% 
  filter(regulator == gene)

# How many genes have ChIP binding?
tf_brig_protein_true <- brig_protein %>% 
  filter(regulator == gene) %>%
  filter(target %in% fea4)

tf_brig_rna_true <- brig_rna %>% 
  filter(regulator == gene) %>% 
  filter(target %in% fea4)

tf_brig_three_true <- brig_three %>% 
  filter(regulator == gene) %>% 
  filter(target %in% fea4)

# Print the result.
print("Below is the result for FEA4.")
print(paste("In Briggs protein only network,",
            dim(tf_brig_protein_true)[1], "out of",
            dim(tf_brig_protein)[1], "confirmed by ChIP"))

print(paste("In Briggs RNA only network,",
            dim(tf_brig_rna_true)[1], "out of",
            dim(tf_brig_rna)[1], "confirmed by ChIP"))

print(paste("In Briggs three network,",
            dim(tf_brig_three_true)[1], "out of",
            dim(tf_brig_three)[1], "confirmed by ChIP"))


# 2.3 O2
gene <- "GRMZM2G015534" # O2

# How many edges in each network?
tf_brig_protein <- brig_protein %>% 
  filter(regulator == gene)

tf_brig_rna <- brig_rna %>% 
  filter(regulator == gene)

tf_brig_three <- brig_three %>% 
  filter(regulator == gene)

# How many genes have ChIP binding?
tf_brig_protein_true <- brig_protein %>% 
  filter(regulator == gene) %>%
  filter(target %in% o2)

tf_brig_rna_true <- brig_rna %>% 
  filter(regulator == gene) %>% 
  filter(target %in% o2)

tf_brig_three_true <- brig_three %>% 
  filter(regulator == gene) %>% 
  filter(target %in% o2)

# Print the result.
print("Below is the result for O2.")
print(paste("In Briggs protein only network,",
            dim(tf_brig_protein_true)[1], "out of",
            dim(tf_brig_protein)[1], "confirmed by ChIP"))

print(paste("In Briggs RNA only network,",
            dim(tf_brig_rna_true)[1], "out of",
            dim(tf_brig_rna)[1], "confirmed by ChIP"))

print(paste("In Briggs three network,",
            dim(tf_brig_three_true)[1], "out of",
            dim(tf_brig_three)[1], "confirmed by ChIP"))


################################################################################

# Also, choose the top 100,000 from Briggs networks (already sorted).
# Do and compare it again.

brig_protein <- brig_protein[1:100000,]
brig_rna <- brig_rna[1:100000,]
brig_three <- brig_three[1:100000,]

# Then re-run the previous code

