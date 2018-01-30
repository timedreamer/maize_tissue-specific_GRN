# This script is to use positive dataset: KN1, O2 and FEA4 to calculate 
# AUROC and AUPR values for each tissue network, as well as random networks.
#
# File in: 1) Four weigh matrix from GENIE3: wm_leaf_tf.RData, wm_root_tf.RData
#             wm_sam_tf.RData and wm_seed_tf.RData
#          2) KN1, FEA4 and O2 binding targets: KN1_chip_binding_10k
#             FEA4_binding_10k.txt and o2_bind_gene.txt.
# File out: 1) AUROC and AUPR numbers were recored while running script.
#
##############################################################################



# setwd and load libraries
setwd("~/projects/NTWK/tissue_ntwk/")
library(tidyverse)
library(PRROC)

###############################################################################
## THIS PART ONLY NEED TO RUN ONCE.
###############################################################################
# load Weigh Matrix calculated by GENIE3. 
# Then only save four TFs(KN1, FEA4, P1 and O2) weigh matrix in each tissue.
load("data/right_seed/wm_leaf_tf.RData")
load("data/right_seed/wm_root_tf.RData")
load("data/right_seed/wm_sam_tf.RData")
load("data/right_seed/wm_seed_tf.RData")


wm_leaf_4tf <- wm_leaf[c("GRMZM2G133331", "GRMZM5G895313"),]
wm_root_4tf <- wm_root[c("GRMZM2G133331", "GRMZM5G895313"),]
wm_sam_4tf <- wm_sam[c("GRMZM2G017087", "GRMZM2G133331"),]
wm_seed_4tf <- wm_seed[c("GRMZM2G017087", "GRMZM2G133331", "GRMZM2G015534"),]

# save wm_tissue_4tf.RData.
save(wm_leaf_4tf, wm_root_4tf, wm_sam_4tf, wm_seed_4tf,
     file = "data/right_seed/wm_tissue_4tf.RData")

rm(wm_leaf, wm_root, wm_sam, wm_seed)


###############################################################################
## LOAD DATA AND RUN FROM HERE
################################################################################

# Load data and get four TF's target.
load("data/right_seed/wm_tissue_4tf.RData")

kn1_tgt <- read_tsv("./data/KN1_chip_binding_10k.txt", col_names = "geneid")
kn1_tgt <- as.vector(kn1_tgt$geneid)

fea4_tgt <- read_tsv("./data/FEA4_binding_10k.txt", col_names = "geneid")
fea4_tgt <- as.vector(fea4_tgt$geneid)

p1_tgt <- read_tsv("./data/P1_chipped.txt", col_names = "geneid")
p1_tgt <- as.vector(p1_tgt$geneid)

o2_tgt <- read_tsv("./data/o2_bind_gene.txt", col_names = "geneid")
o2_tgt <- as.vector(o2_tgt$geneid)




## TRUE BINDINGS. These are targets with changes in expression level (RNA-seq)
kn1_tgt <- read_tsv("data/kn1_true.txt", col_names = "geneid")
kn1_tgt <- as.vector(kn1_tgt$geneid)

fea4_tgt <- read_tsv("data/FEA4_true_binding.txt", col_names = "geneid")
fea4_tgt <- as.vector(fea4_tgt$geneid)

kn1_tgt <- read_tsv("data/kn1_true_withleaf.txt", col_names = "geneid")
kn1_tgt <- as.vector(kn1_tgt$geneid)

## Define a function called get_topEdges().
## You can get top x edges from the TFs in from four tissue weigh matrix.
## I can change the `wm_leaf_4tf` to `wm_leaf` to get whatever TF.

get_topEdges <- function(gene, tissue, topEdges){
  if (tissue == "leaf") {
    tf_leaf <- wm_leaf_4tf[gene,]
    tf_leaf <- rank(tf_leaf,ties.method = "random")/length(tf_leaf)
    tf_leaf <- tf_leaf[order(tf_leaf, decreasing = T)][1:topEdges]
    return(tf_leaf) 
  } 
  
  else if (tissue == "root") {
    tf_root <- wm_root_4tf[gene,]
    tf_root <- rank(tf_root,ties.method = "random")/length(tf_root)
    tf_root <- tf_root[order(tf_root, decreasing = T)][1:topEdges]
    return(tf_root)
  }
  
  else if (tissue == "sam") {
    tf_sam <- wm_sam_4tf[gene,]
    tf_sam <- rank(tf_sam,ties.method = "random")/length(tf_sam)
    tf_sam <- tf_sam[order(tf_sam, decreasing = T)][1:topEdges]
    return(tf_sam)
  } 
  
  else if (tissue == "seed") {
    tf_seed <- wm_seed_4tf[gene,]
    tf_seed <- rank(tf_seed,ties.method = "random")/length(tf_seed)
    tf_seed <- tf_seed[order(tf_seed, decreasing = T)][1:topEdges]
    return(tf_seed)
  }
}

## Set how many top edges you want. I set to 1000, 4000, and 10,000.
number <- 4000

# Loops four tissues for KN1 (GRMZM2G017087).
query <- "GRMZM2G017087" #KN1

for (i in c("sam", "seed")) {
  tf <- get_topEdges(query, tissue = i,number)
  tf_bind <- names(tf) %in% kn1_tgt
  tf_df <- as.data.frame(cbind(tf, tf_bind))
  roc <- roc.curve(scores.class0 = tf_df$tf, weights.class0 = tf_df$tf_bind)
  print(paste("AUROC for KN1 in", i, "is",
              round(roc$auc, digits = 3)))
  pr <- pr.curve(scores.class0 = tf_df$tf, weights.class0 = tf_df$tf_bind)
  print(paste("AUPR for KN1 in", i, "is", 
              round(pr$auc.davis.goadrich, digits = 3)))
}


# For FEA4 (GRMZM2G133331)
query <- "GRMZM2G133331" # FEA4

for (i in c("leaf", "root", "sam", "seed")) {
  tf <- get_topEdges(query, tissue = i,number)
  tf_bind <- names(tf) %in% fea4_tgt
  tf_df <- as.data.frame(cbind(tf, tf_bind))
  roc <- roc.curve(scores.class0 = tf_df$tf, weights.class0 = tf_df$tf_bind)
  print(paste("AUROC for FEA4 in", i, "is",
              round(roc$auc, digits = 3)))
  pr <- pr.curve(scores.class0 = tf_df$tf, weights.class0 = tf_df$tf_bind)
  print(paste("AUPR for FEA4 in", i, "is", 
              round(pr$auc.davis.goadrich, digits = 3)))
}



# For O2 (GRMZM2G015534)
query <- "GRMZM2G015534" # O2

for (i in c("seed")) {
  tf <- get_topEdges(query, tissue = i,number)
  tf_bind <- names(tf) %in% o2_tgt
  tf_df <- as.data.frame(cbind(tf, tf_bind))
  roc <- roc.curve(scores.class0 = tf_df$tf, weights.class0 = tf_df$tf_bind)
  print(paste("AUROC for O2 in", i, "is",
              round(roc$auc, digits = 3)))
  pr <- pr.curve(scores.class0 = tf_df$tf, weights.class0 = tf_df$tf_bind)
  print(paste("AUPR for O2 in", i, "is", 
              round(pr$auc.davis.goadrich, digits = 3)))
}

# Then I basically saved those numbers reported in Excel 
# (ChIP_AUC_tissue.xlsx in dropbox).



################################################################################
## Plot AUROC for KN1 and FEA4
################################################################################

# FOR KN1
i <- "seed"
kn1 <- get_topEdges(query, tissue = i,number)
kn1_bind <- names(kn1) %in% kn1_tgt
kn1_df <- as.data.frame(cbind(kn1, kn1_bind))
pred <- prediction(kn1_df$kn1, kn1_df$kn1_bind)
kn1_auc <- performance(pred,"tpr", "fpr")
plot(kn1_auc, col = "#E69F00")
plot(kn1_auc, add = T, col = "#56B4E9")
plot(kn1_auc, add = T, col = "#D55E00")
plot(kn1_auc, add = T, col = "#CC79A7")
abline(0,1)
title(main = "KN1 AUROC at Four tissu")

################################################################################
## CALCULATE AUROC AND AUPR FROM RANDOM GENES
################################################################################

# only need "gene_name_all" to selecet genes
load("data/cpm1_more_than_ten_percent/cpm_tissue_10%.RData")

## Random ROC and AUPR for KN1
kn1 <- c(10000:1)/10000
rand_roc_result <- c()
rand_aupr_result <- c()
for (i in c(1:10000)) {
  set.seed(i)
  names(kn1) <- sample(gene_name_all, size = 10000)
  kn1_bind <- names(kn1) %in% kn1_tgt
  kn1_df <- as.data.frame(cbind(kn1, kn1_bind))
  roc <- roc.curve(scores.class0 = kn1_df$kn1, 
                   weights.class0 = kn1_df$kn1_bind)
  rand_roc_result <- append(rand_roc_result, roc$auc)
  pr <- pr.curve(scores.class0 = kn1_df$kn1, 
                 weights.class0 = kn1_df$kn1_bind)
  rand_aupr_result <- append(rand_aupr_result, pr$auc.davis.goadrich)
}

mean(rand_roc_result)
mean(rand_aupr_result)

## Random ROC and AUPR for FEA4
fea4 <- c(10000:1)/10000
rand_roc_result <- c()
rand_aupr_result <- c()
for (i in c(1:10000)) {
  set.seed(i)
  names(fea4) <- sample(gene_name_all, size = 10000)
  fea4_bind <- names(fea4) %in% fea4_tgt
  fea4_df <- as.data.frame(cbind(fea4, fea4_bind))
  roc <- roc.curve(scores.class0 = fea4_df$fea4, 
                   weights.class0 = fea4_df$fea4_bind)
  rand_roc_result <- append(rand_roc_result, roc$auc)
  pr <- pr.curve(scores.class0 = fea4_df$fea4, 
                 weights.class0 = fea4_df$fea4_bind)
  rand_aupr_result <- append(rand_aupr_result, pr$auc.davis.goadrich)
}

mean(rand_roc_result)
mean(rand_aupr_result)

## Random ROC and AUPR for O2
o2 <- c(10000:1)/10000
rand_roc_result <- c()
rand_aupr_result <- c()
for (i in c(1:10000)) {
  set.seed(i)
  names(o2) <- sample(gene_name_all, size = 10000)
  o2_bind <- names(o2) %in% o2_tgt
  o2_df <- as.data.frame(cbind(o2, o2_bind))
  roc <- roc.curve(scores.class0 = o2_df$o2, 
                   weights.class0 = o2_df$o2_bind)
  rand_roc_result <- append(rand_roc_result, roc$auc)
  pr <- pr.curve(scores.class0 = o2_df$o2, 
                 weights.class0 = o2_df$o2_bind)
  rand_aupr_result <- append(rand_aupr_result, pr$auc.davis.goadrich)
}

mean(rand_roc_result)
mean(rand_aupr_result)