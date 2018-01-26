# This script is to model TF connectivity (how many targets each TF) with their
# expression. The question I want to answer is: Do highy expressed TF have more 
# targets?
# The figure was used for Figure 5 and Supplemental Figure 3.
# Also, a table with coefficient of variance (CV) and difference in TF
# connectivity was saved.
#
# File in: 1) 1409 TFs that expressed in all tissues: TF_in_4tissues_1409.txt
#          2) expression table: PlantGenomeS1_fourTissue_only.txt
#          3) top 1 million edges for four tissue GRNs:
#              link_list_four_tissue_1million.RData
#
# File out: 1) Scatterplot for Figure 5 (log2): plot_TF_interactions_scatter.png 
#           2) Scatterplot for Supplemental Figure 3 (no log2):
#                plot_TF_interactions_scatter_nolog2.png
#           3) linear regression models: record as code runs
#           4) A table with CV and difference in centrality: 
#                TF_interaction_diff.txt
#
################################################################################

# load library and setwd
library(tidyverse)
library(gridExtra)
library(broom)
library(stringr)


################################################################################
# 1. Prepare Expression Matrix
################################################################################

# read 1409 TFs that expressed in 4 tissues.
tf_all <- read_tsv("data/right_seed/TF_in_4tissues_1409.txt", 
                   col_names = "geneID")

# read expression table
pg_expr <- read_tsv("data/PlantGenomeS1_fourTissue_only.txt", col_names = T) %>% 
  filter(geneid %in% tf_all$geneID)

# calculate mean value and prepare matrix.
mydata <- cbind(MeanI = rowMeans(pg_expr[2:19],na.rm = T))
mydata <- cbind(mydata, MeanI = rowMeans(pg_expr[20:38],na.rm = T))
mydata <- cbind(mydata, MeanI = rowMeans(pg_expr[39:40],na.rm = T))
mydata <- cbind(mydata, MeanI = rowMeans(pg_expr[41:64],na.rm = T))
colnames(mydata) <- c("root", "leaf", "SAM", "seed")

# remove TFs with sd = 0.
expr_sd <- apply(mydata,1, sd)
expr <- as_tibble(mydata[(expr_sd != 0),]) %>% 
  add_column(geneID = pg_expr$geneid, .before = "root") %>% 
  arrange(geneID)

################################################################################
# 2. Prepare Connectivity data.
################################################################################

load("data/right_seed/link_list_four_tissue_1million.RData")

# combine four link_list into one.
ll_leaf <- mutate(ll_leaf, tissue = "leaf")
ll_root <- mutate(ll_root, tissue = "root")
ll_sam <- mutate(ll_sam, tissue = "sam")
ll_seed <- mutate(ll_seed, tissue = "seed")

ll_four <- bind_rows(ll_leaf, ll_root, ll_sam, ll_seed)

# make sure the number is correct
ll_four %>% 
  group_by(tissue) %>% 
  summarise(n = n())

rm(ll_leaf, ll_root, ll_sam, ll_seed)

# Find TFs that inlcude in all four tissues. 1406 TFs have at least one
# interactions in one tissue.

leaf <- ll_four %>%
  filter(tissue == "leaf") %>% 
  distinct(regulatory.gene)

root <-  ll_four %>%
  filter(tissue == "root") %>% 
  distinct(regulatory.gene)

sam <-  ll_four %>%
  filter(tissue == "sam") %>% 
  distinct(regulatory.gene)

seed <-  ll_four %>%
  filter(tissue == "seed") %>% 
  distinct(regulatory.gene)

four_tf <- Reduce(f = intersect, 
               list(leaf$regulatory.gene, root$regulatory.gene,
                    sam$regulatory.gene, seed$regulatory.gene) )

# Find TFs with number of connections in each tissue

count_leaf <- ll_four %>% 
  filter(tissue == "leaf") %>% 
  count(regulatory.gene) 

count_root <- ll_four %>% 
  filter(tissue == "root") %>% 
  count(regulatory.gene)

count_sam <- ll_four %>% 
  filter(tissue == "sam") %>% 
  count(regulatory.gene)

count_seed <- ll_four %>% 
  filter(tissue == "seed") %>% 
  count(regulatory.gene) 

count_four <- inner_join(count_leaf, count_root, 
                         by = "regulatory.gene", suffix = c("leaf", "root")) %>% 
  inner_join(count_sam, by = "regulatory.gene") %>% 
  inner_join(count_seed, by = "regulatory.gene", suffix = c("sam", "seed"))


# Join expression table with count table
expr_count <- inner_join(expr, count_four, by = c("geneID" = "regulatory.gene"))

# rm everything except expr_count
rm(count_four, count_leaf, count_root, count_sam, count_seed)
rm(expr, leaf, root, sam, seed)
rm(ll_four, mydata, tf_all, expr_sd, four_tf, pg_expr)

################################################################################
# 3. Plot scatterplots for four tissues and save.
################################################################################

# Log2 transformed.
p_leaf <- ggplot(expr_count, aes(x = log2(leaf + 1), y = nleaf)) + 
  geom_point() + geom_smooth(method = 'lm') +
  labs(title = "leaf", x = "log2(gene expression)", 
       y = "Number of Interactions") +
  theme(text = element_text(size = 15), plot.title = element_text(hjust = 0.5))

p_root <- ggplot(expr_count, aes(x = log2(root + 1), y = nroot)) + 
  geom_point() + geom_smooth(method = 'lm') +
  labs(title = "root", x = "log2(gene expression)", 
       y = "Number of Interactions") +
  theme(text = element_text(size = 15), plot.title = element_text(hjust = 0.5))

p_sam <- ggplot(expr_count, aes(x = log2(SAM + 1), y = nsam)) + 
  geom_point() + geom_smooth(method = 'lm') +
  labs(title = "SAM", x = "log2(gene expression)", 
       y = "Number of Interactions") +
  theme(text = element_text(size = 15), plot.title = element_text(hjust = 0.5))

p_seed <- ggplot(expr_count, aes(x = log2(seed + 1), y = nseed)) + 
  geom_point() + geom_smooth(method = 'lm') +
  labs(title = "seed", x = "log2(gene expression)", 
       y = "Number of Interactions") +
  theme(text = element_text(size = 15), plot.title = element_text(hjust = 0.5))

p_all <- grid.arrange(p_leaf, p_root, p_sam, p_seed, nrow = 2)

ggsave(p_all, filename = "results/plot_TF_interactions_scatter.png", dpi = 600,
       width = 8, height = 8)

# no log2 transform.
p_leaf <- ggplot(expr_count, aes(x = leaf, y = nleaf)) + 
  geom_point() + geom_smooth(method = 'lm') +
  labs(title = "leaf", x = "(gene expression)", 
       y = "Number of Interactions") +
  theme(text = element_text(size = 15), plot.title = element_text(hjust = 0.5))

p_root <- ggplot(expr_count, aes(x = root, y = nroot)) + 
  geom_point() + geom_smooth(method = 'lm') +
  labs(title = "root", x = "gene expression", 
       y = "Number of Interactions") +
  theme(text = element_text(size = 15), plot.title = element_text(hjust = 0.5))

p_sam <- ggplot(expr_count, aes(x = SAM, y = nsam)) + 
  geom_point() + geom_smooth(method = 'lm') +
  labs(title = "SAM", x = "gene expression", 
       y = "Number of Interactions") +
  theme(text = element_text(size = 15), plot.title = element_text(hjust = 0.5))

p_seed <- ggplot(expr_count, aes(x = seed, y = nseed)) + 
  geom_point() + geom_smooth(method = 'lm') +
  labs(title = "seed", x = "gene expression", 
       y = "Number of Interactions") +
  theme(text = element_text(size = 15), plot.title = element_text(hjust = 0.5))

p_all <- grid.arrange(p_leaf, p_root, p_sam, p_seed, nrow = 2)

ggsave(p_all, filename = "results/plot_TF_interactions_scatter_nolog2.png",
       dpi = 600, width = 8, height = 8)

################################################################################
# 4. Linear Regression.
################################################################################

# Linear regression for log2 transformed.
lm_leaf <- glance(lm(formula = nleaf ~ log2(leaf + 1), data = expr_count))
lm_root <- glance(lm(formula = nroot ~ log2(root + 1), data = expr_count))
lm_sam <- glance(lm(formula = nsam ~ log2(SAM + 1), data = expr_count))
lm_seed <- glance(lm(formula = nseed ~ log2(seed + 1), data = expr_count))

print("Below is the linear regression result for log2 expression data")
print(paste("R2 for leaf, root, SAM and seed is:", 
            round(lm_leaf$adj.r.squared,3), round(lm_root$adj.r.squared, 3),
            round(lm_sam$adj.r.squared, 3), round(lm_seed$adj.r.squared, 3)))
print(paste("P-value for leaf, root, SAM and seed is:", 
            lm_leaf$p.value, lm_root$p.value, lm_sam$p.value, 
            lm_seed$p.value))

# Linear regression for original expression data (no log2 transformed).
lm_leaf <- glance(lm(formula = nleaf ~ leaf, data = expr_count))
lm_root <- glance(lm(formula = nroot ~ root, data = expr_count))
lm_sam <- glance(lm(formula = nsam ~ SAM, data = expr_count))
lm_seed <- glance(lm(formula = nseed ~ seed, data = expr_count))

print("Below is the linear regression result for original expression data")
print(paste("R2 for leaf, root, SAM and seed is:", 
            round(lm_leaf$adj.r.squared,3), round(lm_root$adj.r.squared, 3),
            round(lm_sam$adj.r.squared, 3), round(lm_seed$adj.r.squared, 3)))
print(paste("P-value for leaf, root, SAM and seed is:", 
            lm_leaf$p.value, lm_root$p.value, lm_sam$p.value, 
            lm_seed$p.value))


################################################################################
# 5. Save a table
# Table contains TF expression, centrality, CV, 
# the tissue with higest centrality and the difference between 
# highest and lowest centrality. 
################################################################################

# calculate diff in TF interactions and write the table.
# diff in degreee centrality
count_diff <- expr_count %>% 
  mutate(ndiff = pmax(nleaf, nroot, nsam, nseed) - 
           pmin(nleaf, nroot, nsam, nseed)) %>% 
  arrange(desc(ndiff))

# find the tissue with max interactions.
max_tissue <- as_tibble(colnames(count_diff[,c(6:9)])[max.col(count_diff[,c(6:9)],
                                           ties.method = "first")]) %>% 
  rename(max_tissue = value)

# delete "n" in the tissue name and change to factors.
count_diff <- bind_cols(count_diff, max_tissue) %>% 
  mutate(max_tissue = str_replace(max_tissue, "^n", "")) %>% 
  mutate(max_tissue = as.factor(max_tissue))

# add CV
CV <- function(x){
  (sd(x)/mean(x))*100
}

cv_table <- as_tibble(apply(count_diff[,c(6:9)], 1, CV)) %>% 
  rename(cv = value)


count_diff <- count_diff %>% 
  bind_cols(., cv_table) %>% 
  arrange(desc(cv)) %>% 
  filter(., ndiff > 500)

summary(count_diff[1:20,]$max_tissue)
summary(count_diff[1:50,]$max_tissue)
summary(count_diff[1:100,]$max_tissue)
# save the result.
write_tsv(count_diff, path = "results/TF_interaction_diff.txt", col_names = T)


# clean
rm(list = ls())
gc()
