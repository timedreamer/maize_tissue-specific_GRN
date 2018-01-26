# This script is to plot Supplemental Figure 4, the degree of centrality plots.
# And Figure 6A, the number of connected targets in quantile color scale.

# File in: 1) Top 1 million edges for four tissues:
#             link_list_four_tissue_1million.RData
#
# File out: 1) Plot Supplemental Fig 4: TF_degree_centrality.pdf
#           2) Plot Figure 6A: TF_degree_centrality_quantile_heatmap.pdf
#
################################################################################

# Load data and libraries.
setwd("~/projects/NTWK/tissue_ntwk/")
load("data/right_seed/link_list_four_tissue_1million.RData")
library(igraph)
library(tidyverse)
library(gridExtra)
library(pheatmap)
library(RColorBrewer)


################################################################################
## Plot FOUR COLUMN CHARTS.
# Supplemental Figure 4
################################################################################

# calclate degree for four tissues

cal_deg <- function(link_list){
  link_list %>% 
    select(-weight) %>% 
    group_by(regulatory.gene) %>% 
    summarise(n = n()) %>% 
    arrange(desc(n))
}

deg_sam <- cal_deg(ll_sam)
deg_seed <- cal_deg(ll_seed)
deg_leaf <- cal_deg(ll_leaf)
deg_root <- cal_deg(ll_root)

# I choose TF connected with more than 2000 genes as key TFs.

key_n_leaf <- dim(filter(deg_leaf, n > 2000))[1] # 110
key_n_root <- dim(filter(deg_root, n > 2000))[1] # 53
key_n_sam <- dim(filter(deg_sam, n > 2000))[1] # 88
key_n_seed <- dim(filter(deg_seed, n > 2000))[1] # 56


# plot barplot at four tissues
p_leaf <- ggplot(data = deg_leaf) +
  geom_col(mapping = aes(x = reorder(regulatory.gene, -n), y = n)) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.title = element_text(hjust = 0.5)) +
  labs(x = "TFs in Leaf", y = "Number of targets", 
       title = "TF Degree Centrality in Leaf") +
  geom_segment(aes(x = key_n_leaf, xend = key_n_leaf, y = 0, yend = 5200),
               color = "red", linetype = "dashed",size = 0.3)

p_root <- ggplot(data = deg_root) +
  geom_col(mapping = aes(x = reorder(regulatory.gene, -n), y = n)) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.title = element_text(hjust = 0.5)) +
  labs(x = "TFs in Root", y = "Number of targets", 
       title = "TF Degree Centrality in Root") +
  geom_segment(aes(x = key_n_root, xend = key_n_root, y = 0, yend = 3000),
               color = "red", linetype = "dashed",size = 0.3)

p_sam <- ggplot(data = deg_sam) +
  geom_col(mapping = aes(x = reorder(regulatory.gene, -n), y = n)) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.title = element_text(hjust = 0.5)) +
  labs(x = "TFs in SAM", y = "Number of targets",
       title = "TF Degree Centrality in SAM") +
  geom_segment(aes(x = key_n_sam, xend = key_n_sam, y = 0, yend = 5000),
               color = "red", linetype = "dashed",size = 0.3)

p_seed <- ggplot(data = deg_seed) +
  geom_col(mapping = aes(x = reorder(regulatory.gene, -n), y = n)) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.title = element_text(hjust = 0.5)) +
  labs(x = "TFs in Seed", y = " Number of targets", 
       title = "TF Degree Centrality in Seed") +
  geom_segment(aes(x = key_n_seed, xend = key_n_seed, y = 0, yend = 4000),
               color = "red", linetype = "dashed",size = 0.3)



# Put them together
p_total <- grid.arrange(p_leaf, p_root, p_sam, p_seed, ncol = 2)

# save plot as 8in*8in
ggsave("results/TF_degree_centrality.pdf", 
       plot = p_total, width = 8, height = 8, units = "in")


################################################################################
# PLOT QUANTILE HEATMAP
# Figure 6A
################################################################################

# Because the highest number for TF connection is about 4000-5000, but most of 
# connections are very low. So if I used even color break, most regions will be 
# blue and I do not want to use `scale row`.
# So I tried a script found online using quantile break. Works great.

# define a function to calculate quantile breaks.
quantile_breaks <- function(xs, n = 10) {
  breaks <- quantile(xs, probs = seq(0, 1, length.out = n))
  breaks[!duplicated(breaks)]
}


# Prepare the count number matrix for four tissues.
deg_four <- full_join(deg_leaf, deg_root, by = "regulatory.gene") %>% 
  full_join(., deg_sam, by = "regulatory.gene") %>% 
  full_join(., deg_seed, by = "regulatory.gene") %>% 
  rename(., leaf = n.x, root = n.y, SAM = n.x.x, seed = n.y.y) %>% 
  replace_na(., 
             replace = list(leaf = 0, root = 0, SAM = 0, seed = 0)) %>% 
  select(-regulatory.gene)

# Plot the heatmap 
mat_breaks <- quantile_breaks(as.matrix(deg_four), n = 11)

pdf("results/TF_degree_centrality_quantile_heatmap.pdf",width = 4,height = 6)
pheatmap(
  mat               = as.matrix(deg_four),
  color             = colorRampPalette(rev(brewer.pal(n = 10, name =
                                                        "RdYlBu")))(10),
  # color = inferno(length(mat_breaks) - 1)
  breaks            = mat_breaks,
  border_color      = NA,
  show_colnames     = TRUE,
  show_rownames     = FALSE,
  drop_levels       = TRUE,
  fontsize          = 14,
  fontsize_col      = 20,
  main              = "Quantile Color Scale"
)
dev.off()

