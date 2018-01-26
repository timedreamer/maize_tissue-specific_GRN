# This script is to plot increased positive percentage when use smaller network.
# So from top 10 million edges --> top 1 million edges --> top 100k edges,
# we can observe a increased positive rate.
# Used for Figure 4.
#
# File in: 1) positive_percentage_ChIP_right_seed_10million.txt.
#
# File out: 1) Plot figure 4: enrich_TF_postive_rate.svg
#
################################################################################

# Load libary.
library(tidyverse)


# 1. Read data.
# Large network: top 10 million edges; medium: top 1 million, small: top 100,000
enrich <- read_tsv("results/positive_percentage_ChIP_right_seed_10million.txt", 
                   col_names = T)
head(enrich)
# Combine TF with tissue 
enrich <- enrich %>% 
  unite("TF_tissue", TF, tissue)

# Change levels of x-axis 
enrich$level <- factor(enrich$level, levels = c("large", "medium", "small"))

# Plot lines with points.
p <- ggplot(data = enrich, aes(y = percentage, x = level,
                               color = TF_tissue, stat = "identity")) +
  geom_point(aes(group = TF_tissue), size = 4) + 
  geom_line(aes(group = TF_tissue), size = 2)

# Change axis size and add title.
p + theme(axis.text = element_text(size = 18), 
          axis.title = element_text(size = 18, face = "bold")) + 
  ggtitle("enrich TF postive rate with network size")

# Save the graph as pdf.
ggsave(filename = "results/enrich_TF_postive_rate.svg", width = 7, height = 6)

# clean
rm(list = ls())
gc()

