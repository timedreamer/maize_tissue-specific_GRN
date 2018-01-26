# This script is to plot number of libraries in each tissue.
# Used for FIgure 1A.
# 
# Files in : 1) tissue library info:
#                  Four_tissue_specific_libraries_correct_seed.txt
# Files out: 2) a plot: tissue_lib_count.pdf.


################################################################################
# load libraries and setwd.
setwd("~/projects/NTWK/tissue_ntwk/")
library(ggplot2)
library(tidyverse)

# Read tissue with libraries.
tissue_lib <- read_tsv(file = "data/Four_tissue_specific_libraries_correct_seed.txt",
                       col_names = T)

# count number of libs each tissue

lib_leaf <- sum(!is.na(tissue_lib$leaf))
lib_root <- sum(!is.na(tissue_lib$root))
lib_sam <- sum(!is.na(tissue_lib$sam))
lib_seed <- sum(!is.na(tissue_lib$seed))

# prepare a tibble

lib_tissue <- tibble(
  Tissue = c("leaf", "root", "sam", "seed"),
  Count = c(lib_leaf, lib_root, lib_sam, lib_seed)
)

# Plot the barplot using ggplot2
p1 <- ggplot(data = lib_tissue, 
       mapping = aes(x = Tissue, y = Count)) +
  geom_col() +
  labs(y = "Library Count") +
  theme(axis.text = element_text(size = 25),
        axis.title = element_text(size = 25, face = "bold"))

# save the file
ggsave(filename = "results/tissue_lib_count.pdf", plot = p1,
       width = 7, height = 7, units = "in")
