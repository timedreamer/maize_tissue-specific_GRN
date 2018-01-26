# This script is to normalize four tissue expression matrix.
# Method is TMM + CPM + Log2.

# File in: 1. Four matrix for each tissues gene expression.
#             (expr_leaf, expr_root, expr_sam, expr_seed).
#             From 1_prep_tissue_expr.R
# File out: 1. Four matix of normalized log2(CPM+1) values for each tissue.;
#              (cpm_leaf, cpm_root, cpm_sam and cpm_seed)

################################################################################

# Import libraries.
library(edgeR)


# Normalize for leaf.

y_leaf <- DGEList(counts = expr_leaf,genes = gene_name_all)
keep <- rowSums(cpm(y_leaf) > 1) >= 39 # 39 is 10% of total libraries
y_leaf <- y_leaf[keep,keep.lib.sizes = F]
y_leaf <- calcNormFactors(y_leaf)
gene_name <- as.vector(y_leaf$genes$genes) # updated gene names
cpm_leaf <- cpm(y_leaf, 
                normalized.lib.sizes = TRUE, log = TRUE, prior.count = 1)
rownames(cpm_leaf) <- gene_name

# rm unwanted data
rm(gene_name,keep, y_leaf)


# Normalize for root.

y_root <- DGEList(counts = expr_root,genes = gene_name_all)
keep <- rowSums(cpm(y_root) > 1) >= 18
y_root <- y_root[keep,keep.lib.sizes = F]
y_root <- calcNormFactors(y_root)
gene_name <- as.vector(y_root$genes$genes) # updated gene names
cpm_root <- cpm(y_root,
                normalized.lib.sizes = TRUE, log = TRUE, prior.count = 1)
rownames(cpm_root) <- gene_name

# rm unwanted data
rm(gene_name,keep, y_root)

# Normalize for sam.

y_sam <- DGEList(counts = expr_sam,genes = gene_name_all)
keep <- rowSums(cpm(y_sam) > 1) >= 41
y_sam <- y_sam[keep,keep.lib.sizes = F]
y_sam <- calcNormFactors(y_sam)
gene_name <- as.vector(y_sam$genes$genes) # updated gene names
cpm_sam <- cpm(y_sam,
               normalized.lib.sizes = TRUE, log = TRUE, prior.count = 1)
rownames(cpm_sam) <- gene_name

# rm unwanted data
rm(gene_name,keep, y_sam)

# Normalize for seed.

y_seed <- DGEList(counts = expr_seed,genes = gene_name_all)
keep <- rowSums(cpm(y_seed) > 1) >= 16
y_seed <- y_seed[keep,keep.lib.sizes = F]
y_seed <- calcNormFactors(y_seed)
gene_name <- as.vector(y_seed$genes$genes) # updated gene names
cpm_seed <- cpm(y_seed,
                normalized.lib.sizes = TRUE, log = TRUE, prior.count = 1)
rownames(cpm_seed) <- gene_name

# rm unwanted data
rm(gene_name,keep, y_seed)

# remove raw count data
# cpm_leaf, cpm_root, cpm_sam and cpm_seed were kept
rm(expr_leaf, expr_sam, expr_seed, expr_root)

