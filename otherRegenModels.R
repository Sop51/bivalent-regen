library(readr)
library(biomaRt)
library(dplyr)

# read in the list of bivalent mice genes
bivalent_mice_genes <- read_csv('/Users/sophiemarcotte/Desktop/bivalent_regen_genes.csv')

# read in the conversion file
conversion <- read_tsv('/Users/sophiemarcotte/Desktop/mart_export.txt')

# join the tables together to map orthologs
bivalent_orthologs <- left_join(bivalent_mice_genes, conversion, by = c("x" = "Gene stable ID"))
