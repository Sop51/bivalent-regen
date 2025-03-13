library(readr)
library(biomaRt)
library(dplyr)
library(edgeR)

# read in the list of bivalent mice genes
bivalent_mice_genes <- read_csv('/Users/sophiemarcotte/Desktop/bivalent_regen_genes.csv')

# read in the conversion file
conversion <- read_tsv('/Users/sophiemarcotte/Desktop/mart_export.txt')

# join the tables together to map orthologs
bivalent_orthologs <- left_join(bivalent_mice_genes, conversion, by = c("x" = "Gene stable ID"))

# -------------------- cryoinjury model --------------------#
# read in the normalized counts file
cryo_counts <- read_csv('/Users/sophiemarcotte/Desktop/GSE245878_adultliverzbregeneration_normalisedcounts.csv')

# inner join with the bivalent genes
bivalent_cryo <- inner_join(cryo_counts, bivalent_orthologs, by = c("LLgeneAbbrev" = "Zebrafish gene name"))

bivalent_cryo <- bivalent_cryo %>% distinct(LLgeneAbbrev, .keep_all = TRUE)

# create colData 
samples <- c("R63-SHAM CPM", "R62-SHAM CPM", "R61-SHAM CPM", "R64-SHAM CPM", "R1-SHAM CPM",
             "R33-3-dpci CPM", "R29-3-dpci CPM", "R31-3-dpci CPM", "R27-3-dpci CPM",
             "R17-1dpci", "R21-1dpci", "R22-1dpci", "R18-1dpci",
             "R9-7-dpci", "31-7-dpci", "R11-7-dpci", "30-7-dpci")

groups <- as.factor(c(ifelse(grepl("SHAM", samples), "SHAM",
          ifelse(grepl("3-dpci", samples), "3dpci",
          ifelse(grepl("1dpci", samples), "1dpci",
          ifelse(grepl("7-dpci", samples), "7dpci", "Unknown"))))))

countData <- bivalent_cryo[,c(21, 3:19)]
countData <- as.data.frame(countData)
rownames(countData) <- countData[,1]  
countData <- countData[,-1]

# create the DGElist object
y <- DGEList(counts=countData,group=groups)
design <- model.matrix(~0+groups)
# estimate dispersion
y <- estimateDisp(y, design)
colnames(design) <- levels(groups)
# fit the model
fit <- glmQLFit(y, design, robust=TRUE)

qlf_timepoint1 <- glmQLFTest(fit, contrast = c(1, 0, 0,-1))
qlf_timepoint2 <- glmQLFTest(fit, contrast = c(0, 1, 0,-1))
qlf_timepoint3 <- glmQLFTest(fit, contrast = c(0, 0, 1,-1))

dcpi1 <- as.data.frame(qlf_timepoint1)
dcpi3 <- as.data.frame(qlf_timepoint2)
dcpi7 <- as.data.frame(qlf_timepoint3)

# filter for significant genes

