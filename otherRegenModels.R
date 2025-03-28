library(readr)
library(biomaRt)
library(dplyr)
library(edgeR)

# read in the list of bivalent mice genes
bivalent_mice_genes <- read_csv('/Users/sm2949/Desktop/bivalent_regen_genes.csv')
# read in the conversion file
conversion <- read_tsv('/Users/sm2949/Desktop/mart_export.txt')

# join the tables together to map orthologs
bivalent_orthologs <- left_join(bivalent_mice_genes, conversion, by = c("x" = "Gene stable ID"))

# -------------------- cryoinjury model --------------------#
# read in the normalized counts file
cryo_counts <- read_csv('/Users/sophiemarcotte/Desktop/GSE245878_adultliverzbregeneration.csv')

# inner join with the bivalent genes
bivalent_cryo <- inner_join(cryo_counts, bivalent_orthologs, by = c("LLgeneAbbrev" = "Zebrafish gene name"))

bivalent_cryo <- bivalent_cryo %>% distinct(LLgeneAbbrev, .keep_all = TRUE)

# create colData 
samples <- c("R63-SHAM", "R62-SHAM", "R61-SHAM", "R64-SHAM", "R1-SHAM",
             "R33-3-dpci", "R29-3-dpci", "R31-3-dpci", "R27-3-dpci",
             "R17-1dpci", "R21-1dpci", "R22-1dpci", "R18-1dpci",
             "R9-7-dpci", "31-7-dpci", "R11-7-dpci", "30-7-dpci")

groups <- as.factor(c(ifelse(grepl("SHAM", samples), "SHAM",
          ifelse(grepl("3-dpci", samples), "3dpci",
          ifelse(grepl("1dpci", samples), "1dpci",
          ifelse(grepl("7-dpci", samples), "7dpci", "Unknown"))))))

# create count matrix
countData <- bivalent_cryo[,c(21, 3:19)]
countData <- as.data.frame(countData)
rownames(countData) <- countData[,1]  
countData <- countData[,-1]

# create the DGElist object
y <- DGEList(counts=countData,group=groups)
design <- model.matrix(~0+groups)
# normalize the data
y <- calcNormFactors(y)
# estimate dispersion
y <- estimateDisp(y, design)
colnames(design) <- levels(groups)
# fit the model
fit <- glmQLFit(y, design, robust=TRUE)

# fit test
qlf_timepoint1 <- glmQLFTest(fit, contrast = c(1, 0, 0,-1))
qlf_timepoint2 <- glmQLFTest(fit, contrast = c(0, 1, 0,-1))
qlf_timepoint3 <- glmQLFTest(fit, contrast = c(0, 0, 1,-1))

# save results as a dataframe
dcpi1 <- as.data.frame(qlf_timepoint1)
dcpi3 <- as.data.frame(qlf_timepoint2)
dcpi7 <- as.data.frame(qlf_timepoint3)

# filter for significant genes
dcpi1_significant <- dcpi1[dcpi1$PValue < 0.05, ]
dcpi3_significant <- dcpi3[dcpi3$PValue < 0.05, ]
dcpi7_significant <- dcpi7[dcpi7$PValue < 0.05, ]

# split to up and down reg at each timepoint
up1dpi <- dcpi1_significant[dcpi1_significant$logFC > 1, ]
down1dpi <- dcpi1_significant[dcpi1_significant$logFC < -1, ]

up3dpi <- dcpi3_significant[dcpi3_significant$logFC > 1, ]
down3dpi <- dcpi3_significant[dcpi3_significant$logFC < -1, ]

up7dpi <- dcpi7_significant[dcpi7_significant$logFC > 1, ]
down7dpi <- dcpi7_significant[dcpi7_significant$logFC < -1, ]


# DCPI1
# convert rownames of dcpi1_significant into a column
dcpi1_significant <- dcpi1_significant %>% 
  tibble::rownames_to_column(var = "gene_id")
# perform left join
sig_dcpi1 <- left_join(dcpi1_significant, bivalent_cryo, 
                       by = c("gene_id" = "Zebrafish gene stable ID"))
     
# DCPI3
# convert rownames of dcpi1_significant into a column
dcpi3_significant <- dcpi3_significant %>% 
  tibble::rownames_to_column(var = "gene_id")
# perform left join
sig_dcpi3 <- left_join(dcpi3_significant, bivalent_cryo, 
                       by = c("gene_id" = "Zebrafish gene stable ID"))

# DCPI7
# convert rownames of dcpi1_significant into a column
dcpi3_significant <- dcpi3_significant %>% 
  tibble::rownames_to_column(var = "gene_id")
# perform left join
sig_dcpi3 <- left_join(dcpi3_significant, bivalent_cryo, 
                       by = c("gene_id" = "Zebrafish gene stable ID"))

# ---------------------- APAP Injury Model -------------------------- #
# read in the data
apap_counts <- read_tsv('/Users/sophiemarcotte/Desktop/GSE230829_RNAseq_counts_raw.txt')
apap_filtered <- apap_counts[,c(1, 18:33)]

apap_samples <- c("WA-12-1_S29","WA-12-2_S30","WA-12-3_S31","WA-12-4_S32","WA-24-1_S45",
                  "WA-24-2_S46","WA-24-3_S47","WA-24-4_S48","WD-12-1_S25","WD-12-2_S26",
                  "WD-12-3_S27","WD-12-4_S28","WD-24-1_S41","WD-24-2_S42","WD-24-3_S43","WD-24-4_S44")

groups_apap <- as.factor(c(ifelse(grepl("WA-12", apap_samples), "APAP-12",
                           ifelse(grepl("WA-24", apap_samples), "APAP-24",
                           ifelse(grepl("WD-12", apap_samples), "Control-12",
                           ifelse(grepl("WD-24", apap_samples), "Control-24", "Unknown"))))))

# create count matrix
countData_apap <-apap_filtered[,c(1:17)]
countData_apap <- countData_apap %>% distinct(...1, .keep_all = TRUE)
countData_apap <- as.data.frame(countData_apap)
rownames(countData_apap) <- countData_apap[,1]  
countData_apap <- countData_apap[,-1]

# create the DGElist object
y <- DGEList(counts=countData_apap, group=groups_apap)
# filter out low expressed genes
keep <- filterByExpr(y)
y <- y[keep, , keep.lib.sizes=FALSE]
# create the design matrix
design <- model.matrix(~0+groups_apap)
# normalize the data
y <- calcNormFactors(y)
# estimate dispersion
y <- estimateDisp(y, design)
colnames(design) <- levels(groups_apap)
# fit the model
fit <- glmQLFit(y, design, robust=TRUE)

# fit test
qlf_timepoint12 <- glmQLFTest(fit, contrast = c(1, 0, -1, 0))
qlf_timepoint24 <- glmQLFTest(fit, contrast = c(0, 1, 0,-1))

# save results as a dataframe
apap12 <- as.data.frame(qlf_timepoint12)
apap24 <- as.data.frame(qlf_timepoint24)

# filter for significant genes
apap12_significant <- apap12[apap12$PValue < 0.05, ]
apap24_significant <- apap24[apap24$PValue < 0.05, ]

# split to up and down reg at each timepoint
up.apap12 <- apap12_significant[apap12_significant$logFC > 1, ]
up.apap24 <- apap24_significant[apap24_significant$logFC > 1, ]

# filter for only bivalent genes
up.apap12 <- up.apap12[rownames(up.apap12) %in% zebrafish_bivalent, ]
up.apap24 <- up.apap24[rownames(up.apap24) %in% zebrafish_bivalent, ]

# ---------------------------------- PHx injury model ----------------------------------#
phx_counts <- read_tsv('/Users/sophiemarcotte/Desktop/fragmentCounts.txt', skip = 1)
# fix the col names 
colnames(phx_counts) <- gsub("^.*/(.*)\\.sorted\\.bam$", "\\1", colnames(phx_counts))
# only select relevant columns 
phx_counts <- phx_counts[,c(1,7:24)]
# set the gene names as the col names
phx_counts <- as.data.frame(phx_counts)
rownames(phx_counts) <- phx_counts[,1]  
phx_counts <- phx_counts[,-1]

# create the metadata
phx_samples <- colnames(phx_counts)
groups_phx <- as.factor(c(ifelse(grepl("120", phx_samples), "120hr",
                         ifelse(grepl("168", phx_samples), "168hr",
                         ifelse(grepl("24", phx_samples), "24hr",
                         ifelse(grepl("6hr", phx_samples), "6hr",
                         ifelse(grepl("72", phx_samples), "72hr", 
                         ifelse(grepl("sham", phx_samples), "Control", "Unknown"))))))))


# create the DGElist object
y <- DGEList(counts=phx_counts, group=groups_phx)
# filter out low expressed genes
keep <- filterByExpr(y)
y <- y[keep, , keep.lib.sizes=FALSE]
# create the design model
design <- model.matrix(~0+groups_phx)
# normalize the data
y <- calcNormFactors(y)
# estimate dispersion
y <- estimateDisp(y, design)
colnames(design) <- levels(groups_phx)
# fit the model
fit <- glmQLFit(y, design, robust=TRUE)

# fit test
qlf_timepoint120 <- glmQLFTest(fit, contrast = c(1, 0, 0, 0, 0, -1))
qlf_timepoint168 <- glmQLFTest(fit, contrast = c(0, 1, 0, 0, 0, -1))
qlf_timepoint24 <- glmQLFTest(fit, contrast = c(0, 0, 1, 0, 0, -1))
qlf_timepoint6 <- glmQLFTest(fit, contrast = c(0, 0, 0, 1, 0, -1))
qlf_timepoint72 <- glmQLFTest(fit, contrast = c(0, 0, 0, 0, 1, -1))

# save results as a dataframe
phx6 <- as.data.frame(qlf_timepoint6)
phx24 <- as.data.frame(qlf_timepoint24)
phx72 <- as.data.frame(qlf_timepoint72)
phx120 <- as.data.frame(qlf_timepoint120)
phx168 <- as.data.frame(qlf_timepoint168)

# define the filtering function
filter_significant_genes <- function(df, gene_list) {
  df <- df[df$PValue < 0.05 & df$logFC > 1, ]  # filter for significance and logFC
  df <- df[rownames(df) %in% gene_list, ]  # only bivalent genes
  return(df)
}

# apply the function to each dataset
phx6_filtered <- filter_significant_genes(phx6, zebrafish_bivalent)
phx24_filtered <- filter_significant_genes(phx24, zebrafish_bivalent)
phx72_filtered <- filter_significant_genes(phx72, zebrafish_bivalent)
phx120_filtered <- filter_significant_genes(phx120, zebrafish_bivalent)
phx168_filtered <- filter_significant_genes(phx168, zebrafish_bivalent)
