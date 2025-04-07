library(ChIPseeker)
library(rtracklayer)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(org.Mm.eg.db)
library(dplyr)
library(tidyverse)
library(genomation)
library(gUtils)

# load the TxDb object for mouse
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene

# Define column names 
col_names <- c("chr", "start", "end", "ChIP_island_read_count", "CONTROL_island_read_count", "p_value", "fold_change", "FDR_threshold")

# read in the broad peak files for each timepoint of H3K27me3
hr0_broadpeak <- read_delim(
  '/Users/sophiemarcotte/Desktop/0hr_WT_K27_S26_mdup_Addchr-W1000-G4000-islands-summary.txt',
  delim = "\t", 
  col_names = col_names,
  trim_ws = TRUE
)
hr30_broadpeak <- read_delim(
  '/Users/sophiemarcotte/Desktop/GSM3753290_WT-K27_30h_S2_chr-W1000-G4000-islands-summary.txt',
  delim = "\t", 
  col_names = col_names,
  trim_ws = TRUE
) 

hr40_broadpeak <- read_delim(
  '/Users/sophiemarcotte/Desktop/GSM3753291_WT-K27_40h_S3_chr-W1000-G4000-islands-summary.txt',
  delim = "\t", 
  col_names = col_names,
  trim_ws = TRUE
) 

hr96_broadpeak <- read_delim(
  '/Users/sophiemarcotte/Desktop/GSM3753292_WT-K27_96h_S4_chr-W1000-G4000-islands-summary.txt',
  delim = "\t", 
  col_names = col_names,
  trim_ws = TRUE
) 

# turn each dataframe to a granges object
hr0 <- GRanges(
  seqnames = hr0_broadpeak$chr,
  ranges = IRanges(
    start = hr0_broadpeak$start,
    end = hr0_broadpeak$end
  ),
  ChIP_island_read_count = hr0_broadpeak$ChIP_island_read_count,
  CONTROL_island_read_count = hr0_broadpeak$CONTROL_island_read_count,
  p_value = hr0_broadpeak$p_value,
  fold_change = hr0_broadpeak$fold_change
)

hr30 <- GRanges(
  seqnames = hr30_broadpeak$chr,
  ranges = IRanges(
    start = hr30_broadpeak$start,
    end = hr30_broadpeak$end
  ),
  ChIP_island_read_count = hr30_broadpeak$ChIP_island_read_count,
  CONTROL_island_read_count = hr30_broadpeak$CONTROL_island_read_count,
  p_value = hr30_broadpeak$p_value,
  fold_change = hr30_broadpeak$fold_change
)

hr40 <- GRanges(
  seqnames = hr40_broadpeak$chr,
  ranges = IRanges(
    start = hr40_broadpeak$start,
    end = hr40_broadpeak$end
  ),
  ChIP_island_read_count = hr40_broadpeak$ChIP_island_read_count,
  CONTROL_island_read_count = hr40_broadpeak$CONTROL_island_read_count,
  p_value = hr40_broadpeak$p_value,
  fold_change = hr40_broadpeak$fold_change
)

hr96 <- GRanges(
  seqnames = hr96_broadpeak$chr,
  ranges = IRanges(
    start = hr96_broadpeak$start,
    end = hr96_broadpeak$end
  ),
  ChIP_island_read_count = hr96_broadpeak$ChIP_island_read_count,
  CONTROL_island_read_count = hr96_broadpeak$CONTROL_island_read_count,
  p_value = hr96_broadpeak$p_value,
  fold_change = hr96_broadpeak$fold_change
)

# read in the H3K4me3 peaks
H3K4me3_peaks = readNarrowPeak('/Users/sophiemarcotte/Desktop/GSM3561048_Sample_WTB-H3K4me3_peaks.narrowPeak')
# rename the chromosomes
H3K4me3_peaks <- gr.chr(H3K4me3_peaks)

# --------------- extracting a list of bivalent genes across all timepoints ------------#
############# join hr0 27me3 with 4me3 ###############
hr0_bivalent_overlapping <- findOverlaps(hr0, H3K4me3_peaks)
hr0_bivalent_overlap <- H3K4me3_peaks[subjectHits(hr0_bivalent_overlapping)]

# annotate the peaks again
hr0_bivalent_annotation <- annotatePeak(
  hr0_bivalent_overlap,
  tssRegion = c(-5000, 5000),
  TxDb = txdb,
  assignGenomicAnnotation = TRUE,
  annoDb = "org.Mm.eg.db"
)

# convert to a dataframe
hr0_bivalent_df <- as.data.frame(hr0_bivalent_annotation)

############# join hr30 27me3 with 4me3 ###############
hr30_bivalent_overlapping <- findOverlaps(hr30, H3K4me3_peaks)
hr30_bivalent_overlap <- H3K4me3_peaks[subjectHits(hr30_bivalent_overlapping)]

# annotate the peaks again
hr30_bivalent_annotation <- annotatePeak(
  hr30_bivalent_overlap,
  tssRegion = c(-5000, 5000),
  TxDb = txdb,
  assignGenomicAnnotation = TRUE,
  annoDb = "org.Mm.eg.db"
)

# convert to a dataframe
hr30_bivalent_df <- as.data.frame(hr30_bivalent_annotation)

############# join hr40 27me3 with 4me3 ###############
hr40_bivalent_overlapping <- findOverlaps(hr40, H3K4me3_peaks)
hr40_bivalent_overlap <- H3K4me3_peaks[subjectHits(hr40_bivalent_overlapping)]

# annotate the peaks again
hr40_bivalent_annotation <- annotatePeak(
  hr40_bivalent_overlap,
  tssRegion = c(-5000, 5000),
  TxDb = txdb,
  assignGenomicAnnotation = TRUE,
  annoDb = "org.Mm.eg.db"
)

# convert to a dataframe
hr40_bivalent_df <- as.data.frame(hr40_bivalent_annotation)

############# join hr96 27me3 with 4me3 ###############
hr96_bivalent_overlapping <- findOverlaps(hr96, H3K4me3_peaks)
hr96_bivalent_overlap <- H3K4me3_peaks[subjectHits(hr96_bivalent_overlapping)]

# annotate the peaks again
hr96_bivalent_annotation <- annotatePeak(
  hr96_bivalent_overlap,
  tssRegion = c(-5000, 5000),
  TxDb = txdb,
  assignGenomicAnnotation = TRUE,
  annoDb = "org.Mm.eg.db"
)

# convert to a dataframe
hr96_bivalent_df <- as.data.frame(hr96_bivalent_annotation)

# find gene dynamics that are changing overtime
# get unique ENSEMBL gene lists for each time point
genes_hr0 <- unique(hr0_bivalent_df$ENSEMBL)
genes_hr30 <- unique(hr30_bivalent_df$ENSEMBL)
genes_hr40 <- unique(hr40_bivalent_df$ENSEMBL)
genes_hr96 <- unique(hr96_bivalent_df$ENSEMBL)

# get ALL bivalent genes in one list
all_genes <- unique(c(genes_hr0, genes_hr30, genes_hr40, genes_hr96))
# save to CSV
write.csv(data.frame(ENSEMBL = all_genes), "/Users/sophiemarcotte/Desktop/unique_ensembl_bivalent.csv", row.names = FALSE)

# to only find those losing the methyl mark
# find genes present in hr0 but NOT in hr30
lost_hr30 <- setdiff(genes_hr0, genes_hr30)

# find genes present in hr30 but NOT in hr40
lost_hr40 <- setdiff(genes_hr30, genes_hr40)

# find genes present in hr40 but NOT in hr96
lost_hr96 <- setdiff(genes_hr40, genes_hr96)

# combine all lost genes into one list
all_bivalent_regen_genes <- c(lost_hr30, lost_hr40, lost_hr96)

# get unique genes from the combined list
unique_bivalent_genes <- unique(all_bivalent_regen_genes)

# save to CSV
write.csv(unique_bivalent_genes, "bivalent_regen_genes.csv", row.names = FALSE)
