library(dplyr)
library(purrr)
library(readxl)
library(clusterProfiler)
library(org.Dr.eg.db)

# read in all data sets for different injury models
cryoinjury <- read_excel('/Users/sophiemarcotte/Desktop/CrossRegen.xlsx', sheet = 'Cryoinjury')
apap <- read_excel('/Users/sophiemarcotte/Desktop/CrossRegen.xlsx', sheet = 'apap')
phx <- read_excel('/Users/sophiemarcotte/Desktop/CrossRegen.xlsx', sheet = 'phx')
mtz <- read_excel('/Users/sophiemarcotte/Desktop/CrossRegen.xlsx', sheet = 'mtz')
heart_valve <- read_excel('/Users/sophiemarcotte/Desktop/CrossRegen.xlsx', sheet = 'heart valve')
caudal_fin <- read_excel('/Users/sophiemarcotte/Desktop/CrossRegen.xlsx', sheet = 'caudal fin')
retina <- read_excel('/Users/sophiemarcotte/Desktop/CrossRegen.xlsx', sheet = 'retina')
spinal_cord <- read_excel('/Users/sophiemarcotte/Desktop/CrossRegen.xlsx', sheet = 'spinal cord')
ventricular_apex <- read_excel('/Users/sophiemarcotte/Desktop/CrossRegen.xlsx', sheet = 'ventricular apex')

# different injury models DOWN ----
# define a function to pull out the down regulated genes from all time points
pull_out_down_genes <- function(df){
  # find columns with "down" in their name
  down_columns <- df %>%
    dplyr::select(contains("down")) %>%
    names()
  # create a list of all elements from those columns
  down_values_list <- df %>%
    dplyr::select(all_of(down_columns)) %>%
    unlist() %>%
    na.omit() 
  return(down_values_list)
}

# run the functions on all injury models
cryoinjury_down <- pull_out_down_genes(cryoinjury)
apap_down <- pull_out_down_genes(apap)
phx_down <- pull_out_down_genes(phx)
mtz_down <- pull_out_down_genes(mtz)

# overlap the lists to get the final lists
overlap_down_genes <- Reduce(intersect, list(cryoinjury_down, apap_down, phx_down, mtz_down))

# different injury models UP ----
# define a function to pull out the up regulated genes from all time points
pull_out_up_genes <- function(df){
  # find columns with "up" in their name
  up_columns <- df %>%
    dplyr::select(contains("up")) %>%
    names()
  # create a list of all elements from those columns
  up_values_list <- df %>%
    dplyr::select(all_of(up_columns)) %>%
    unlist() %>%
    na.omit() 
  return(up_values_list)
}

# run the functions on all injury models
cryoinjury_up <- pull_out_up_genes(cryoinjury)
apap_up <- pull_out_up_genes(apap)
phx_up <- pull_out_up_genes(phx)
mtz_up <- pull_out_up_genes(mtz)

# Combine all gene lists into one vector
all_genes <- c(cryoinjury_up, apap_up, phx_up, mtz_up)
# Count how many times each gene appears
gene_counts <- table(all_genes)
# Keep genes that appear in at least 3 datasets
genes_at_least_3 <- names(gene_counts[gene_counts >= 3])
genes_up_unique <- data.frame(gene_id = unique(genes_at_least_3))


# convert to dataframes
up_df <- as.data.frame(overlap_up_genes)
down_df <- as.data.frame(overlap_down_genes)

# combine both lists with gene symbols
symbol_up <- left_join(genes_up_unique, onlyZeb_conversion, by = c('gene_id' = 'Gene stable ID'))
symbol_down <- left_join(down_df, onlyZeb_conversion, by = c('overlap_down_genes' = 'Gene stable ID'))
symbol_down <- symbol_down[!duplicated(symbol_down$overlap_down_genes), ]

# different organ models DOWN ----
# run the functions on all organ models
heart_valve_down <- pull_out_down_genes(heart_valve)
caudal_fin_down <- pull_out_down_genes(caudal_fin)
retina_down <- pull_out_down_genes(retina)
spinal_cord_down <- pull_out_down_genes(spinal_cord)
ventricular_apex_down <- pull_out_down_genes(ventricular_apex)

# overlap the lists to get the final lists
organ_down_genes <- Reduce(intersect, list(heart_valve_down, caudal_fin_down, retina_down, spinal_cord_down, ventricular_apex_down))

# different organ models UP ----
# run the functions on all organ models
heart_valve_up <- pull_out_up_genes(heart_valve)
caudal_fin_up <- pull_out_up_genes(caudal_fin)
retina_up <- pull_out_up_genes(retina)
spinal_cord_up <- pull_out_up_genes(spinal_cord)
ventricular_apex_up <- pull_out_up_genes(ventricular_apex)

# overlap the lists to get the final lists
organ_up_genes <- Reduce(intersect, list(heart_valve_up, caudal_fin_up, retina_up, spinal_cord_up, ventricular_apex_up))

down_organ_df <- as.data.frame(organ_down_genes)

# combine both lists with gene symbols
symbol_down_organ <- left_join(down_organ_df, conversion, by = c('organ_down_genes' = 'Zebrafish gene stable ID'))

