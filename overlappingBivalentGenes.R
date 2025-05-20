library(dplyr)
library(purrr)
library(readxl)
library(clusterProfiler)
library(org.Dr.eg.db)
library(ggplot2)
library(biomaRt)

# read in all data sets for different injury models
cryoinjury <- read_excel('/Users/sm2949/Desktop/BivalentRegen.xlsx', sheet = 'Cryoinjury')
apap <- read_excel('/Users/sm2949/Desktop/BivalentRegen.xlsx', sheet = 'APAP')
phx <- read_excel('/Users/sm2949/Desktop/BivalentRegen.xlsx', sheet = 'PHx')
mtz <- read_excel('/Users/sm2949/Desktop/BivalentRegen.xlsx', sheet = 'Mtz')
heart_valve <- read_excel('/Users/sm2949/Desktop/BivalentRegen.xlsx', sheet = 'Heart Valve')
caudal_fin <- read_excel('/Users/sm2949/Desktop/BivalentRegen.xlsx', sheet = 'Caudal Fin')
retina <- read_excel('/Users/sm2949/Desktop/BivalentRegen.xlsx', sheet = 'Retina')
spinal_cord <- read_excel('/Users/sm2949/Desktop/BivalentRegen.xlsx', sheet = 'Spinal Cord')
ventricular_apex <- read_excel('/Users/sm2949/Desktop/BivalentRegen.xlsx', sheet = 'Ventricular Apex')

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

# overlap the lists to get the final lists of ones in ALL lists
overlap_down_genes <- Reduce(intersect, list(cryoinjury_down, apap_down, phx_down, mtz_down))

# overlap lists to get genes in at least 3 of the models 
unique_down_lists <- list(
  unique(cryoinjury_down),
  unique(apap_down),
  unique(phx_down),
  unique(mtz_down)
)
# combine all deduplicated genes into one vector
all_down_genes <- unlist(unique_down_lists)
# count how many times each gene appears
down_gene_counts <- table(all_down_genes)
# get genes that appear in at least 3 models
down_genes_least_3 <- names(down_gene_counts[down_gene_counts >= 3])

# save list
writeLines(down_genes_least_3, "/Users/sm2949/Desktop/acrossInjuryDownGenes.txt")

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

# overlap the lists to get the final lists of ONLY phx and mtz
overlap_up_genes <- Reduce(intersect, list(phx_up, mtz_up, cryoinjury_up, apap_up))
overlap_up_unique <- unique(overlap_up_genes)
# save list of ONLY phx and mtz overalps
writeLines(overlap_up_unique, "/Users/sm2949/Desktop/MtzPhxUpGenes.txt")


# overlap lists to get genes in at least 3 of the models 
unique_up_lists <- list(
  unique(cryoinjury_up),
  unique(apap_up),
  unique(phx_up),
  unique(mtz_up)
)
# combine all deduplicated genes into one vector
all_up_genes <- unlist(unique_up_lists)
# count how many times each gene appears
up_gene_counts <- table(all_up_genes)
# get genes that appear in at least 3 models
up_genes_least_3 <- names(up_gene_counts[up_gene_counts >= 3])

# save list
writeLines(up_genes_least_3, "/Users/sm2949/Desktop/acrossInjuryUpGenes.txt")

# different organ models DOWN ----
# run the functions on all organ models
heart_valve_down <- pull_out_down_genes(heart_valve)
caudal_fin_down <- pull_out_down_genes(caudal_fin)
retina_down <- pull_out_down_genes(retina)
spinal_cord_down <- pull_out_down_genes(spinal_cord)
ventricular_apex_down <- pull_out_down_genes(ventricular_apex)

# overlap the lists to get the final lists
organ_down_genes <- Reduce(intersect, list(heart_valve_down, caudal_fin_down, retina_down, spinal_cord_down, ventricular_apex_down))

# overlap lists to get genes in at least 3 of the models 
unique_down_organ_lists <- list(
  unique(heart_valve_down),
  unique(caudal_fin_down),
  unique(retina_down),
  unique(spinal_cord_down),
  unique(ventricular_apex_down)
)
# combine all deddownlicated genes into one vector
all_down_organ_genes <- unlist(unique_down_organ_lists)
# count how many times each gene appears
down_organ_gene_counts <- table(all_down_organ_genes)
# get genes that appear in at least 3 models
down_organ_genes_least_3 <- names(down_organ_gene_counts[down_organ_gene_counts >= 3])

# save list
writeLines(down_organ_genes_least_3, "/Users/sm2949/Desktop/acrossOrganDownGenes.txt")

# different organ models UP ----
# run the functions on all organ models
heart_valve_up <- pull_out_up_genes(heart_valve)
caudal_fin_up <- pull_out_up_genes(caudal_fin)
retina_up <- pull_out_up_genes(retina)
spinal_cord_up <- pull_out_up_genes(spinal_cord)
ventricular_apex_up <- pull_out_up_genes(ventricular_apex)

# overlap the lists to get the final lists
organ_up_genes <- Reduce(intersect, list(heart_valve_up, caudal_fin_up, retina_up, spinal_cord_up, ventricular_apex_up))

# overlap lists to get genes in at least 3 of the models 
unique_up_organ_lists <- list(
  unique(heart_valve_up),
  unique(caudal_fin_up),
  unique(retina_up),
  unique(spinal_cord_up),
  unique(ventricular_apex_up)
)
# combine all deduplicated genes into one vector
all_up_organ_genes <- unlist(unique_up_organ_lists)
# count how many times each gene appears
up_organ_gene_counts <- table(all_up_organ_genes)
# get genes that appear in at least 3 models
up_organ_genes_least_3 <- names(up_organ_gene_counts[up_organ_gene_counts >= 3])

# save list
writeLines(up_organ_genes_least_3, "/Users/sm2949/Desktop/acrossOrganUpGenes.txt")

# create a data frame for plotting ----
gene_counts <- data.frame(
  Regulation = c("Upregulated", "Downregulated"),
  Count = c(length(overlap_up_genes), length(overlap_down_genes))
)

# define custom colors
custom_colors <- c("Upregulated" = "#800020",   # Burgundy red
                   "Downregulated" = "#000080") # Navy blue

# plot
ggplot(gene_counts, aes(x = Regulation, y = Count, fill = Regulation)) +
  geom_col(width = 0.5, show.legend = FALSE) +
  scale_fill_manual(values = custom_colors) +
  geom_text(aes(label = Count),
            vjust = -0.3,
            size = 5,
            color = "black") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +  # Add space at top
  labs(
    title = "Shared Bivalent Gene Regulation Across Injury Models",
    x = NULL,
    y = "Number of Genes"
  ) +
  theme_minimal(base_size = 12) +
  theme(    
    panel.grid.minor = element_blank(),                                                                
    axis.text.x = element_text(color = "black", size = 12),
    axis.text.y = element_text(size = 12, color = "black"),
    axis.title.y = element_text(size = 12),
    plot.title = element_text(face = "bold", hjust = 0.5, size = 12),
    panel.border = element_rect(color = "black", fill = NA, size = 1)   # black box around panel
  )

# convert to gene symbols to create a table ---
ensembl <- useEnsembl(biomart = "ensembl", dataset = "drerio_gene_ensembl")
converted_up_genes <- getBM(
  attributes = c("ensembl_gene_id", "external_gene_name", "description"),
  filters = "ensembl_gene_id",
  values = overlap_up_genes,
  mart = ensembl
)

# clean up
converted_up_genes$description <- gsub("\\ \\[.*\\]", "", converted_up_genes$description)
# rename columns for clarity
colnames(converted_up_genes) <- c("Ensembl Gene ID", "Gene Symbol", "Description")

# Sort alphabetically by gene symbol
converted_genes <- converted_up_genes[order(converted_up_genes$`Gene Symbol`), ]
rownames(converted_genes) <- NULL

formattable(converted_genes, list(
  area(row = which(converted_genes$`Gene Symbol` == "cbx7a")) ~ 
    formatter("span", 
              style = ~ style(background = "#ffefcc", font.weight = "bold"))
))
