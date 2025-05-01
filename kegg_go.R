library(clusterProfiler)
library(org.Dr.eg.db)
library(AnnotationDbi)
library(GOplot)
library(biomaRt)
library(readr)

# helper functions ----
# convert to entrez ids
convert_to_entrez <- function(gene_list){
  entrez <- bitr(gene_list, fromType = "ENSEMBL",
                  toType = "ENTREZID",
                  OrgDb = org.Dr.eg.db)
  return(entrez)
}

# run go BP
go_bp <- function(entrez_ids){
  go_enrich <- enrichGO(gene = entrez_ids$ENTREZID,
                        OrgDb = org.Dr.eg.db,
                        ont = "CC",  
                        pAdjustMethod = "BH",
                        pvalueCutoff = 0.05,
                        readable = TRUE)
  # create a barplot or dotplot
  go_plot <- dotplot(go_enrich, showCategory = 10, font.size = 10) 
  # return both result and plot
  return(list(results = go_enrich@result, plot = go_plot))
}

# run kegg 
run_kegg <- function(entriz_ids){
  kegg_enrich <- enrichKEGG(gene = entriz_ids$ENTREZID,
                            organism = 'dre',  
                            pvalueCutoff = 0.05)
  kegg_enrich <- setReadable(kegg_enrich, OrgDb = org.Dr.eg.db, keyType = "ENTREZID")
  # create a plot
  kegg_plot <- dotplot(kegg_enrich, showCategory = 20, font.size = 10) + ggtitle("KEGG Pathways")
  # return both result and plot
  return(list(results = kegg_enrich@result, plot = kegg_plot))
}

# get gene symbols from kegg results for a pathway
get_symbols_from_description <- function(df, description) {
  # filter for the row that matches the desired Description
  row <- df %>% filter(Description == description)
  # extract ENTREZ IDs from 'geneID' column
  entrez_ids <- unlist(strsplit(row$geneID, "/"))
  # map ENTREZ to gene symbols
  symbols_df <- bitr(entrez_ids,
                     fromType = "ENTREZID",
                     toType = "SYMBOL",
                     OrgDb = org.Dr.eg.db)
  
  return(symbols_df$SYMBOL)
}

# mtz/phx injury upregulated
phxmtx_entrez_up <- convert_to_entrez(overlap_up_unique)
phxmtx_go_up <- go_bp(phxmtx_entrez_up)
phxmtx_go_up$plot 
phxmtx_kegg_up <- run_kegg(phxmtx_entrez_up)
phxmtx_kegg_up$plot 

# all injury up
injury_entrez_up <- convert_to_entrez(up_genes_least_3)
injury_go_up <- go_bp(injury_entrez_up)
injury_go_up$plot # none significant
injury_kegg_up <- run_kegg(injury_entrez_up)
injury_kegg_up$plot # p53 pathway

# all injury down
injury_entrez_down <- convert_to_entrez(down_genes_least_3)
injury_go_down <- go_bp(injury_entrez_down)
injury_go_down$plot
injury_kegg_down <- run_kegg(injury_entrez_down)
injury_kegg_down$plot

# all organ up
organ_entrez_up <- convert_to_entrez(up_organ_genes_least_3)
organ_go_up <- go_bp(organ_entrez_up)
organ_go_up$plot 
organ_kegg_up <- run_kegg(organ_entrez_up)
organ_kegg_up$plot # none

# all organ down
organ_entrez_down <- convert_to_entrez(down_organ_genes_least_3)
organ_go_down <- go_bp(organ_entrez_down)
organ_go_down$plot 
organ_kegg_down <- run_kegg(organ_entrez_down)
organ_kegg_down$plot # none

# convert a list of genes to human orthologs ----
orthologs <- read_tsv('/Users/sm2949/Desktop/mart_export.txt')
zebrafish_genes_df_phxmtz <- data.frame(Zebrafish_Ensembl = overlap_up_unique)
annotated_genes <- left_join(zebrafish_genes_df_phxmtz, orthologs, by = c("Zebrafish_Ensembl"="Gene stable ID"))
writeLines(annotated_genes$`Human gene name`, "/Users/sm2949/Desktop/MtzPhxUpGenesHumanOrthologs.txt")
