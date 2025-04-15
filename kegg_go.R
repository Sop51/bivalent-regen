library(clusterProfiler)
library(org.Dr.eg.db)
library(AnnotationDbi)
library(GOplot)

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
                        ont = "BP",  
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

# cryo injury ----
# up
cryo_entrez_up <- convert_to_entrez(cryoinjury_up)
cryo_go_up <- go_bp(cryo_entrez_up)
cryo_kegg_up <- run_kegg(cryo_entrez_up)

# down
cryo_entrez_down <- convert_to_entrez(cryoinjury_down)
cryo_go_down <- go_bp(cryo_entrez_down)
cryo_kegg_down <- run_kegg(cryo_entrez_down)

# apap injury ----
apap_entrez_up <- convert_to_entrez(apap_up)
apap_go_up <- go_bp(apap_entrez_up)
apap_kegg_up <- run_kegg(apap_entrez_up)

# down
apap_entrez_down <- convert_to_entrez(apap_down)
apap_go_down <- go_bp(apap_entrez_down)
apap_kegg_down <- run_kegg(apap_entrez_down)

# phx injury ----
phx_entrez_up <- convert_to_entrez(phx_up)
phx_go_up <- go_bp(phx_entrez_up)
phx_kegg_up <- run_kegg(phx_entrez_up)

# down
phx_entrez_down <- convert_to_entrez(phx_down)
phx_go_down <- go_bp(phx_entrez_down)
phx_kegg_down <- run_kegg(phx_entrez_down)

# mtz injury ----
mtz_entrez_up <- convert_to_entrez(mtz_up)
mtz_go_up <- go_bp(mtz_entrez_up)
mtz_kegg_up <- run_kegg(mtz_entrez_up)

# down
mtz_entrez_down <- convert_to_entrez(mtz_down)
mtz_go_down <- go_bp(mtz_entrez_down)
mtz_kegg_down <- run_kegg(mtz_entrez_down)

# heart valve organ ----
heart_valve_entrez_up <- convert_to_entrez(heart_valve_up)
heart_valve_go_up <- go_bp(heart_valve_entrez_up)
heart_valve_kegg_up <- run_kegg(heart_valve_entrez_up)

# down
heart_valve_entrez_down <- convert_to_entrez(heart_valve_down)
heart_valve_go_down <- go_bp(heart_valve_entrez_down)
heart_valve_kegg_down <- run_kegg(heart_valve_entrez_down)

# caudal fin organ ----
caudal_fin_entrez_up <- convert_to_entrez(caudal_fin_up)
caudal_fin_go_up <- go_bp(caudal_fin_entrez_up)
caudal_fin_kegg_up <- run_kegg(caudal_fin_entrez_up)

# down
caudal_fin_entrez_down <- convert_to_entrez(caudal_fin_down)
caudal_fin_go_down <- go_bp(caudal_fin_entrez_down)
caudal_fin_kegg_down <- run_kegg(caudal_fin_entrez_down)

# retina organ ----
retina_entrez_up <- convert_to_entrez(retina_up)
retina_go_up <- go_bp(retina_entrez_up)
retina_kegg_up <- run_kegg(retina_entrez_up)

# down
retina_entrez_down <- convert_to_entrez(retina_down)
retina_go_down <- go_bp(retina_entrez_down)
retina_kegg_down <- run_kegg(retina_entrez_down)

# ventricular apex organ ----
ventricular_apex_entrez_up <- convert_to_entrez(ventricular_apex_up)
ventricular_apex_go_up <- go_bp(ventricular_apex_entrez_up)
ventricular_apex_kegg_up <- run_kegg(ventricular_apex_entrez_up)

# down
ventricular_apex_entrez_down <- convert_to_entrez(ventricular_apex_down)
ventricular_apex_go_down <- go_bp(ventricular_apex_entrez_down)
ventricular_apex_kegg_down <- run_kegg(ventricular_apex_entrez_down)

# spinal cord organ ----
spinal_cord_entrez_up <- convert_to_entrez(spinal_cord_up)
spinal_cord_go_up <- go_bp(spinal_cord_entrez_up)
spinal_cord_kegg_up <- run_kegg(spinal_cord_entrez_up)

# down
spinal_cord_entrez_down <- convert_to_entrez(spinal_cord_down)
spinal_cord_go_down <- go_bp(spinal_cord_entrez_down)
spinal_cord_kegg_down <- run_kegg(spinal_cord_entrez_down)

