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
  go_df <- go_enrich@result
  return(go_df)
}

# run kegg 
run_kegg <- function(entriz_ids){
  kegg_enrich <- enrichKEGG(gene = entriz_ids$ENTREZID,
                            organism = 'dre',  
                            pvalueCutoff = 0.05)
  kegg_df <- kegg_enrich@result
  return(kegg_df)
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
injury_entrez_up <- convert_to_entrez(injury_up_genes)
injury_go_up <- go_bp(injury_entrez_up)
injury_kegg_up <- run_kegg(injury_entrez_up)

# all injury down

# all organ up

# all organ down

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

# overlapping kegg terms
overlap_up_kegg <- Reduce(intersect, list(mtz_kegg_up$Description, apap_kegg_up$Description, 
                                          phx_kegg_up$Description, cryo_kegg_up$Description,
                                          heart_valve_kegg_up$Description, caudal_fin_kegg_up$Description,
                                          retina_kegg_up$Description, spinal_cord_kegg_up$Description,
                                          ventricular_apex_kegg_up$Description))


# looking at pathway genes from the overlapping pathways
# apelin signaling pathway
apelin_signaling_pathway_cryo <- get_symbols_from_description(cryo_kegg_up, "Apelin signaling pathway")
apelin_signaling_pathway_apap <- get_symbols_from_description(apap_kegg_up, "Apelin signaling pathway")
apelin_signaling_pathway_phx <- get_symbols_from_description(phx_kegg_up, "Apelin signaling pathway")
apelin_signaling_pathway_mtz <- get_symbols_from_description(mtz_kegg_up, "Apelin signaling pathway")
# combine
apelin_signaling_combined <- unique(c(
  apelin_signaling_pathway_cryo,
  apelin_signaling_pathway_apap,
  apelin_signaling_pathway_phx,
  apelin_signaling_pathway_mtz
))

# ecm-receptor interaction
ecm_receptor_interaction_cryo <- get_symbols_from_description(cryo_kegg_up, "ECM-receptor interaction")
ecm_receptor_interaction_apap <- get_symbols_from_description(apap_kegg_up, "ECM-receptor interaction")
ecm_receptor_interaction_phx <- get_symbols_from_description(phx_kegg_up, "ECM-receptor interaction")
ecm_receptor_interaction_mtz <- get_symbols_from_description(mtz_kegg_up, "ECM-receptor interaction")
# combine
ecm_receptor_combined <- unique(c(
  ecm_receptor_interaction_cryo,
  ecm_receptor_interaction_apap,
  ecm_receptor_interaction_phx,
  ecm_receptor_interaction_mtz
))

# salmonella infection
salmonella_infection_cryo <- get_symbols_from_description(cryo_kegg_up, "Salmonella infection")
salmonella_infection_apap <- get_symbols_from_description(apap_kegg_up, "Salmonella infection")
salmonella_infection_phx <- get_symbols_from_description(phx_kegg_up, "Salmonella infection")
salmonella_infection_mtz <- get_symbols_from_description(mtz_kegg_up, "Salmonella infection")
# combine
salmonella_combined <- unique(c(
  salmonella_infection_cryo,
  salmonella_infection_apap,
  salmonella_infection_phx,
  salmonella_infection_mtz
))



