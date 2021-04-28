#!/usr/bin/Rscript

# --------------------------------------------------------------------------
#! Perform patwhay enrichment in the human chondrocyte signature
# --------------------------------------------------------------------------

# Working directory
if(dir.exists("/Users/mac/cloudstor/")) {
  setwd("/Users/mac/cloudstor/scott_projects/osteochondro_signature/") # Uses practice data (5% of cells from each sample) if running locally
  place <- "local"
} else {
  setwd("/share/ScratchGeneral/scoyou/sarah_projects/osteochondro_signature/")
  place <- "wolfpack"
}

# Libraries
library('dplyr')
library('tidyr')
library('ggplot2')
library('readr')
library('Matrix')
library('clusterProfiler')
library('ReactomePA')
library('org.Hs.eg.db')
library('gplots')
library('msigdbr')
library('gtools')
library('ggcorrplot')

# Load gene annotation info
geneActivity <- read.csv('project_results/humanchondro_bulk/Chondrocyte_signature_limma.csv', row.names = 1, header = TRUE)

# Perform enrichment analysis on each element of list
universe <- as.character(geneActivity[geneActivity$Biotype == "protein_coding", "GeneSymbol"])
universe <- universe[!is.na(universe)]
universe_entrez <- as.character(geneActivity[geneActivity$Biotype == "protein_coding", "EntrezID"])
universe_entrez <- universe_entrez[!is.na(universe_entrez)]

# Makes HPO geneset for enrichment testing
h_df <- msigdbr(species = "Homo sapiens", category = "C5", subcategory = "HPO")
h_t2g <- h_df %>% dplyr::select(gs_name, gene_symbol) %>% as.data.frame()

# Makes TF geneset for enrichment testing
htf_df <- msigdbr(species = "Homo sapiens", category = "C3", subcategory = 'TFT:GTRD')
htf_t2g <- htf_df %>% dplyr::select(gs_name, gene_symbol) %>% as.data.frame()

# Makes PID geneset for enrichment testing
hPID_df <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:PID")
hPID_t2g <- hPID_df %>% dplyr::select(gs_name, gene_symbol) %>% as.data.frame()

# Make list of genes of interest
humansig <- geneActivity[which(geneActivity$Signature & geneActivity$Biotype == "protein_coding"), ]
mousesig <- geneActivity[which(geneActivity$Mouse_Signature & geneActivity$Biotype == "protein_coding"), ]
conserved <- geneActivity[which(geneActivity$Signature_conserved & geneActivity$Biotype == "protein_coding"), ]
active <- geneActivity[which(geneActivity$humanchondro_all & geneActivity$Biotype == "protein_coding"), ]

DGE_list <- list("HumanSig" = humansig, "MouseSig" = mousesig, "ChondroConserved" = conserved, "Chondrocyte_any" = active)

# Perform pathway enrichment analysis
for(i in names(DGE_list)){
  interesting <- DGE_list[[i]]
  
  BPenrich <- enrichGO(interesting[, "GeneSymbol"],
                       OrgDb = org.Hs.eg.db,
                       keyType = "SYMBOL",
                       ont = "BP",
                       pvalueCutoff = 0.05,
                       pAdjustMethod = "bonferroni",
                       minGSSize = 10,
                       maxGSSize = 500,
                       # readable = TRUE,
                       universe = universe)
  GOBP_GOI <- as.data.frame(BPenrich@result)
  write.csv(GOBP_GOI, paste0("project_results/humanchondro_bulk/Enriched_GOBP_", i, ".csv"))
  
  MFenrich <- enrichGO(interesting[, "GeneSymbol"],
                       OrgDb = org.Hs.eg.db,
                       keyType = "SYMBOL",
                       ont = "MF",
                       pvalueCutoff = 0.05,
                       pAdjustMethod = "bonferroni",
                       minGSSize = 10,
                       maxGSSize = 500,
                       # readable = TRUE,
                       universe = universe)
  GOMF_GOI <- as.data.frame(MFenrich@result)
  write.csv(GOMF_GOI, paste0("project_results/humanchondro_bulk/Enriched_GOMF_", i, ".csv"))
  
  CCenrich <- enrichGO(interesting[, "GeneSymbol"],
                       OrgDb = org.Hs.eg.db,
                       keyType = "SYMBOL",
                       ont = "CC",
                       pvalueCutoff = 0.05,
                       pAdjustMethod = "bonferroni",
                       minGSSize = 10,
                       maxGSSize = 500,
                       # readable = TRUE,
                       universe = universe)
  GOCC_GOI <- as.data.frame(CCenrich@result)
  write.csv(GOCC_GOI, paste0("project_results/humanchondro_bulk/Enriched_GOCC_", i, ".csv"))
  
  # Perform HPO enrichment analysis
  HPOEnrich <- enricher(gene = interesting[, "GeneSymbol"], 
                        TERM2GENE = h_t2g,
                        pvalueCutoff = 1,
                        qvalueCutoff = 1,
                        pAdjustMethod = "bonferroni",
                        minGSSize = 15,
                        maxGSSize = 500,
                        universe = universe)
  HPO_GOI <- as.data.frame(HPOEnrich@result)
  write.csv(HPO_GOI, paste0("project_results/humanchondro_bulk/Enriched_HPO_", i, ".csv"))
  
  TFEnrich <- enricher(gene = interesting[, "GeneSymbol"], 
                       TERM2GENE = htf_t2g,
                       pvalueCutoff = 1,
                       qvalueCutoff = 1,
                       pAdjustMethod = "bonferroni",
                       minGSSize = 15,
                       maxGSSize = 500,
                       universe = universe)
  TF_GOI <- as.data.frame(TFEnrich@result)
  write.csv(TF_GOI, paste0("project_results/humanchondro_bulk/Enriched_TF_", i, ".csv"))
  
  PIDEnrich <- enricher(gene = interesting[, "GeneSymbol"], 
                        TERM2GENE = hPID_t2g,
                        pvalueCutoff = 1,
                        qvalueCutoff = 1,
                        pAdjustMethod = "bonferroni",
                        minGSSize = 15,
                        maxGSSize = 500,
                        universe = universe)
  PID_GOI <- as.data.frame(PIDEnrich@result)
  write.csv(PID_GOI, paste0("project_results/humanchondro_bulk/Enriched_PID_", i, ".csv"))
  
  
  # Perform KEGG enrichment analysis
  KEGGenricMmig <- enrichKEGG(interesting[, "EntrezID"],
                              organism = "hsa",
                              keyType = "kegg",
                              pvalueCutoff = 0.05,
                              pAdjustMethod = "bonferroni",
                              universe = universe_entrez)
  write.csv(KEGGenricMmig, paste0("project_results/humanchondro_bulk/Enriched_KEGG_", i, ".csv"))
  
  # Perform REACTOME enrichment analysis
  REACTOMEenricMmig <- enrichPathway(interesting[, "EntrezID"],
                                     organism = "human",
                                     pvalueCutoff = 0.05,
                                     pAdjustMethod = "bonferroni",
                                     # readable = TRUE,
                                     universe = universe_entrez)
  write.csv(REACTOMEenricMmig, paste0("project_results/humanchondro_bulk/Enriched_REACTOME_", i, ".csv"))
}


