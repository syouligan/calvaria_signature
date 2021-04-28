#!/usr/bin/Rscript

# --------------------------------------------------------------------------
# Make gene expression table for human chondro and organ data
# --------------------------------------------------------------------------

# Set adaptive working directory
if(dir.exists("/Users/mac/cloudstor/")) {
  setwd("/Users/mac/cloudstor/scott_projects/osteochondro_signature/") # Uses practice data (5% of cells from each sample) if running locally
  place <- "local"
} else {
  setwd("/share/ScratchGeneral/scoyou/sarah_projects/osteochondro_signature/")
  place <- "wolfpack"
}

library("limma")
library("edgeR")
library("ggplot2")
library("reshape")
library('biomaRt')

# Make dataframe of counts for human organs and chondrocytes
# --------------------------------------------------------------------------

# Add ensembl IDs
combined_counts <- read.csv("data/humanchondro_bulk/Human_chondro_organ_FPKM.csv", header = TRUE)

ensembl <- useMart(biomart = "ensembl", dataset="hsapiens_gene_ensembl", host = "may2015.archive.ensembl.org")
human_geneinfo <- getBM(attributes=c('ensembl_gene_id', 'refseq_mrna'),
                        filters = 'refseq_mrna',
                        values = combined_counts$Refseq.id,
                        mart = ensembl)

idx <- match(combined_counts$Refseq.id, human_geneinfo$refseq_mrna)
combined_counts$Ensembl <- human_geneinfo$ensembl_gene_id [idx]
combined_counts <- combined_counts[,3:114]

# Keep transcript with highest expression in chondrocytes
combined_counts <- data.frame(combined_counts %>%
  group_by(Ensembl) %>%
  summarize_all(sum))
  
combined_counts <- combined_counts[! (combined_counts$Ensembl == "" | is.na(combined_counts$Ensembl)), ]

write.csv(combined_counts, 'project_results/humanchondro_bulk/Humanchondro_organs_FPKM_genes.csv', row.names = FALSE)





