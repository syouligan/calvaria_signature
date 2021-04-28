#!/usr/bin/Rscript

# --------------------------------------------------------------------------
# Import human OA sc dataset to examine gene expression with OA severity
# --------------------------------------------------------------------------

# Set adaptive working directory
if(dir.exists("/Users/mac/cloudstor/")) {
  setwd("/Users/mac/cloudstor/scott_projects/osteochondro_signature/") # Uses practice data (5% of cells from each sample) if running locally
  place <- "local"
} else {
  setwd("/share/ScratchGeneral/scoyou/sarah_projects/osteochondro_signature/")
  place <- "wolfpack"
}

library('DropletUtils')
library('dplyr')
library('tidyr')
library('scater')
library('Matrix')
library('SingleCellExperiment')
library('Seurat')
library('phateR')
library('Rmagic')
library('readr')
library('viridis')
library('biomaRt')

# Load counts data
counts <- as(t(read.table("data/humanchondroOA_sc/GSE104782_allcells_UMI_count.txt", header = TRUE, row.names = 1)), "dgCMatrix")

# Load cell information
obs <- read.csv("data/humanchondroOA_sc/GSE104782_Table_Cell_quality_information_and_clustering_information.csv", header = TRUE)
obs$Sample <- gsub( "\\_S[0-9].*", "", obs$Cell)
obs$Grade <- gsub(".*[_]([^.]+)[.].*", "\\1", obs$Cell)

# Add gene descriptions
ensembl <- useMart(biomart = "ensembl", dataset="hsapiens_gene_ensembl", host = "may2015.archive.ensembl.org")
human_geneinfo <- getBM(attributes=c('ensembl_gene_id', 'description', 'entrezgene', ''),
                        filters = 'ensembl_gene_id',
                        values = gene_annotations$Ensembl,
                        mart = ensembl)


