#!/usr/bin/Rscript

# --------------------------------------------------------------------------
# Use Baccin dataset to define genes in chondrocytes
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

# Load data on bone marrow niche (https://github.com/veltenlab/rnamagnet), Reference: C. Baccin, J. Al-Sabah, L. Velten, et al.: Combined single-cell and spatial transcriptomics reveals the molecular, cellular and spatial organization of bone marrow niches.
# --------------------------------------------------------------------------

# Load data
load("data/baccin_sc/NicheData10x.rda")
genes <- read.table("data/baccin_sc/Baccin_genes.tsv")
colnames(genes) <- c("Ensembl", "GeneSymbol")
genes$GeneSymbol.unique <- make.unique(genes$GeneSymbol)

# Load chondrocyte gene activity data
geneActivity <- read.csv('project_results/growthplate_bulk/Chondrocyte_signature_limma.csv', row.names = 1, header = TRUE)
geneActivity$Ensembl <- rownames(geneActivity)

# Make sce object
NicheData10x$Cell_type_baccin <- Idents(NicheData10x)
NicheData10x_sce <- as.SingleCellExperiment(NicheData10x)

# Add gene names to row info
idx <- match(rownames(NicheData10x_sce), genes$GeneSymbol.unique)
rowData(NicheData10x_sce)$Ensembl <- genes$Ensembl [idx]
rowData(NicheData10x_sce)$GeneSymbol.unique <- genes$GeneSymbol.unique [idx]

# Add chondrocyte signature data
all_rowdata <- data.frame(base::merge(x = rowData(NicheData10x_sce), y = geneActivity, by = "Ensembl", all.x = TRUE, sort = FALSE))

idx <- match(rownames(NicheData10x_sce), all_rowdata$GeneSymbol.unique)
all_rowdata <- all_rowdata[idx, ]
rownames(all_rowdata) <- rownames(NicheData10x_sce)
all_rowdata$Signature[is.na(all_rowdata$Signature)] <- FALSE
rowData(NicheData10x_sce) <- all_rowdata

# Save dataset components
writeMM(counts(NicheData10x_sce), "data/baccin_sc/Baccin_filtered_counts.mtx")
write.csv(colData(NicheData10x_sce), "data/baccin_sc/Baccin_colData.csv")
write.csv(rowData(NicheData10x_sce), "data/baccin_sc/Baccin_rowData.csv", row.names = TRUE)
