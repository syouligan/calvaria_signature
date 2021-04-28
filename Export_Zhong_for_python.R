#!/usr/bin/Rscript

# --------------------------------------------------------------------------
# Use zhong_sc dataset to define genes in chondrocytes
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

# Load data on Zhong L, Yao L, Tower RJ, Wei Y et al. Single cell transcriptomics identifies a unique adipose lineage cell population that regulates bone marrow environment. Elife 2020 Apr 14;9. PMID: 32286228
# --------------------------------------------------------------------------

# Load data
genes <- read.table("data/zhong_sc/Zhong_genes.tsv")
colnames(genes) <- c("Ensembl", "GeneSymbol")
genes$GeneSymbol.unique <- make.unique(genes$GeneSymbol)

cell_mat_1 <- Matrix::Matrix(as.matrix(read.table("data/zhong_sc/GSM4318799_1M_matrix.txt", sep = ' ')), sparse = TRUE)
cell_mat_1.5 <- Matrix::Matrix(as.matrix(read.table("data/zhong_sc/GSM4318800_1.5M_matrix.txt", sep = ' ')), sparse = TRUE)
cell_mat_3 <- Matrix::Matrix(as.matrix(read.table("data/zhong_sc/GSM4318801_3M_matrix.txt", sep = ' ')), sparse = TRUE)
cell_mat_16 <- Matrix::Matrix(as.matrix(read.table("data/zhong_sc/GSM4318802_16M_matrix.txt", sep = ' ')), sparse = TRUE)

# Find and filter to common genes between all datasets
common_genes <- Reduce(intersect, list(rownames(cell_mat_1), rownames(cell_mat_1.5), rownames(cell_mat_3), rownames(cell_mat_16)))

cell_mat_1 <- cell_mat_1[rownames(cell_mat_1) %in% common_genes,]
cell_mat_1.5 <- cell_mat_1.5[rownames(cell_mat_1.5) %in% common_genes,]
cell_mat_3 <- cell_mat_3[rownames(cell_mat_3) %in% common_genes,]
cell_mat_16 <- cell_mat_16[rownames(cell_mat_16) %in% common_genes,]

# Create unified data matrix
cell_mat <- cbind(cell_mat_1, cell_mat_1.5, cell_mat_3, cell_mat_16)

# Annotate with chondrocyte data
# --------------------------------------------------------------------------

# Load chondrocyte gene activity data
geneActivity <- read.csv('project_results/humanchondro_bulk/Chondrocyte_signature_limma.csv', row.names = 1, header = TRUE)
geneActivity$Ensembl <- rownames(geneActivity)

# Make sce object
zhong_sce <- SingleCellExperiment(assays=list(counts=cell_mat))

# Add gene names to row info
idx <- match(rownames(zhong_sce), genes$GeneSymbol.unique)
rowData(zhong_sce)$Ensembl <- genes$Ensembl [idx]
rowData(zhong_sce)$GeneSymbol.unique <- genes$GeneSymbol.unique [idx]

# Add chondrocyte signature data
all_rowdata <- data.frame(base::merge(x = rowData(zhong_sce), y = geneActivity, by.x= 'Ensembl', by.y = "MouseEnsemblID", all.x = TRUE, sort = FALSE))

idx <- match(rownames(zhong_sce), all_rowdata$GeneSymbol.unique)
all_rowdata <- all_rowdata[idx, ]
rownames(all_rowdata) <- rownames(zhong_sce)
all_rowdata$Signature[is.na(all_rowdata$Signature)] <- FALSE
all_rowdata$Mouse_Signature[is.na(all_rowdata$Mouse_Signature)] <- FALSE
all_rowdata$Signature_conserved[is.na(all_rowdata$Signature_conserved)] <- FALSE
all_rowdata$Nosology2019[is.na(all_rowdata$Nosology2019)] <- FALSE

rowData(zhong_sce) <- all_rowdata

# Remove tomato genes
zhong_sce <- zhong_sce[! rownames(zhong_sce) %in% 'tomato',]

# Add vector of ages to colData
zhong_sce$Ages <- c(rep("1", ncol(cell_mat_1)), rep("1.5", ncol(cell_mat_1.5)), rep("3", ncol(cell_mat_3)), rep("16", ncol(cell_mat_16)))

# Save dataset components
writeMM(counts(zhong_sce), "data/zhong_sc/Zhong_filtered_counts.mtx")
write.csv(colData(zhong_sce), "data/zhong_sc/Zhong_colData.csv")
write.csv(rowData(zhong_sce), "data/zhong_sc/Zhong_rowData.csv", row.names = TRUE)
