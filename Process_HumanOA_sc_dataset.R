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

# Load chondrocyte gene activity data
geneActivity <- read.csv('project_results/humanchondro_bulk/Chondrocyte_signature_limma.csv', row.names = 1, header = TRUE)
geneActivity$Ensembl <- rownames(geneActivity)

# Load counts data
counts <- as(t(read.table("data/humanchondroOA_sc/GSE104782_allcells_UMI_count.txt", header = TRUE, row.names = 1)), "dgCMatrix")

# Load cell information
obs <- read.csv("data/humanchondroOA_sc/GSE104782_Table_Cell_quality_information_and_clustering_information.csv", header = TRUE)
obs$Sample <- gsub( "\\_S[0-9].*", "", obs$Cell)
obs$Grade <- gsub(".*[_]([^.]+)[.].*", "\\1", obs$Cell)
obs$Score <- gsub("S", "", obs$Grade)

# Add gene descriptions
ensembl <- useMart(biomart = "ensembl", dataset="hsapiens_gene_ensembl", host = "may2015.archive.ensembl.org")
human_geneinfo <- getBM(attributes=c('ensembl_gene_id', 'description', 'entrezgene', 'hgnc_symbol'),
                        filters = 'hgnc_symbol',
                        values = colnames(counts),
                        mart = ensembl)

var <- data.frame("Gene_name" = colnames(counts))

idx <- match(var$Gene_name, human_geneinfo$hgnc_symbol)
var$Ensembl <- human_geneinfo$ensembl_gene_id [idx]
var$Description <- human_geneinfo$description [idx]
var$EntrezID <- human_geneinfo$entrezgene [idx]
var$GeneSymbol <- human_geneinfo$hgnc_symbol [idx]

# Merge with geneActivity data
all_rowdata <- data.frame(base::merge(x = var, y = geneActivity, by = "Ensembl", all.x = TRUE, sort = FALSE))

all_rowdata$Signature[is.na(all_rowdata$Signature)] <- FALSE
all_rowdata$Mouse_Signature[is.na(all_rowdata$Mouse_Signature)] <- FALSE
all_rowdata$Signature_conserved[is.na(all_rowdata$Signature_conserved)] <- FALSE
all_rowdata$Nosology2019[is.na(all_rowdata$Nosology2019)] <- FALSE

# Set Col11a2 TRUE for signature manually - issues with mapping from humanchondro_bulk dataset
all_rowdata$Signature[all_rowdata$Ensembl == 'ENSG00000223699'] <- TRUE
all_rowdata$Mouse_Signature[all_rowdata$Ensembl == 'ENSG00000223699'] <- TRUE
all_rowdata$Signature_conserved[all_rowdata$Ensembl == 'ENSG00000223699'] <- TRUE
all_rowdata$Nosology2019[all_rowdata$Ensembl == 'ENSG00000223699'] <- TRUE
rownames(all_rowdata) <- all_rowdata$Gene_name

# Make single cell experiment object and export for python
# --------------------------------------------------------------------------
filtered_exp_sce <- SingleCellExperiment(assays = list('counts' = t(counts)), colData = obs, rowData = all_rowdata)

# Remove filtered cells
filtered_exp_sce <- filtered_exp_sce[,!filtered_exp_sce$Cluster == ""]

# Save dataset components
writeMM(counts(filtered_exp_sce), "data/humanchondroOA_sc/OAsc_filtered_counts.mtx")
write.csv(colData(filtered_exp_sce), "data/humanchondroOA_sc/OAsc_colData.csv", row.names = FALSE)
write.csv(rowData(filtered_exp_sce), "data/humanchondroOA_sc/OAsc_rowData.csv", row.names = TRUE)
