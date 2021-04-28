#!/usr/bin/Rscript

# --------------------------------------------------------------------------
# Import he_sc dataset to define genes in chondrocytes
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

# Read in gene information
var <- read.csv("data/humanembyonic_sc/GRCh38_genes_var.csv", header = TRUE, row.names = 1)
var$GeneSymbol.unique <- rownames(var$GeneSymbol)

# Load counts data and creat cell info
rm(counts, obs)
for(samp in c("GSM4274193_CS22_2", "GSM4274192_CS22", "GSM4274191_CS20")){
  if(samp == "GSM4274193_CS22_2") {
    counts <- as(readMM(paste0("data/humanembyonic_sc/GSE143753_RAW/", samp, "_longbone_rawdata.mtx")), "dgCMatrix")
    print(nrow(counts))
    obs <- data.frame('Cell_ID' = paste0(samp, "_", 1:nrow(counts)), "Sample" = samp)
  } else {
    tmp <- as(readMM(paste0("data/humanembyonic_sc/GSE143753_RAW/", samp, "_longbone_rawdata.mtx")), "dgCMatrix")
    print(nrow(tmp))
    counts <- rbind(counts, tmp)
    obs <- rbind(obs, data.frame('Cell_ID' = paste0(samp, 1:nrow(tmp)), "Sample" = samp))
  }
}

# Perform differential gene expression
filtered_exp_sce <- SingleCellExperiment(assays = list('counts' = t(counts)), colData = obs, rowData = var)

# Annotate with chondrocyte data
# --------------------------------------------------------------------------

# Load chondrocyte gene activity data
geneActivity <- read.csv('project_results/humanchondro_bulk/Chondrocyte_signature_limma.csv', row.names = 1, header = TRUE)
geneActivity$Ensembl <- rownames(geneActivity)

# Merge gene info and reorder
all_rowdata <- data.frame(base::merge(x = rowData(filtered_exp_sce), y = geneActivity, by = "Ensembl", all.x = TRUE, sort = FALSE))
idx <- match(rowData(filtered_exp_sce)$Ensembl, all_rowdata$Ensembl)
all_rowdata <- all_rowdata[idx, ]
rownames(all_rowdata) <- rownames(filtered_exp_sce)

all_rowdata$Signature[is.na(all_rowdata$Signature)] <- FALSE
all_rowdata$Mouse_Signature[is.na(all_rowdata$Mouse_Signature)] <- FALSE
all_rowdata$Signature_conserved[is.na(all_rowdata$Signature_conserved)] <- FALSE
all_rowdata$Nosology2019[is.na(all_rowdata$Nosology2019)] <- FALSE

rowData(filtered_exp_sce) <- all_rowdata

# Save dataset components
writeMM(counts(filtered_exp_sce), "data/humanembyonic_sc/He_filtered_counts.mtx")
write.csv(colData(filtered_exp_sce), "data/humanembyonic_sc/He_colData.csv", row.names = FALSE)
write.csv(rowData(filtered_exp_sce), "data/humanembyonic_sc/He_rowData.csv", row.names = TRUE)

