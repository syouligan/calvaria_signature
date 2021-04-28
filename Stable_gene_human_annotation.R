#!/usr/bin/Rscript

# --------------------------------------------------------------------------
#! Annotate signature with human gene ids (hg19 and hg38)
# --------------------------------------------------------------------------

# Working directory
if(dir.exists("/Users/mac/cloudstor/")) {
  setwd("/Users/mac/cloudstor/scott_projects/osteochondro_signature/") # Uses practice data (5% of cells from each sample) if running locally
  place <- "local"
} else {
  setwd("/share/ScratchGeneral/scoyou/sarah_projects/osteochondro_signature/")
  place <- "wolfpack"
}

library("biomaRt")
library("openxlsx")
library("org.Hs.eg.db")
library('edgeR')

# Load and prepare data
# --------------------------------------------------------------------------

# Load chondrocyte gene activity data
Gene_annotation <- read.csv('data/growthplate_bulk/Gene_annotation.csv', header = TRUE, row.names = 1)
Gene_annotation <- Gene_annotation[,-c(5, 6, 7)]

# Annotate with human symbols
# --------------------------------------------------------------------------

# Subset to gene universe (protein-coding genes expressed in any tissues)
genes <- rownames(Gene_annotation[Gene_annotation$Biotype == 'protein_coding', ])

# Find human ensembl ids for protein coding genes
ensembl <- useMart(biomart = "ensembl", dataset="mmusculus_gene_ensembl", host = "feb2014.archive.ensembl.org")
human_orths <- getBM(attributes=c('ensembl_gene_id', 'hsapiens_homolog_ensembl_gene', 'hsapiens_homolog_orthology_confidence'),
                     filters = 'ensembl_gene_id',
                     values = genes,
                     mart = ensembl)

# Annotate mouse ENSEMBL ids with human ENSEMBL Ids
ensembl <- useMart(biomart = "ensembl", dataset="mmusculus_gene_ensembl", host = "feb2021.archive.ensembl.org")
human_orths_new <- getBM(attributes=c('ensembl_gene_id', 'hsapiens_homolog_ensembl_gene', 'hsapiens_homolog_associated_gene_name', 'hsapiens_homolog_orthology_confidence'),
                     filters = 'ensembl_gene_id',
                     values = genes,
                     mart = ensembl)

# Add human gene symbols to Gene_annotation dataframe
idx <- match(rownames(Gene_annotation), human_orths$ensembl_gene_id)
Gene_annotation$Human_Ensembl_ID_hg19 <- human_orths$hsapiens_homolog_ensembl_gene [idx]

idx <- match(rownames(Gene_annotation), human_orths_new$ensembl_gene_id)
Gene_annotation$Human_Ensembl_ID_hg38 <- human_orths_new$hsapiens_homolog_ensembl_gene [idx]
Gene_annotation$Human_geneSymbol_hg38 <- human_orths_new$hsapiens_homolog_associated_gene_name [idx]

# Annotate with human entrez ids
ensembl <- useMart(biomart = "ensembl", dataset="hsapiens_gene_ensembl", host = "feb2021.archive.ensembl.org")
human_entrez <- getBM(attributes=c('ensembl_gene_id', 'entrezgene_id'),
                         filters = 'ensembl_gene_id',
                         values = Gene_annotation$Human_Ensembl_ID_hg38,
                         mart = ensembl)

idx <- match(Gene_annotation$Human_Ensembl_ID_hg38, human_entrez$ensembl_gene_id)
Gene_annotation$Human_ENTREZ_ID <- human_entrez$entrezgene_id [idx]

# Remove unwanted columns and write data
write.csv(Gene_annotation, "data/growthplate_bulk/Gene_annotation_human.csv")
