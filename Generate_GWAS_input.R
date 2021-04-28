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
geneActivity <- read.csv('project_results/growthplate_bulk/Gene_annotation_w_activity.csv', row.names = 1, header = TRUE)
signature <- read.csv('project_results/growthplate_bulk/Chondrocyte_signature_limma.csv', row.names = 1, header = TRUE)

idx <- match(rownames(geneActivity), rownames(signature))
geneActivity$Signature <- signature$Signature [idx]
geneActivity$Universe_signature <- signature$Chondrocyte_any [idx]

# Load data files and calculate mean FPKM in chondrocytes and other tissues
# --------------------------------------------------------------------------

# Mean expression in chondrocytes
chondro_count <- read.csv('data/growthplate_bulk/RSEM_count.csv', row.names = 1, header = TRUE)
rownames(chondro_count) <- sub("(.*?)\\..*", "\\1", rownames(chondro_count))
geneActivity$Log2_cpm_chondrocytes <- rowMeans(cpm(chondro_count, log = TRUE))

# Mean expression in other tissues
organ_count <- read.csv('data/organstissues_bulk/RSEM_count.csv', row.names = 1, header = TRUE)
rownames(organ_count) <- sub("(.*?)\\..*", "\\1", rownames(organ_count))

marrow_count <- read.csv('data/boneymarrow_bulk/RSEM_count.csv', row.names = 1, header = TRUE)
rownames(marrow_count) <- sub("(.*?)\\..*", "\\1", rownames(marrow_count))
marrow_count <- marrow_count[,6:10] # Just marrow samples because bone wasnt used in the signature definition.
combined_count <- cbind(organ_FPKM, marrow_count)
geneActivity$Log2_cpm_other <- rowMeans(cpm(combined_count, log = TRUE))

# Annotate with human symbols
# --------------------------------------------------------------------------

# Subset to gene universe (protein-coding genes expressed in any tissues)
geneActivity <- geneActivity[geneActivity$Biotype == 'protein_coding', ]

# Find human ensembl ids for protein coding genes
ensembl <- useMart(biomart = "ensembl", dataset="mmusculus_gene_ensembl", host = "feb2014.archive.ensembl.org")
human_orths <- getBM(attributes=c('ensembl_gene_id', 'hsapiens_homolog_ensembl_gene', 'hsapiens_homolog_orthology_confidence'),
                     filters = 'ensembl_gene_id',
                     values = rownames(geneActivity),
                     mart = ensembl)

# Annotate mouse ENSEMBL ids with human ENSEMBL Ids
ensembl <- useMart(biomart = "ensembl", dataset="mmusculus_gene_ensembl", host = "feb2021.archive.ensembl.org")
human_orths_new <- getBM(attributes=c('ensembl_gene_id', 'hsapiens_homolog_ensembl_gene', 'hsapiens_homolog_associated_gene_name', 'hsapiens_homolog_orthology_confidence'),
                     filters = 'ensembl_gene_id',
                     values = rownames(geneActivity),
                     mart = ensembl)

# Add human gene symbols to geneActivity dataframe
idx <- match(rownames(geneActivity), human_orths$ensembl_gene_id)
geneActivity$Human_Ensembl_ID_hg19 <- human_orths$hsapiens_homolog_ensembl_gene [idx]

idx <- match(rownames(geneActivity), human_orths_new$ensembl_gene_id)
geneActivity$Human_Ensembl_ID_hg38 <- human_orths_new$hsapiens_homolog_ensembl_gene [idx]
geneActivity$Human_geneSymbol_hg38 <- human_orths_new$hsapiens_homolog_associated_gene_name [idx]

# Remove unwanted columns and write data
geneActivity <- geneActivity[,c(1:4, 64:71)]
write.csv(geneActivity, "project_results/growthplate_bulk/Chondrocyte_signature_hg37_hg38_GWAS_input.csv")
