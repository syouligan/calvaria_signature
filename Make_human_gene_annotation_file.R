#!/usr/bin/Rscript

# --------------------------------------------------------------------------
# Create human gene annotation file
# --------------------------------------------------------------------------

# Set adaptive working directory
if(dir.exists("/Users/mac/cloudstor/")) {
  setwd("/Users/mac/cloudstor/scott_projects/osteochondro_signature/") # Uses practice data (5% of cells from each sample) if running locally
  place <- "local"
} else {
  setwd("/share/ScratchGeneral/scoyou/sarah_projects/osteochondro_signature/")
  place <- "wolfpack"
}

library('rtracklayer')
library('biomaRt')
library('dplyr')

# Load gene annotation info
gene_annotations <- import.gff('data/humanchondro_bulk/Homo_sapiens.GRCh38.79.gtf.gz')
gene_annotations <- gene_annotations[gene_annotations$type == 'gene',]
gene_annotations <- data.frame('Ensembl' = gene_annotations$gene_id, 'GeneSymbol' = gene_annotations$gene_name, 'Biotype' = gene_annotations$gene_biotype)

# Add gene descriptions and mouse ensembl IDs
ensembl <- useMart(biomart = "ensembl", dataset="hsapiens_gene_ensembl", host = "may2015.archive.ensembl.org")
human_geneinfo <- getBM(attributes=c('ensembl_gene_id', 'description', 'mmusculus_homolog_ensembl_gene', 'rnorvegicus_homolog_ensembl_gene'),
                     filters = 'ensembl_gene_id',
                     values = gene_annotations$Ensembl,
                     mart = ensembl)

idx <- match(gene_annotations$Ensembl, human_geneinfo$ensembl_gene_id)
gene_annotations$Description <- human_geneinfo$description [idx]
gene_annotations$MouseEnsemblID <- human_geneinfo$mmusculus_homolog_ensembl_gene [idx]
gene_annotations$RatEnsemblID <- human_geneinfo$rnorvegicus_homolog_ensembl_gene [idx]

# Add entrez ID
human_geneinfo <- getBM(attributes=c('ensembl_gene_id', 'entrezgene'),
                        filters = 'ensembl_gene_id',
                        values = gene_annotations$Ensembl,
                        mart = ensembl)

idx <- match(gene_annotations$Ensembl, human_geneinfo$ensembl_gene_id)
gene_annotations$EntrezID <- human_geneinfo$entrezgene [idx]

# Add Nosology of genetic skeletal disorders info
nosology <- read.csv("~/cloudstor/scott_projects/osteochondro_signature/data/nosology_2019/Nosology_ensembl_annotated.csv", header = TRUE)
gene_annotations$Nosology2019 <- gene_annotations$Ensembl %in% nosology$Human_Ensembl_ID[!nosology$Human_Ensembl_ID == ""] | gene_annotations$GeneSymbol %in% nosology$Gene[!nosology$Gene == ""]

# Reorder and save
gene_annotations <- gene_annotations %>%
  arrange(Ensembl) %>%
  dplyr::select(Ensembl, GeneSymbol, Description, EntrezID, Biotype, Nosology2019, MouseEnsemblID, RatEnsemblID)

write.csv(gene_annotations, "data/humanchondro_bulk/Gene_annotation.csv", row.names = FALSE)

