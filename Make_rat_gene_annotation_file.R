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
gene_annotations <- import.gff('data/ratchondro_bulk/Rattus_norvegicus.Rnor_6.0.103.gtf.gz')
gene_annotations <- gene_annotations[gene_annotations$type == 'gene',]
gene_annotations <- data.frame('Ensembl' = gene_annotations$gene_id, 'GeneSymbol' = gene_annotations$gene_name, 'Biotype' = gene_annotations$gene_biotype)

# Add gene descriptions and mouse ensembl IDs
ensembl <- useMart(biomart = "ensembl", dataset="rnorvegicus_gene_ensembl")
rat_geneinfo <- getBM(attributes=c('ensembl_gene_id', 'description', 'mmusculus_homolog_ensembl_gene', 'mmusculus_homolog_associated_gene_name', 'hsapiens_homolog_ensembl_gene', 'hsapiens_homolog_associated_gene_name'),
                        filters = 'ensembl_gene_id',
                        values = gene_annotations$Ensembl,
                        mart = ensembl)

idx <- match(gene_annotations$Ensembl, rat_geneinfo$ensembl_gene_id)
gene_annotations$Description <- rat_geneinfo$description [idx]
gene_annotations$MouseEnsemblID <- rat_geneinfo$mmusculus_homolog_ensembl_gene [idx]
gene_annotations$MouseGeneSymbol <- rat_geneinfo$mmusculus_homolog_associated_gene_name [idx]
gene_annotations$HumanEnsemblID <- rat_geneinfo$hsapiens_homolog_ensembl_gene [idx]
gene_annotations$HumanGeneSymbol <- rat_geneinfo$hsapiens_homolog_associated_gene_name [idx]

# Add entrez ID
rat_geneinfo <- getBM(attributes=c('ensembl_gene_id', 'entrezgene_id'),
                        filters = 'ensembl_gene_id',
                        values = gene_annotations$Ensembl,
                        mart = ensembl)

idx <- match(gene_annotations$Ensembl, rat_geneinfo$ensembl_gene_id)
gene_annotations$EntrezID <- rat_geneinfo$entrezgene_id [idx]

# Add Nosology of genetic skeletal disorders info
nosology <- read.csv("~/cloudstor/scott_projects/osteochondro_signature/data/nosology_2019/Nosology_ensembl_annotated.csv", header = TRUE)
gene_annotations$Nosology2019 <- gene_annotations$HumanEnsemblID %in% nosology$Human_Ensembl_ID
gene_annotations$Nosology2019[gene_annotations$HumanGeneSymbol %in% c('IKBKG', 'RBM8A', 'FAM58A')] <- TRUE # Genes with new ENSEMBL IDs
gene_annotations$Nosology2019[gene_annotations$HumanEnsemblID == ""] <- FALSE # Genes with new ENSEMBL IDs

# Reorder and save
gene_annotations <- gene_annotations %>%
  arrange(Ensembl) %>%
  select(Ensembl, GeneSymbol, Description, EntrezID, Biotype, Nosology2019, MouseEnsemblID, MouseGeneSymbol, HumanEnsemblID, HumanGeneSymbol)

write.csv(gene_annotations, "data/ratchondro_bulk/Gene_annotation.csv", row.names = FALSE)

