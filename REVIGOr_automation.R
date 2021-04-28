#!/usr/bin/Rscript

# --------------------------------------------------------------------------
#! Perform patwhay enrichment in the chondrocyte signature
# --------------------------------------------------------------------------

# Working directory
if(dir.exists("/Users/mac/cloudstor/")) {
  setwd("/Users/mac/cloudstor/scott_projects/osteochondro_signature/") # Uses practice data (5% of cells from each sample) if running locally
  place <- "local"
} else {
  setwd("/share/ScratchGeneral/scoyou/sarah_projects/osteochondro_signature/")
  place <- "wolfpack"
}

library("clusterProfiler")
library("org.Mm.eg.db")
library("mclust")
library("tm")
library("ggplot2")
library("revigoR")
library("reticulate")

# Load chondrocyte gene activity data
geneActivity <- read.csv('project_results/growthplate_bulk/Chondrocyte_signature_limma.csv', row.names = 1, header = TRUE)

# Load significantly enriched GO terms
BPsig <- read.csv("project_results/growthplate_bulk/Enriched_GOBP_ChondroSig.csv", header = TRUE)
rownames(BPsig) <- BPsig$ID
input <- BPsig[BPsig$p.adjust < 0.05, "p.adjust", drop = FALSE]

# Submit significant GO terms and scrape output from revigo
dir.create("project_results/growthplate_bulk/REVIGO")
reticulate::virtualenv_create('revigo')
reticulate::virtualenv_install(envname = 'revigo',
                               packages = c('numpy', 'werkzeug', 'robobrowser'),
                               ignore_installed = TRUE)
scrape_revigo("project_results/growthplate_bulk/REVIGO", input)




