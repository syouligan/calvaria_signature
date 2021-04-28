#!/usr/bin/Rscript

# --------------------------------------------------------------------------
# Find genes differentially expressed in hypertrophic and proliferative chondrocytes 
# --------------------------------------------------------------------------

# Set adaptive working directory
if(dir.exists("/Users/mac/cloudstor/")) {
  setwd("/Users/mac/cloudstor/scott_projects/osteochondro_signature/") # Uses practice data (5% of cells from each sample) if running locally
  place <- "local"
} else {
  setwd("/share/ScratchGeneral/scoyou/sarah_projects/osteochondro_signature/")
  place <- "wolfpack"
}

library("limma")
library("edgeR")
library("ggplot2")
library("reshape")
library("RUVSeq")
library("scMerge")

# Load data files (make sure working directory is osteocyte_signature directory) and remove trailing decimal from gene ids
chondro_count <- read.csv('data/growthplate_bulk/RSEM_count.csv', row.names = 1, header = TRUE)
rownames(chondro_count) <- sub("(.*?)\\..*", "\\1", rownames(chondro_count))

# Load gene annotation info
geneActivity <- read.csv('project_results/growthplate_bulk/Gene_annotation_w_activity.csv', row.names = 1, header = TRUE)

# Subset to active genes
genemap <- geneActivity[geneActivity$HZ_all | geneActivity$PZ_all, ]
expression <- chondro_count[rownames(chondro_count) %in% rownames(genemap),]

# Make a design table
sampleName <- colnames(expression)
sampleType <- factor(c(gsub( "\\_[0-9]*$", "", sampleName)))
sampleRep <- factor(c(gsub( "^.*?_","", sampleName)))
sampleBone <- ifelse(grepl("Phalanx", sampleName), "Phalanx", "Tibia")
sampleWeek <- ifelse(grepl("4wk", sampleName), "4wk", "1wk")
sampleZone <- ifelse(grepl("PZ", sampleName), "PZ", "HZ")
sampleTable <- data.frame(sampleName, sampleType, sampleRep, sampleBone, sampleWeek, sampleZone)

# Make DGEList object and normalise using TMM
# Find differentially expressed genes with TMM normalisation
design <- model.matrix(~sampleZone + sampleWeek + sampleBone, data = sampleTable)
y <- DGEList(counts = expression, group=sampleZone, genes = genemap)
y <- calcNormFactors(y, method="TMM")
y <- estimateGLMCommonDisp(y, design)
y <- estimateGLMTagwiseDisp(y, design)
fit <- glmFit(y, design)
lrt <- glmLRT(fit, coef=2)
topTags(lrt)
res <- data.frame(topTags(lrt, n = nrow(expression)))

idx <- match(rownames(geneActivity), rownames(res))
geneActivity[,"PZvsHZ-LFC"] <- res$logFC [idx]
geneActivity[,"PZvsHZ-LR"] <- res$LR [idx]
geneActivity[,"PZvsHZ-FDR"] <- res$FDR [idx]
geneActivity[,"PZvsHZ-significant"] <- geneActivity[,"PZvsHZ-FDR"] < 0.05

# Write gene annotation csv with gene activity information
write.csv(geneActivity, "project_results/growthplate_bulk/Gene_annotation_w_activity_PZvsHZ.csv")

