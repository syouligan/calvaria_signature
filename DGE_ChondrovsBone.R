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

chondro_count <- read.csv('data/growthplate_bulk/RSEM_count.csv', row.names = 1, header = TRUE)
rownames(chondro_count) <- sub("(.*?)\\..*", "\\1", rownames(chondro_count))

marrow_count <- read.csv('data/boneymarrow_bulk/RSEM_count.csv', row.names = 1, header = TRUE)
rownames(marrow_count) <- sub("(.*?)\\..*", "\\1", rownames(marrow_count))

combined_counts <- cbind(chondro_count, marrow_count)

# Load gene annotation info
geneActivity <- read.csv('project_results/growthplate_bulk/Gene_annotation_w_activity_PZvsHZ.csv', row.names = 1, header = TRUE)

# Subset to active genes
genemap <- geneActivity[geneActivity$HZ_all | geneActivity$PZ_all | geneActivity$Clean_activity == 5, ]
expression <- combined_counts[rownames(combined_counts) %in% rownames(genemap),]

# Make a design table
sampleName <- colnames(expression)
sampleType <- factor(c(gsub( "\\_[0-9]*$", "", sampleName)))
sampleRep <- factor(c(gsub( "^.*?_","", sampleName)))
sampleTissue <- ifelse(grepl("Z", sampleName), "Chondrocyte", "Bone")
sampleTable <- data.frame(sampleName, sampleType, sampleRep, sampleTissue)

# Make DGEList object and normalise using TMM
# Find differentially expressed genes with TMM normalisation
design <- model.matrix(~sampleTissue, data = sampleTable)
y <- DGEList(counts = expression, group=sampleTissue, genes = genemap)
y <- calcNormFactors(y, method="TMM")
y <- estimateGLMCommonDisp(y, design)
y <- estimateGLMTagwiseDisp(y, design)
fit <- glmFit(y, design)
lrt <- glmLRT(fit, coef=2)
topTags(lrt)
res <- data.frame(topTags(lrt, n = nrow(expression)))

idx <- match(rownames(geneActivity), rownames(res))
geneActivity[,"ChondrovsBone-LFC"] <- res$logFC [idx]
geneActivity[,"ChondrovsBone-LR"] <- res$LR [idx]
geneActivity[,"ChondrovsBone-FDR"] <- res$FDR [idx]
geneActivity[,"ChondrovsBone-significant"] <- geneActivity[,"ChondrovsBone-FDR"] < 0.05

# Write gene annotation csv with gene activity information
write.csv(geneActivity, "project_results/growthplate_bulk/Gene_annotation_w_activity_PZvsHZ_ChondrovsBone.csv")
