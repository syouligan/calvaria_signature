#!/usr/bin/Rscript

# --------------------------------------------------------------------------
# Find genes specifically expressed in chondrocytes
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
library("scran")
library("scater")

# Load data files (make sure working directory is osteocyte_signature directory) and remove trailing decimal from gene ids
chondro_count <- read.csv('data/growthplate_bulk/RSEM_count.csv', row.names = 1, header = TRUE)
rownames(chondro_count) <- sub("(.*?)\\..*", "\\1", rownames(chondro_count))

organ_count <- read.csv('data/organstissues_bulk/RSEM_count.csv', row.names = 1, header = TRUE)
rownames(organ_count) <- sub("(.*?)\\..*", "\\1", rownames(organ_count))

combined_counts <- cbind(chondro_count, organ_count)

# Load gene annotation info
geneActivity <- read.csv('project_results/growthplate_bulk/Gene_annotation_w_activity.csv', row.names = 1, header = TRUE)

# Subset to active genes
genemap <- geneActivity[geneActivity$HZ_all | geneActivity$PZ_all, ]
expression <- combined_counts[rownames(combined_counts) %in% rownames(genemap),]

# Make a design table
sampleName <- colnames(expression)
sampleType <- c(gsub( "\\_[0-9]*$", "", sampleName))
sampleType[grepl("PZ|HZ", sampleType)] <- 'Chondro'
sampleType <- factor(sampleType)
sampleRep <- factor(c(gsub( "^.*?_","", sampleName)))
sampleChondro <- ifelse(grepl("PZ|HZ", sampleName), "Chondro", "Other")
sampleType[grepl("PZ|HZ", sampleType)] <- 'Chondro'
sampleTable <- data.frame(sampleName, sampleType, sampleRep, sampleChondro)

# Normalise counts based on House keeping genes to merge with other organs
data("segList_ensemblGeneID", package = "scMerge")

set <- newSeqExpressionSet(as.matrix(round(expression)), phenoData = data.frame(sampleTable, row.names=colnames(expression)))
plotPCA(set, cex=1)
set <- betweenLaneNormalization(set, which="upper")
plotPCA(set, cex=1)
hks <- segList_ensemblGeneID$mouse$bulkRNAseqHK[segList_ensemblGeneID$mouse$bulkRNAseqHK %in% rownames(genemap)]
set1 <- RUVg(set, hks, k=1)
plotPCA(set1, cex=1)

# Make DGEList object and normalise using TMM
# Find differentially expressed genes with TMM normalisation
design <- model.matrix(~0 + sampleType, data=pData(set1))
CONTRASTS <- makeContrasts( vsAdr = sampleTypeChondro-sampleTypeAdr,
                            vsAor = sampleTypeChondro-sampleTypeAor,
                            vsBFat = sampleTypeChondro-sampleTypeBFat,
                            vsBstm = sampleTypeChondro-sampleTypeBstm,
                            vsCer = sampleTypeChondro-sampleTypeCer,
                            vsHrt = sampleTypeChondro-sampleTypeHrt,
                            vsHyp = sampleTypeChondro-sampleTypeHyp,
                            vsKid = sampleTypeChondro-sampleTypeKid,
                            vsLiv = sampleTypeChondro-sampleTypeLiv,
                            vsLun = sampleTypeChondro-sampleTypeLun,
                            vsMus = sampleTypeChondro-sampleTypeMus,
                            vsWFat = sampleTypeChondro-sampleTypeWFat,
                            levels = design )
y <- DGEList(counts = expression, group=sampleType, genes = genemap)
y <- calcNormFactors(y, method="TMM")
y <- estimateGLMCommonDisp(y, design)
y <- estimateGLMTagwiseDisp(y, design)
fit <- glmFit(y, design)
lrt <- glmLRT(fit, contrast = CONTRASTS)
topTags(lrt)
res <- data.frame(topTags(lrt, n = nrow(expression)))

res$Significant <- res$FDR < 0.05
res$Mean_LFC <- rowMeans(res[,grepl("logFC", colnames(res))])
res$All_up <- rowSums(res[,grepl("logFC", colnames(res))] > 0) == 12
res$Signature <- res$All_up & res$Significant
sum(res$Signature)

# Write gene annotation csv with gene activity information
write.csv(res, "project_results/growthplate_bulk/Chondrocyte_signature_edgeR.csv")

