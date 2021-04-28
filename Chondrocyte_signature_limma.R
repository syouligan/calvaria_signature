#!/usr/bin/Rscript

# --------------------------------------------------------------------------
# Find genes specifically expressed in chondrocytes (chondrocyte transcriptome signature)
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
library("scran")
library("scater")
library('ggcorrplot')

# Load data files (make sure working directory is osteocyte_signature directory) and remove trailing decimal from gene ids
# --------------------------------------------------------------------------

chondro_count <- read.csv('data/growthplate_bulk/RSEM_count.csv', row.names = 1, header = TRUE)
rownames(chondro_count) <- sub("(.*?)\\..*", "\\1", rownames(chondro_count))

organ_count <- read.csv('data/organstissues_bulk/RSEM_count.csv', row.names = 1, header = TRUE)
rownames(organ_count) <- sub("(.*?)\\..*", "\\1", rownames(organ_count))

marrow_count <- read.csv('data/boneymarrow_bulk/RSEM_count.csv', row.names = 1, header = TRUE)
rownames(marrow_count) <- sub("(.*?)\\..*", "\\1", rownames(marrow_count))

combined_counts <- cbind(chondro_count, organ_count, marrow_count)

# Load gene annotation info
geneActivity <- read.csv('project_results/growthplate_bulk/Gene_annotation_w_activity_PZvsHZ_ChondrovsBone.csv', row.names = 1, header = TRUE)

# Process data ready for differential expression analysis
# --------------------------------------------------------------------------

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

# Define experiemental design and contrasts
design <- model.matrix(~0 + sampleType, data=sampleTable)
contrasts <- makeContrasts( vsAdr = sampleTypeChondro-sampleTypeAdr,
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
                            # vsClean = sampleTypeChondro-sampleTypeClean,
                            vsMarrow = sampleTypeChondro-sampleTypeMarrow,
                            levels = design )

# Define signature based on differentially expressed genes in each comparison
# --------------------------------------------------------------------------

# Make DGEList object and normalise using TMM
y <- DGEList(counts = expression, group=sampleType, genes = genemap)
y <- calcNormFactors(y, method="TMM")

# Voom transform, fit linear models and caluclate differential expression
y <- voom(y, design)
fit <- lmFit(y, design)
fit <- contrasts.fit(fit, contrasts)
fit <- treat(fit, lfc = 0)
DEG1 <- decideTests(fit)
summary(DEG1)

# Append LFC and p-value for each comparison to gene Activity table
for (cont in colnames(contrasts)) {
    res <- data.frame(topTreat(fit, n = Inf, p.value = Inf, coef = cont))
    write.csv(res, paste0("project_results/growthplate_bulk/DEGs_", cont, ".csv"))
    
    idx <- match(rownames(geneActivity), rownames(res))
    geneActivity[,paste0(cont, "_LFC")] <- res$logFC [idx]
    geneActivity[,paste0(cont, "_t")] <- res$t [idx]
    geneActivity[,paste0(cont, "_adj.P.Val")] <- res$adj.P.Val [idx]
    geneActivity[,paste0(cont, "_significant")] <- geneActivity[,paste0(cont, "_adj.P.Val")] < 0.05
}

# Find genes significantly differentially upregulated in comprarisons with each tissue
geneActivity$Up_all <- rowSums(geneActivity[grepl("_LFC", colnames(geneActivity))] > 0) == ncol(contrasts)
geneActivity$Mean_LFC <- rowMeans(geneActivity[grepl("_LFC", colnames(geneActivity))])
geneActivity$Significant_all <- rowSums(geneActivity[grepl("_adj.P.Val", colnames(geneActivity))] < 0.05) == ncol(contrasts)
geneActivity$Signature <- geneActivity$Up_all & geneActivity$Significant_all & geneActivity$Biotype == "protein_coding"
geneActivity$Signature[is.na(geneActivity$Signature)] <- FALSE

# Number of genes in signature
sum(geneActivity$Signature)

# Write gene annotation csv with gene activity information
write.csv(geneActivity, "project_results/growthplate_bulk/Chondrocyte_signature_limma.csv")

