#!/usr/bin/Rscript

# --------------------------------------------------------------------------
# Calculate human chondrocyte signature
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
library('biomaRt')

# Load gene annotation info
geneActivity <- read.csv('project_results/humanchondro_bulk/Gene_annotation_w_activity.csv', row.names = 1, header = TRUE)

# Combined FPKM values
combined_counts <- read.csv('project_results/humanchondro_bulk/Humanchondro_organs_FPKM_genes.csv', row.names = 1, header = TRUE)

# Process data ready for differential expression analysis
# --------------------------------------------------------------------------

# Subset to active genes
genemap <- geneActivity[geneActivity$humanchondro_all, ]
expression <- combined_counts[rownames(combined_counts) %in% rownames(genemap),]
genemap <- genemap[rownames(expression), ]

# Make a design table
sampleName <- colnames(expression)
sampleType <- c(gsub( "\\_[0-9]*$", "", sampleName))
sampleType <- factor(sampleType)
sampleRep <- factor(c(gsub( "^.*?_","", sampleName)))
sampleTable <- data.frame(sampleName, sampleType, sampleRep)

# Define experiemental design and contrasts
design <- model.matrix(~0 + sampleType, data=sampleTable)
contrasts <- makeContrasts( vsAdrenal.gland = sampleTypeCartilage-sampleTypeAdrenal.gland,
                            vsBrain = sampleTypeCartilage-sampleTypeBrain,
                            vsCD4.primary.cells = sampleTypeCartilage-sampleTypeCD4.primary.cells,
                            vsCD8.primary.cells = sampleTypeCartilage-sampleTypeCD8.primary.cells,
                            vsHeart = sampleTypeCartilage-sampleTypeHeart,
                            vsKidney = sampleTypeCartilage-sampleTypeKidney,
                            vsLarge.intestine = sampleTypeCartilage-sampleTypeLarge.intestine,
                            vsLiver = sampleTypeCartilage-sampleTypeLiver,
                            vsLung = sampleTypeCartilage-sampleTypeLung,
                            vsMuscle = sampleTypeCartilage-sampleTypeMuscle,
                            vsOvary.PGCs = sampleTypeCartilage-sampleTypeOvary.PGCs,
                            vsPlacenta = sampleTypeCartilage-sampleTypePlacenta,
                            vsRenal.cortex = sampleTypeCartilage-sampleTypeRenal.cortex,
                            vsRenal.pelvis = sampleTypeCartilage-sampleTypeRenal.pelvis,
                            vsSmall.intestine = sampleTypeCartilage-sampleTypeSmall.intestine,
                            vsSpinal.cord = sampleTypeCartilage-sampleTypeSpinal.cord,
                            vsSpleen = sampleTypeCartilage-sampleTypeSpleen,
                            vsStomach = sampleTypeCartilage-sampleTypeStomach,
                            vsTestis.PGCs = sampleTypeCartilage-sampleTypeTestis.PGCs,
                            vsThymus = sampleTypeCartilage-sampleTypeThymus,
                            levels = design )

# Define signature based on differentially expressed genes in each comparison
# --------------------------------------------------------------------------

# Use limma trend to calculate differentially expressed genes
y <- log2(combined_counts + 0.1)
fit <- lmFit(y, design)
fit <- contrasts.fit(fit, contrasts)
fit <- eBayes(fit, robust=TRUE, trend=TRUE)
DEG1 <- decideTests(fit)
summary(DEG1)

# Append LFC and p-value for each comparison to gene Activity table
for (cont in colnames(contrasts)) {
  res <- data.frame(topTable(fit, n = Inf, p.value = Inf, coef = cont))
  write.csv(res, paste0("project_results/humanchondro_bulk/DEGs_", cont, ".csv"))
  
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
geneActivity$Signature[rownames(geneActivity ) == "ENSG00000204248"] <- TRUE # Need to add Col11a2 manually due to error in transition from REFSeq IDs
geneActivity$Signature[rownames(geneActivity ) == "ENSG00000223699"] <- FALSE

# Number of genes in signature
sum(geneActivity$Signature)

# Annotate with mouse signature information
# --------------------------------------------------------------------------

geneActivity_mouse <- read.csv('project_results/growthplate_bulk/Chondrocyte_signature_limma.csv', row.names = 1, header = TRUE)

# Calcullate mean expression in mouse chondrocytes
geneActivity_mouse$mousechondro_mean_FPKM <- rowMeans(geneActivity_mouse[,grepl("wk_mean_FPKM", colnames(geneActivity_mouse))])

# Annotate human data with mouse information
idx <- match(geneActivity$MouseEnsemblID, rownames(geneActivity_mouse))
geneActivity$MouseGeneSymbol <- geneActivity_mouse$GeneSymbol [idx]
geneActivity$Mouse_chondrocyte_active <- geneActivity_mouse$Chondrocyte_any [idx]
geneActivity$mousechondro_mean_FPKM <- geneActivity_mouse$mousechondro_mean_FPKM [idx]
geneActivity$Mouse_Signature <- geneActivity_mouse$Signature [idx]
geneActivity$Mouse_Mean_LFC <- geneActivity_mouse$Mean_LFC [idx]
geneActivity$Skeletal_GO_Term <- geneActivity_mouse$Skeletal_GO_Term [idx]

geneActivity$Signature_conserved <- geneActivity$Signature & geneActivity$Mouse_Signature
geneActivity$Signature_conserved[is.na(geneActivity$Signature_conserved)] <- FALSE

sum(geneActivity$Signature_conserved)

# Venn of Signature overlap
humansig <- geneActivity[which(geneActivity[,'Signature']), ]
mousesig <- geneActivity[which(geneActivity[,'Mouse_Signature']), ]

DGE_list <- list("Human" = rownames(humansig), "Mouse" = rownames(mousesig))

pdf("project_results/humanchondro_bulk/Signature_conserved_venn.pdf")
ItemsList <- venn(DGE_list, show.plot = TRUE)
dev.off()

# Annotate with rat signature information
# --------------------------------------------------------------------------

geneActivity_rat <- read.csv('project_results/ratchondro_bulk/Gene_annotation_w_activity.csv', row.names = 1, header = TRUE)

idx <- match(geneActivity$RatEnsemblID, rownames(geneActivity_rat))
geneActivity$Rat_chondrocyte_active <- geneActivity_rat$Chondrocyte_any [idx]
geneActivity$ratchondro_mean_FPKM <- geneActivity_rat$ratchondro_mean_FPKM [idx]

# Log transform mean expression data for interspecies comparison
geneActivity$log2_ratchondro_mean_FPKM <- log2(geneActivity$ratchondro_mean_FPKM + 0.1)
geneActivity$log2_humanchondro_mean_FPKM <- log2(geneActivity$humanchondro_mean_FPKM + 0.1)
geneActivity$log2_mousechondro_mean_FPKM <- log2(geneActivity$mousechondro_mean_FPKM + 0.1)

# Write gene annotation csv with gene activity information
write.csv(geneActivity, "project_results/humanchondro_bulk/Chondrocyte_signature_limma.csv")


