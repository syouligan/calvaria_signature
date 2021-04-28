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

# Load chondrocyte gene activity data
geneActivity <- read.csv('project_results/humanchondro_bulk/Chondrocyte_signature_limma.csv', row.names = 1, header = TRUE)

# Load REVIGO output (calculatd from GO ID and p-values from enrichGO function, C=0.9, database = Homo sapiens, similarity measure = SimRel)
revigoBP <- read.csv("project_results/humanchondro_bulk/Revigo_GOBP.csv", header = TRUE)

# Load significantly enriched GO terms
BPsig <- read.csv("project_results/humanchondro_bulk/Enriched_GOBP_ChondroConserved.csv", header = TRUE)

# Identify semantically similar clusters of significantly enriched GO terms
# Removal of redundant terms and MDS scaling was performed using http://revigo.irb.hr/ (Supek, Fran, et al. "REVIGO summarizes and visualizes long lists of gene ontology terms." PloS one 6.7 (2011): e21800.)
# --------------------------------------------------------------------------

# Cluster terms based on REVIGO semantic space values and plot BIC, classification and uncertainty plots
revigoBP <- revigoBP[! revigoBP$PlotX == " null",]
semValsBP <- revigoBP[,c("PlotX", "PlotY")]
semValsBP$PlotX <- as.numeric( as.character(semValsBP$PlotX) )
semValsBP$PlotY <- as.numeric( as.character(semValsBP$PlotY) )

clusterBP <- Mclust(semValsBP, modelNames = c("VII", "VVI", "VVV"), G = 1:20)

# Make clustering QC-plots
pdf("project_results/humanchondro_bulk/REVIGO_BP_Signature_Conserved_clusters.pdf")
par(mfrow=c(2,2))
plot(clusterBP, what = "BIC")
plot(clusterBP, what = "classification")
plot(clusterBP, what = "density")
plot(clusterBP, what = "density", type = "persp", col = "dodgerblue3")
dev.off()

# Make dataframe with cluster classification data
classBP <- data.frame(revigoBP)
classBP$classification <- clusterBP$classification
classgenesBP <- classBP[ ,c("TermID", "Name", "classification", "Value", "PlotX", "PlotY")]
idx <- match(classgenesBP$TermID, BPsig$ID)
classgenesBP$Genes <- BPsig$geneID [ idx ]

# Order and save dataframes
classgenesBP <- classgenesBP[order(classgenesBP[,3], classgenesBP[,4]),]
write.csv(classgenesBP, "project_results/humanchondro_bulk/REVIGO_BP_Signature_Conserved_clusters.csv", row.names = FALSE)

# Make dataframe of data for plotting
BP_class <- data.frame(classgenesBP)
BP_class$Name <- as.character(BP_class$Name)
BP_class$PlotX <- as.numeric(as.character(BP_class$PlotX))
BP_class$PlotY <- as.numeric(as.character(BP_class$PlotY))
BP_class$classification <- as.factor(as.character(BP_class$classification))
BP_class$classification <- factor(BP_class$classification, sort(as.numeric(levels(BP_class$classification))))

# Plot clusters with 50% CI and size adjusted to uncertainty
ggplot(data = BP_class) +
  stat_ellipse(data = BP_class, aes( PlotX, PlotY, fill = classification), alpha = 0.3, geom = "polygon", type = "norm", level = 0.5, color = "magenta") +
  geom_point(aes(PlotX, PlotY, fill = classification), shape = 21, colour = "magenta", size = 6, alpha = 0.4) +
  scale_fill_brewer("Cluster", type = "qual", palette = "Set3") +
  scale_size("Certainty", range = c(1, 5), limits = c(0.5, 1)) +
  labs(x = "Semantic space X", y = "Semantic space Y") +
  theme_classic() +
  ggsave("project_results/humanchondro_bulk/BP_Signature_Conserved_cluster_classification.pdf")

# Make bar plots for the top 5 GOBP terms in each cluster ranked by p-value (ascending)
BP_class$Value <- abs(as.numeric(as.character(BP_class$Value)))
all_top <- data.frame()
for (i in rev(levels(BP_class$classification))){
  clust_top <- BP_class[BP_class$classification == i,]
  if(nrow(clust_top) > 2){ 
    clust_top <- clust_top[with(clust_top, order(-Value)),]
    all_top <- rbind(clust_top[1:3,], all_top)
  }
  else {
    print("Too few")
  }
}

all_top$Label <- paste0(all_top$TermID, ": ", all_top$Name)
all_top$Label <- factor(all_top$Label, levels = rev(all_top$Label))

ggplot(data = all_top) +
  geom_col(aes(x = Label, y = Value, fill = classification)) +
  scale_fill_brewer("Cluster", type = "qual", palette = "Set3") +
  theme_classic() +
  labs(x = "-log10 p-value") +
  theme(axis.title.y = element_blank(), axis.text.x=element_text(angle=0,hjust=1)) +
  coord_flip() +
  ggsave("project_results/humanchondro_bulk/ClusterBP_Signature_Conserved_GO_barplot.pdf")

