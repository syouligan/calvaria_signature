#!/usr/bin/Rscript

# --------------------------------------------------------------------------
#! Perform hypergeometric enrichment on Nosology Skeletal Genetic Disorder Groups
# --------------------------------------------------------------------------

# Working directory
if(dir.exists("/Users/mac/cloudstor/")) {
  setwd("/Users/mac/cloudstor/scott_projects/osteochondro_signature/") # Uses practice data (5% of cells from each sample) if running locally
  place <- "local"
} else {
  setwd("/share/ScratchGeneral/scoyou/sarah_projects/osteochondro_signature/")
  place <- "wolfpack"
}

library(ggplot2)
library(dplyr)
library(forcats)
library(gplots)

# Load and prepare data
# --------------------------------------------------------------------------

# Load chondrocyte gene activity data
geneActivity <- read.csv('project_results/humanchondro_bulk/Chondrocyte_signature_limma.csv', row.names = 1, header = TRUE)

# Load NSGD annotation
NSGD <- read.csv("data/nosology_2019/Nosology_ensembl_annotated.csv", header = TRUE)
NSGD <- NSGD[! (NSGD$Human_Ensembl_ID == "" | is.na(NSGD$Human_Ensembl_ID)), ]
NSGD_key <- NSGD %>%
  dplyr::select(Nosology_ID, Group) %>%
  distinct()

# Overall hypergeometric enrichment in different enriched lists
# --------------------------------------------------------------------------

for(i in c("Signature", "Mouse_Signature", "Signature_conserved")) {
  #Annotate NSGD with those in the human chondrocyte signature
  chondrosig <- geneActivity[which(geneActivity[,i]), ]
  universe <- geneActivity[geneActivity$Biotype == 'protein_coding', ]
  NSGD$Signature <- NSGD$Human_Ensembl_ID %in% rownames(chondrosig) | NSGD$Gene %in% chondrosig$GeneSymbol
  NSGD$Expressed <- NSGD$Human_Ensembl_ID %in% rownames(geneActivity[geneActivity$humanchondro_all, ])  | NSGD$Gene %in% geneActivity[geneActivity$humanchondro_all, 'GeneSymbol']
  
  # Calculate overall enrichment of Nosology genes in the Chondrocyte Signature
  Overlap <- length(unique(NSGD[NSGD$Signature, 'Human_Ensembl_ID']))
  Group1 <- length(unique(NSGD$Human_Ensembl_ID))
  Group2 <- length(unique(rownames(chondrosig)))
  Universe <- length(unique(rownames(universe)))
  Total_enrichment <- phyper(Overlap-1, Group2, Universe-Group2, Group1, lower.tail= FALSE)
  Total_expected <- Group1*Group2/Universe
  Total_foldChange <- Overlap/Total_expected
  Expressed <- length(unique(NSGD[NSGD$Expressed, 'Human_Ensembl_ID']))
  Percent_expressed <- Expressed/Group1*100
  Percent_nosology <- Overlap/Group1*100
  Percent_signature <- Overlap/Group2*100
  
  if( i == "Signature") {
    Overall_enrichment <- data.frame(row.names = c("Overlap", "Group1", "Group2", "Universe", "Total_enrichment", "Total_expected", "Total_foldChange", "Expressed", "Percent_expressed", "Percent_signature", "Percent_nosology"))
    Overall_enrichment[,i] <- c("Overlap" = Overlap, "Group1" = Group1, "Group2" = Group2, "Universe" = Universe, "Total_enrichment" = Total_enrichment, "Total_expected" = Total_expected, "Total_foldChange" = Total_foldChange, "Expressed" = Expressed, "Percent_expressed" = Percent_expressed, "Percent_signature" = Percent_signature, "Percent_nosology" = Percent_nosology)
  } else {
    Overall_enrichment[,i] <- c("Overlap" = Overlap, "Group1" = Group1, "Group2" = Group2, "Universe" = Universe, "Total_enrichment" = Total_enrichment, "Total_expected" = Total_expected, "Total_foldChange" = Total_foldChange, "Expressed" = Expressed, "Percent_expressed" = Percent_expressed, "Percent_signature" = Percent_signature, "Percent_nosology" = Percent_nosology)
  }

  # Groupwise hypergeometric enrichment
  # --------------------------------------------------------------------------
  
  # Calculate hypergeometric enrichment in each group
  HyperEnrich <- NSGD %>%
    dplyr::select(Nosology_ID, Mouse_Ensembl_ID, Signature, Expressed) %>%
    distinct() %>%
    group_by(Nosology_ID) %>%
    summarise(Overlap = sum(Signature), # Number of signature gene in group
              Group1 = n(), # Total number in group
              Group2 = length(unique(rownames(chondrosig))), # Total number in signature
              Universe = length(unique(rownames(universe))),
              Exp = sum(Expressed)) %>% # Total number of protein coding genes in genome
    mutate(Percent_exp = Exp/Group1*100,
           Percent = Overlap/Group1*100,
           Pvalue = phyper(Overlap-1, Group2, Universe-Group2, Group1, lower.tail= FALSE)) %>% # hypergeometric enrichment
    mutate(AdjustedPValue = p.adjust(Pvalue, method = 'BH'),
           Expected = Group1*Group2/Universe) %>%
    mutate(FoldChange = Overlap/Expected)
  
  # Add group name to enrichment results
  idx <- match(HyperEnrich$Nosology_ID, NSGD_key$Nosology_ID)
  HyperEnrich$Group <- NSGD_key$Group [idx]
  
  write.csv(HyperEnrich, paste0("project_results/humanchondro_bulk/Nosology_groupwise_enrichment_", i, ".csv"))
  
  # Make a dotplot of percentage enriched within dysplasia groups
  HyperEnrich %>%
    mutate(Group = fct_reorder(Group, -log10(AdjustedPValue))) %>%
    ggplot() +
    # geom_bar(aes(x = Percent_exp, y = Group), fill = "light blue", stat="identity", width = 0.3) +
    # geom_vline(xintercept = seq(0, 100, 25), linetype = "dashed", color="black") +
    geom_point(aes(x = -log10(AdjustedPValue), y = Group, fill = Percent, size = Group1), shape = 21, color = "black", alpha = 0.5) +
    scale_fill_gradientn(aes(fill = Percent), colors = c("orange", "purple"), limits = c(0, 100), values = NULL, space = "Lab", na.value = "grey80", guide = "colourbar", aesthetics = "fill") +
    scale_size_continuous("Group size", breaks = seq(5, 35, 5), limits = c(1, 35), range = c(1,9)) +
    # scale_x_continuous(breaks = seq(0, 100, 25), labels = seq(0, 100, 25), position = "top", expand = c(0,0.5)) +
    labs(x = "% of Dysplasia Group", x = NA) +
    theme_classic() +
    ggsave(paste0("project_results/humanchondro_bulk/Nosology_signature_percent_dotplot_", i, ".pdf"), useDingbats = FALSE)
  
}

# Plot piecharts of overall enrichment
Overall_enrichment <- data.frame(t(Overall_enrichment))
Overall_enrichment$Percent_not_signature <- 100-Overall_enrichment$Percent_signature
Overall_enrichment$Signatures <- rownames(Overall_enrichment)

Overall_enrichment %>%
  dplyr::select(Signatures, Percent_signature, Percent_not_signature) %>%
  melt.data.frame(id.vars = "Signatures") %>%
  ggplot(aes(x="", y=value, fill=variable)) +
  geom_bar(stat="identity", width=1) +
  coord_polar("y", start=0) +
  scale_fill_manual(values = c("#69b3a2", "Grey90")) +
  facet_wrap(~Signatures, ncol = 3, nrow = 1) +
  theme_classic() +
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        legend.position = 'none') +
  ggsave("project_results/humanchondro_bulk/Nosology_signature_enrichment_piecharts.pdf")

# Calculate the number of individual diseases caused by genes in each signature
chondrosig <- geneActivity[which(geneActivity[,'Signature']), ]
NSGD$Human_Signature <- NSGD$Human_Ensembl_ID %in% rownames(chondrosig) | NSGD$Gene %in% chondrosig$GeneSymbol

chondrosig <- geneActivity[which(geneActivity[,'Mouse_Signature']), ]
NSGD$Mouse_Signature <- NSGD$Human_Ensembl_ID %in% rownames(chondrosig) | NSGD$Gene %in% chondrosig$GeneSymbol

chondrosig <- geneActivity[which(geneActivity[,'Signature_conserved']), ]
NSGD$Signature_conserved <- NSGD$Human_Ensembl_ID %in% rownames(chondrosig) | NSGD$Gene %in% chondrosig$GeneSymbol

Overall_enrichment$Number_disease <- c(sum(NSGD$Human_Signature), sum(NSGD$Mouse_Signature), sum(NSGD$Signature_conserved))
write.csv(Overall_enrichment, "project_results/humanchondro_bulk/Nosology_overall_enrichment_statistics.csv")

