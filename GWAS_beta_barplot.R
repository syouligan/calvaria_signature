#!/usr/bin/Rscript

# --------------------------------------------------------------------------
# Jenson-shannon analysis of single cell data
# --------------------------------------------------------------------------

# Set adaptive working directory
if(dir.exists("/Users/mac/cloudstor/")) {
  setwd("/Users/mac/cloudstor/scott_projects/osteochondro_signature/") # Uses practice data (5% of cells from each sample) if running locally
  place <- "local"
} else {
  setwd("/share/ScratchGeneral/scoyou/sarah_projects/osteochondro_signature/")
  place <- "wolfpack"
}

library(ggridges)
library(ggplot2)
library(reshape2)

# Load and prepare data
# --------------------------------------------------------------------------

# Read gene annotation information
GWAS_magma <- read.csv('project_results/GWAS/MAGMA/GWAS_MAGMA_geneset_results.csv', header = TRUE)
GWAS_magma$ymin <- GWAS_magma$Beta - GWAS_magma$SE
GWAS_magma$ymax <- GWAS_magma$Beta + GWAS_magma$SE

# Correct p-value
GWAS_magma$p.adjust <- p.adjust(GWAS_magma$P.value, method = "BH")
GWAS_magma$LogPvalue <- -log10(GWAS_magma$p.adjust)

# Make factors
GWAS_magma$Trait <- factor(GWAS_magma$Trait, levels = rev(unique(GWAS_magma$Trait)))
GWAS_magma$Signature <- factor(GWAS_magma$Signature, levels = rev(c("ChondroConSig", "ChondroHumSig", "ChondroMusSig")))

# Plot Beta per trait
ggplot(GWAS_magma, aes(x = Trait, y = Beta, group = Signature)) +
  geom_bar(stat = 'identity', aes(fill = Signature), position = 'dodge') +
  geom_errorbar(aes(ymin = ymin, ymax = ymax), position = 'dodge') +
  geom_hline(yintercept = 0, linetype = 'dashed') +
  coord_flip() +
  theme_classic() +
  theme(axis.line.y = element_blank()) +
  ggsave("project_results/GWAS/Trait_GSenrichment_Beta_plot.pdf")



