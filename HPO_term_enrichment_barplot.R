#!/usr/bin/Rscript

# --------------------------------------------------------------------------
# Bar plot of Human phenotype data
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

# Read HPO enrichment information
HPO <- read.csv('project_results/humanchondro_bulk/Enriched_HPO_ChondroConserved.csv', header = TRUE)

HPO[1:25, ] %>%
  mutate(Description = factor(Description, levels = rev(Description)),
         LogPvalue = -log10(p.adjust)) %>%
  ggplot(aes(y=Description, x=LogPvalue)) +
  geom_bar(stat = 'identity') +
  ggsave("project_results/humanchondro_bulk/HPO_enrichment_barplot.pdf")


