#!/usr/bin/Rscript

# --------------------------------------------------------------------------
#! Plot significant HPO terms
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
library(tidyr)
library(forcats)

# Load and prepare data
# --------------------------------------------------------------------------

# Load chondrocyte gene activity data
geneActivity <- read.csv('project_results/growthplate_bulk/Chondrocyte_signature_limma.csv', row.names = 1, header = TRUE)

# Load HPO enrichment data
HPO <- read.csv("project_results/growthplate_bulk/Enriched_HPO_ChondroSig.csv", header = TRUE)

# Extract number of genes in each group
HPO <- separate(data = HPO, col = BgRatio, into = c("Group", "Total"), sep = "\\/")
HPO$Group <- as.numeric(HPO$Group)
HPO$Percent <- HPO$Count / HPO$Group *100

# Calculate the percent of genes expressed in the HPO group
HPO_any <- read.csv("project_results/growthplate_bulk/Enriched_HPO_Chondrocyte_any.csv", header = TRUE)
HPO_any <- separate(data = HPO_any, col = BgRatio, into = c("Group", "Total"), sep = "\\/")
HPO_any$Group <- as.numeric(HPO_any$Group)
HPO_any$Percent_exp <- HPO_any$Count / HPO_any$Group *100

idx <- match(HPO$ID, HPO_any$ID)
HPO$Percent_exp <- HPO_any$Percent_exp [idx]

# Subset to significantly represented terms
HPO_sig <- HPO[HPO$p.adjust < 0.05,]

HPO_sig %>%
  mutate(Description = fct_reorder(Description, Percent)) %>%
  top_n(50, Percent) %>%
  ggplot() +
  geom_bar(aes(x = Percent_exp, y = Description), fill = "light blue", stat="identity", width = 0.3) +
  geom_vline(xintercept = seq(0, 100, 25), linetype = "dashed", color="black") +
  geom_point(aes(x = Percent, y = Description, fill = -log10(p.adjust), size = Group), shape = 21, color = "white") +
  scale_fill_gradientn(aes(fill = -log10(p.adjust)), colors = c("orange", "purple"), limits = c(1.30103, max(-log10(HPO_sig$p.adjust))), values = NULL, space = "Lab", na.value = "grey80", guide = "colourbar", aesthetics = "fill") +
  scale_size_continuous("Group size", range = c(1,9)) +
  scale_x_continuous(breaks = seq(0, 100, 25), labels = seq(0, 100, 25), position = "top", expand = c(0,0.5)) +
  labs(x = "% of HPO genes", x = NA) +
  theme_classic() +
  ggsave("project_results/growthplate_bulk/HPO_signature_percent_dotplot.pdf", useDingbats = FALSE)


     