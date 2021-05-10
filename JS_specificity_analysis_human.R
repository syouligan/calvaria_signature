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
geneActivity <- read.csv('project_results/humanchondro_bulk/Chondrocyte_signature_limma.csv', header = TRUE)

# Read JS score data for human ESC scRNAseq (1.0 resolution clusters)
JS_1 <- read.csv("~/cloudstor/scott_projects/osteochondro_signature/project_results/humanembyonic_sc/JS_specificity_scores_Cell_types.csv", header = TRUE)
JS_1 <- JS_1[,c(13:16,2:12)]
JS_1$Signature_conserved <- ifelse(JS_1$Signature_conserved == "True", TRUE, FALSE)

# Subset to protein coding genes
JS_1 <- JS_1[JS_1$Biotype == 'protein_coding', ]

# Plot overall distribution of JS specificity scores
# --------------------------------------------------------------------------

# Melt into long format
JS_1_long <- melt(JS_1[,1:15], id.vars = c('GeneSymbol', 'Ensembl', 'Biotype', 'Signature_conserved'))

# Breakdown gene expression based on specificity score (> 95 percentile of JS scores is specific)
js_threshold <- quantile(JS_1_long$value, c(0.95))
JS_1$Above_JS_threshold <- rowSums(JS_1[,5:15] > js_threshold)

ggplot(JS_1_long, aes(x = value)) +
  geom_density(fill = "#69b3a2") +
  theme_classic() +
  geom_vline(xintercept = js_threshold, linetype = 'dashed') +
  ggsave("project_results/humanembyonic_sc/JS_specificity_distribution_95th_percentile.pdf")

# Make ridge plot of JS score for each tissue
ggplot(JS_1_long, aes(x = value, y = variable)) +
  geom_density_ridges(fill = "#69b3a2") +
  geom_vline(xintercept = js_threshold, linetype = 'dashed') +
  ggsave("project_results/humanembyonic_sc/JS_specificity_ridges_all_genes_95th_percentile.pdf")


ggplot(JS_1_long[JS_1_long$Signature_conserved,], aes(x = value, y = variable)) +
  geom_density_ridges(fill = "#69b3a2") +
  geom_vline(xintercept = js_threshold, linetype = 'dashed') +
  ggsave("project_results/humanembyonic_sc/JS_specificity_ridges_signature_genes_95th_percentile.pdf")

# Classify genes based on specificity scores
JS_1 <- JS_1 %>%
  mutate(Specificity = case_when(
    Above_JS_threshold > 1 ~ "Moderate specificity",
    Above_JS_threshold == 1 ~ "High specificity",
    Above_JS_threshold < 1 ~ "Non-specific")
  )

# Make bar charts of cell type of origin for highly specific genes
# --------------------------------------------------------------------------

# Number of genes per cell types with high specifity
JS_sig_high <- JS_1 %>%
  filter(Signature_conserved & Specificity == "High specificity")

JS_sig_genes <- data.frame('High' = colSums(JS_sig_high[,5:15] > js_threshold))

# Number of genes per cell types with moderate specifity (note by definition, each gene will be present in more than 1 group thus total > 101)
JS_sig_mod <- JS_1 %>%
  filter(Signature_conserved & Specificity == "Moderate specificity")

JS_sig_genes$Moderate <- colSums(JS_sig_mod[,5:15] > js_threshold)

# Number of genes with any specificity
JS_sig_genes$Cell_types <- factor(rownames(JS_sig_genes), levels = rownames(JS_sig_genes)[order(JS_sig_genes$High + JS_sig_genes$Moderate)])

# Make a bar plot of number of genes per group
JS_sig_genes %>%
  melt(id.vars = "Cell_types") %>%
  mutate(variable = factor(variable, levels = c('Moderate', 'High'))) %>%
  ggplot(aes(x = Cell_types, y = value, fill = variable)) +
  labs(x = "Number of genes", y = "Cell types") +
  geom_bar(stat = 'identity') +
  scale_fill_manual(values = c("#FFD99F", "#7D1D67")) +
  coord_flip() +
  theme_classic() +
  theme(aspect.ratio = 1) +
  ggsave("project_results/humanembyonic_sc/JS_specificity_breakdown_Signature_conserved_genes_barplots.pdf")

# Make heatmap of JS specificity scores for all genes and signature genes
# --------------------------------------------------------------------------

# Plot specificity of all genes genes
JS_sig <- JS_1

hc_order <- JS_sig$GeneSymbol[order(JS_sig$Chondrogenic)]
JS_sig_long <- melt(JS_sig[, c(1, 5:15)])
colnames(JS_sig_long) <- c('GeneSymbol', 'Cell_type', 'Specificity')
JS_sig_long$Cell_type <- factor(JS_sig_long$Cell_type, levels = rev(rownames(JS_sig_genes)[order(JS_sig_genes$High + JS_sig_genes$Moderate)]))

JS_sig_long %>%
  mutate(GeneSymbol = factor(GeneSymbol, levels = hc_order)) %>%
  ggplot(aes(x=Cell_type, y=GeneSymbol, fill = Specificity)) + 
  geom_tile() +
  theme_classic() +
  # coord_fixed() +
  theme(axis.text.y = element_text(color = 'black', face = 'italic'),
        axis.text.x = element_text(color = 'black', hjust = 1, angle = 45),
        axis.line=element_blank(),
        axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank()) +
  scale_fill_gradientn(colours = hcl.colors(7, 'Temps', rev = FALSE), breaks = seq(0, max(JS_sig_long$Specificity), 0.1), limits = c(min(JS_sig_long$Specificity), max(JS_sig_long$Specificity))) +
  ggsave("project_results/humanembyonic_sc/JS_specificity_heatmap_all_genes.pdf")

# Plot specificity of signature genes
JS_sig <- JS_1 %>%
  filter(Signature_conserved)

hc_order <- JS_sig$GeneSymbol[order(JS_sig$Chondrogenic)]
JS_sig_long <- melt(JS_sig[, c(1, 5:15)])
colnames(JS_sig_long) <- c('GeneSymbol', 'Cell_type', 'Specificity')
JS_sig_long$Cell_type <- factor(JS_sig_long$Cell_type, levels = rev(rownames(JS_sig_genes)[order(JS_sig_genes$High + JS_sig_genes$Moderate)]))

JS_sig_long %>%
  mutate(GeneSymbol = factor(GeneSymbol, levels = hc_order)) %>%
  ggplot(aes(x=Cell_type, y=GeneSymbol, fill = Specificity)) + 
  geom_tile() +
  theme_classic() +
  # coord_fixed() +
  theme(axis.text.y = element_text(color = 'black', face = 'italic'),
        axis.text.x = element_text(color = 'black', hjust = 1, angle = 45),
        axis.line=element_blank(),
        axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank()) +
  scale_fill_gradientn(colours = hcl.colors(7, 'Temps', rev = FALSE), breaks = seq(0, max(JS_sig_long$Specificity), 0.1), limits = c(min(JS_sig_long$Specificity), max(JS_sig_long$Specificity))) +
  ggsave("project_results/humanembyonic_sc/JS_specificity_heatmap_signature_conserved_genes.pdf")
         
# Make pie charts of specificity for all genes and just those in the conserved signature
# --------------------------------------------------------------------------

JS_1 %>%
  dplyr::select(GeneSymbol, Signature_conserved, Above_JS_threshold, Specificity) %>%
  distinct() %>%
  group_by(Specificity) %>%
  summarise(Number = n()) %>%
  mutate(Percent_exp = Number/sum(Number)*100,
         Total = sum(Number)) %>%
  ggplot(aes(x="", y=Percent_exp, fill=Specificity)) +
  geom_bar(stat="identity", width=1) +
  coord_polar("y", start=0) +
  scale_fill_manual(values = c("#7D1D67", "#FFD99F", "Grey90")) +
  theme_classic() +
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank()) +
  ggsave("project_results/humanembyonic_sc/JS_specificity_breakdown_all_genes.pdf")

JS_1 %>%
  dplyr::select(GeneSymbol, Signature_conserved, Above_JS_threshold, Specificity) %>%
  distinct() %>%
  filter(Signature_conserved) %>%
  group_by(Specificity) %>%
  summarise(Number = n()) %>%
  mutate(Percent_exp = Number/sum(Number)*100,
         Total = sum(Number)) %>%
  ggplot(aes(x="", y=Percent_exp, fill=Specificity)) +
  geom_bar(stat="identity", width=1) +
  coord_polar("y", start=0) +
  scale_fill_manual(values = c("#7D1D67", "#FFD99F", "Grey90")) +
  theme_classic() +
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank()) +
  ggsave("project_results/humanembyonic_sc/JS_specificity_breakdown_Signature_conserved_genes.pdf")

