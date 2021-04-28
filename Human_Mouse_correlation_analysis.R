#!/usr/bin/Rscript

# --------------------------------------------------------------------------
# Correlation analysis between human and mouse chondrocytes
# --------------------------------------------------------------------------

# Set adaptive working directory
if(dir.exists("/Users/mac/cloudstor/")) {
  setwd("/Users/mac/cloudstor/scott_projects/osteochondro_signature/") # Uses practice data (5% of cells from each sample) if running locally
  place <- "local"
} else {
  setwd("/share/ScratchGeneral/scoyou/sarah_projects/osteochondro_signature/")
  place <- "wolfpack"
}

library(dplyr)
library(reshape2)
library(ggplot2)
library(ggpubr)
library(gplots)
library(msigdbr)
library(org.Hs.eg.db)
library(ggrepel)

# Load gene annotation info
geneActivity_mouse <- read.csv('project_results/growthplate_bulk/Chondrocyte_signature_limma.csv', header = TRUE)
geneActivity_rat <- read.csv('project_results/ratchondro_bulk/Gene_annotation_w_activity.csv', header = TRUE)
geneActivity_human <- read.csv('project_results/humanchondro_bulk/Chondrocyte_signature_limma.csv', header = TRUE)

# Subset to actively expressed protein coding genes
human_active <- geneActivity_human[geneActivity_human$humanchondro_all & geneActivity_human$Biotype == "protein_coding", ] # human
mouse_active <- geneActivity_human[geneActivity_human$Mouse_chondrocyte_active & geneActivity_human$Biotype == "protein_coding", ] # mouse
rat_active <- geneActivity_human[geneActivity_human$Rat_chondrocyte_active & geneActivity_human$Biotype == "protein_coding", ] # rat

# Compare gene activity in human and mouse chondrocytes
# --------------------------------------------------------------------------

# Number of genes active in humans without mouse or rat orthologs
sum(human_active$MouseEnsemblID == "")
sum(human_active$RatEnsemblID == "")

# Number of protein coding genes active in chondrocytes
DGE_list <- list("Human" = human_active[human_active$MouseEnsemblID != "" & human_active$RatEnsemblID != "", "X"], "Mouse" = mouse_active[mouse_active$MouseEnsemblID != "" & mouse_active$RatEnsemblID != "", "X"], "Rat" = rat_active[rat_active$MouseEnsemblID != "" & rat_active$RatEnsemblID != "", "X"])

pdf("project_results/humanchondro_bulk/Human_Mouse_Rat_activity_venn.pdf")
ItemsList <- venn(DGE_list, show.plot = TRUE)
dev.off()

# Overall correlation in protein-coding genes 
# --------------------------------------------------------------------------

hm_coactive <- geneActivity_human[geneActivity_human$Biotype == "protein_coding" & geneActivity_human$MouseEnsemblID != "" & geneActivity_human$RatEnsemblID != "", ]

COR <- round(cor(hm_coactive$log2_humanchondro_mean_FPKM, hm_coactive$log2_mousechondro_mean_FPKM, method = 'sp', use = "complete.obs"), 2)

# Calculate correlation and plot relationship
ggplot(hm_coactive, aes(x=log2_humanchondro_mean_FPKM, y=log2_mousechondro_mean_FPKM)) +
  geom_point() +
  stat_cor(method = "spearman", cor.coef.name = "rho", alternative = "two.sided") +
  geom_smooth(method = 'lm', se = TRUE) +
  theme_classic() +
  theme(aspect.ratio = 1) +
  ggsave("project_results/humanchondro_bulk/correlations/Mouse_human_correlation_overall.pdf")

COR <- round(cor(hm_coactive$log2_humanchondro_mean_FPKM, hm_coactive$log2_ratchondro_mean_FPKM, method = 'sp', use = "complete.obs"), 2)

# Calculate correlation and plot relationship
ggplot(hm_coactive, aes(x=log2_humanchondro_mean_FPKM, y=log2_ratchondro_mean_FPKM)) +
  geom_point() +
  stat_cor(method = "spearman", cor.coef.name = "rho", alternative = "two.sided") +
  geom_smooth(method = 'lm', se = TRUE) +
  theme_classic() +
  theme(aspect.ratio = 1) +
  ggsave("project_results/humanchondro_bulk/correlations/Rat_human_correlation_overall.pdf")

COR <- round(cor(hm_coactive$log2_mousechondro_mean_FPKM, hm_coactive$log2_ratchondro_mean_FPKM, method = 'sp', use = "complete.obs"), 2)

# Calculate correlation and plot relationship
ggplot(hm_coactive, aes(x=log2_mousechondro_mean_FPKM, y=log2_ratchondro_mean_FPKM)) +
  geom_point() +
  stat_cor(method = "spearman", cor.coef.name = "rho", alternative = "two.sided") +
  geom_smooth(method = 'lm', se = TRUE) +
  theme_classic() +
  theme(aspect.ratio = 1) +
  ggsave("project_results/humanchondro_bulk/correlations/Mouse_rat_correlation_overall.pdf")

# Correlation across key GO biological processes 
# --------------------------------------------------------------------------

# Makes GO genesets for correlation testing
hGOBP_df <- msigdbr(species = "Homo sapiens", category = "C5", subcategory = "GO:BP")

hGOBP_t2g <- hGOBP_df %>% 
  dplyr::select(gs_exact_source, gs_name, gene_symbol) %>%
  group_split(gs_exact_source)

# Caluculate correlation for significant GO terms
rm(cor_df_GOBP)
for (term in 1:length(hGOBP_t2g)) {
  tmp <- hGOBP_t2g[[term]]
  genes <- tmp$gene_symbol
  tmp_2 <- hm_coactive[which(hm_coactive$GeneSymbol %in% genes), ]
  
  if(length(unique(tmp_2$GeneSymbol)) > 25) {
    
    # Calculate spearman correlation between human and mouse
    COR_humanmouse <- round(cor(tmp_2$log2_humanchondro_mean_FPKM, tmp_2$log2_mousechondro_mean_FPKM, method = 'sp', use = "complete.obs"), 2)
    COR_humanrat <- round(cor(tmp_2$log2_humanchondro_mean_FPKM, tmp_2$log2_ratchondro_mean_FPKM, method = 'sp', use = "complete.obs"), 2)
    COR_mouserat <- round(cor(tmp_2$log2_mousechondro_mean_FPKM, tmp_2$log2_ratchondro_mean_FPKM, method = 'sp', use = "complete.obs"), 2)
    
    if(exists("cor_df_GOBP")){
      cor_df_GOBP <- rbind(cor_df_GOBP, data.frame('ID' = unique(tmp$gs_exact_source), 'Term' = unique(tmp$gs_name), 'Genes' = nrow(tmp_2), 'rho_humanmouse' = COR_humanmouse, 'rho_humanrat' = COR_humanrat, 'rho_mouserat' = COR_mouserat))
    } else {
      cor_df_GOBP <- data.frame('ID' = unique(tmp$gs_exact_source), 'Term' = unique(tmp$gs_name), 'Genes' = nrow(tmp_2), 'rho_humanmouse' = COR_humanmouse, 'rho_humanrat' = COR_humanrat, 'rho_mouserat' = COR_mouserat)
    }
    
    # Plot spearman correlation between human and rodents
    if(COR_humanmouse > 0.6) {
      ggplot(tmp_2, aes(x=log2_humanchondro_mean_FPKM, y=log2_mousechondro_mean_FPKM)) +
        geom_text(aes(label = GeneSymbol)) +
        stat_cor(method = "spearman", cor.coef.name = "rho", alternative = "two.sided") +
        geom_smooth(method = 'lm', se = TRUE) +
        theme_classic() +
        theme(aspect.ratio = 1) +
        ggsave(paste0("project_results/humanchondro_bulk/correlations/GOBP/", COR_humanmouse, "_Mouse_human_correlation_", unique(tmp$gs_name), ".pdf"))
    } else {
      print("Low correlation mouse")
    }

    if(COR_humanrat > 0.6) {
      ggplot(tmp_2, aes(x=log2_humanchondro_mean_FPKM, y=log2_ratchondro_mean_FPKM)) +
        geom_text(aes(label = GeneSymbol)) +
        stat_cor(method = "spearman", cor.coef.name = "rho", alternative = "two.sided") +
        geom_smooth(method = 'lm', se = TRUE) +
        theme_classic() +
        theme(aspect.ratio = 1) +
        ggsave(paste0("project_results/humanchondro_bulk/correlations/GOBP/", COR_humanrat, "_Rat_human_correlation_", unique(tmp$gs_name), ".pdf"))
    } else {
      print("Low correlation rat")
    }
    
        
  } else {
    print("Too few genes")
  }
}

# Make labels for key terms
key_terms <- c("GO_CARTILAGE_DEVELOPMENT_INVOLVED_IN_ENDOCHONDRAL_BONE_MORPHOGENESIS", "GO_BONE_GROWTH", "GO_CHONDROCYTE_DIFFERENTIATION")
cor_df_GOBP$KEY <- ""
cor_df_GOBP$key_terms <- cor_df_GOBP$Term %in% key_terms
ix_label <- match(key_terms, cor_df_GOBP$Term)
cor_df_GOBP$KEY[ix_label] <- cor_df_GOBP$Term[ix_label]

# Plot correlation distribution labelling key terms
cor_df_GOBP %>%
  # arrange(rho) %>%
  # mutate(Term = factor(KEY, levels = KEY)) %>%
  ggplot(aes(x = rho_humanmouse, y = rho_humanrat, label = KEY)) +
  geom_point(aes(color = key_terms), alpha = 0.2) +
  scale_colour_manual(values = c("Black", "Red")) +
  geom_text_repel() +
  theme(aspect.ratio = 1) +
  theme_classic() +
  geom_hline(yintercept = c(0.5), linetype = 'dashed') +
  geom_vline(xintercept = c(0.5), linetype = 'dashed') +
  theme(axis.line=element_blank(), aspect.ratio = 1) +
  ggsave("project_results/humanchondro_bulk/correlations/Correlation_distribution_GOBP_rat_mouse.pdf")

# Correlation across key HPO terms
# --------------------------------------------------------------------------

# Makes GO genesets for correlation testing
hHPO_df <- msigdbr(species = "Homo sapiens", category = "C5", subcategory = "HPO")

hHPO_t2g <- hHPO_df %>% 
  dplyr::select(gs_exact_source, gs_name, gene_symbol) %>%
  group_split(gs_exact_source)

# Caluculate correlation for significant HPO terms
rm(cor_df_HPO)
for (term in 1:length(hHPO_t2g)) {
  tmp <- hHPO_t2g[[term]]
  genes <- tmp$gene_symbol
  tmp_2 <- hm_coactive[which(hm_coactive$GeneSymbol %in% genes), ]
  
  if(length(unique(tmp_2$GeneSymbol)) > 25) {
    
    # Calculate spearman correlation between human and mouse
    COR_humanmouse <- round(cor(tmp_2$log2_humanchondro_mean_FPKM, tmp_2$log2_mousechondro_mean_FPKM, method = 'sp', use = "complete.obs"), 2)
    COR_humanrat <- round(cor(tmp_2$log2_humanchondro_mean_FPKM, tmp_2$log2_ratchondro_mean_FPKM, method = 'sp', use = "complete.obs"), 2)
    COR_mouserat <- round(cor(tmp_2$log2_mousechondro_mean_FPKM, tmp_2$log2_ratchondro_mean_FPKM, method = 'sp', use = "complete.obs"), 2)
    
    if(exists("cor_df_HPO")){
      cor_df_HPO <- rbind(cor_df_HPO, data.frame('ID' = unique(tmp$gs_exact_source), 'Term' = unique(tmp$gs_name), 'Genes' = nrow(tmp_2), 'rho_humanmouse' = COR_humanmouse, 'rho_humanrat' = COR_humanrat, 'rho_mouserat' = COR_mouserat))
    } else {
      cor_df_HPO <- data.frame('ID' = unique(tmp$gs_exact_source), 'Term' = unique(tmp$gs_name), 'Genes' = nrow(tmp_2), 'rho_humanmouse' = COR_humanmouse, 'rho_humanrat' = COR_humanrat, 'rho_mouserat' = COR_mouserat)
    }
    
    # Plot spearman correlation between human and rodents
    if(COR_humanmouse > 0.6) {
      ggplot(tmp_2, aes(x=log2_humanchondro_mean_FPKM, y=log2_mousechondro_mean_FPKM)) +
        geom_text(aes(label = GeneSymbol)) +
        stat_cor(method = "spearman", cor.coef.name = "rho", alternative = "two.sided") +
        geom_smooth(method = 'lm', se = TRUE) +
        theme_classic() +
        theme(aspect.ratio = 1) +
        ggsave(paste0("project_results/humanchondro_bulk/correlations/HPO/", COR_humanmouse, "_Mouse_human_correlation_", unique(tmp$gs_name), ".pdf"))
    } else {
      print("Low correlation mouse")
    }
    
    if(COR_humanrat > 0.6) {
      ggplot(tmp_2, aes(x=log2_humanchondro_mean_FPKM, y=log2_ratchondro_mean_FPKM)) +
        geom_text(aes(label = GeneSymbol)) +
        stat_cor(method = "spearman", cor.coef.name = "rho", alternative = "two.sided") +
        geom_smooth(method = 'lm', se = TRUE) +
        theme_classic() +
        theme(aspect.ratio = 1) +
        ggsave(paste0("project_results/humanchondro_bulk/correlations/HPO/", COR_humanrat, "_Rat_human_correlation_", unique(tmp$gs_name), ".pdf"))
    } else {
      print("Low correlation rat")
    }
    
    
  } else {
    print("Too few genes")
  }
}

# Make labels for key terms
key_terms <- c("HP_ABNORMALITY_OF_UPPER_LIMB_EPIPHYSIS_MORPHOLOGY", "HP_SHORT_THORAX", "HP_OSTEOARTHRITIS")
cor_df_HPO$KEY <- ""
cor_df_HPO$key_terms <- cor_df_HPO$Term %in% key_terms
ix_label <- match(key_terms, cor_df_HPO$Term)
cor_df_HPO$KEY[ix_label] <- cor_df_HPO$Term[ix_label]

# Plot correlation distribution labelling key terms
cor_df_HPO %>%
  # arrange(rho) %>%
  # mutate(Term = factor(KEY, levels = KEY)) %>%
  ggplot(aes(x = rho_humanmouse, y = rho_humanrat, label = KEY)) +
  geom_point(aes(color = key_terms), alpha = 0.2) +
  scale_colour_manual(values = c("Black", "Red")) +
  geom_text_repel() +
  theme(aspect.ratio = 1) +
  theme_classic() +
  geom_hline(yintercept = c(0.5), linetype = 'dashed') +
  geom_vline(xintercept = c(0.5), linetype = 'dashed') +
  theme(axis.line=element_blank(), aspect.ratio = 1) +
  ggsave("project_results/humanchondro_bulk/correlations/Correlation_distribution_HPO_rat_mouse.pdf")

# Correlation across REACTOME pathways
# --------------------------------------------------------------------------

# Makes GO genesets for correlation testing
hREACTOME_df <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:REACTOME")

hREACTOME_t2g <- hREACTOME_df %>% 
  dplyr::select(gs_exact_source, gs_name, gene_symbol) %>%
  group_split(gs_exact_source)

# Caluculate correlation for significant REACTOME terms
rm(cor_df_REACTOME)
for (term in 1:length(hREACTOME_t2g)) {
  tmp <- hREACTOME_t2g[[term]]
  genes <- tmp$gene_symbol
  tmp_2 <- hm_coactive[which(hm_coactive$GeneSymbol %in% genes), ]
  
  if(length(unique(tmp_2$GeneSymbol)) > 25) {
    
    # Calculate spearman correlation between human and mouse
    COR_humanmouse <- round(cor(tmp_2$log2_humanchondro_mean_FPKM, tmp_2$log2_mousechondro_mean_FPKM, method = 'sp', use = "complete.obs"), 2)
    COR_humanrat <- round(cor(tmp_2$log2_humanchondro_mean_FPKM, tmp_2$log2_ratchondro_mean_FPKM, method = 'sp', use = "complete.obs"), 2)
    COR_mouserat <- round(cor(tmp_2$log2_mousechondro_mean_FPKM, tmp_2$log2_ratchondro_mean_FPKM, method = 'sp', use = "complete.obs"), 2)
    
    if(exists("cor_df_REACTOME")){
      cor_df_REACTOME <- rbind(cor_df_REACTOME, data.frame('ID' = unique(tmp$gs_exact_source), 'Term' = unique(tmp$gs_name), 'Genes' = nrow(tmp_2), 'rho_humanmouse' = COR_humanmouse, 'rho_humanrat' = COR_humanrat, 'rho_mouserat' = COR_mouserat))
    } else {
      cor_df_REACTOME <- data.frame('ID' = unique(tmp$gs_exact_source), 'Term' = unique(tmp$gs_name), 'Genes' = nrow(tmp_2), 'rho_humanmouse' = COR_humanmouse, 'rho_humanrat' = COR_humanrat, 'rho_mouserat' = COR_mouserat)
    }
    
    # Plot spearman correlation between human and rodents
    if(COR_humanmouse > 0.6) {
      ggplot(tmp_2, aes(x=log2_humanchondro_mean_FPKM, y=log2_mousechondro_mean_FPKM)) +
        geom_text(aes(label = GeneSymbol)) +
        stat_cor(method = "spearman", cor.coef.name = "rho", alternative = "two.sided") +
        geom_smooth(method = 'lm', se = TRUE) +
        theme_classic() +
        theme(aspect.ratio = 1) +
        ggsave(paste0("project_results/humanchondro_bulk/correlations/REACTOME/", COR_humanmouse, "_Mouse_human_correlation_", unique(tmp$gs_name), ".pdf"))
    } else {
      print("Low correlation mouse")
    }
    
    if(COR_humanrat > 0.6) {
      ggplot(tmp_2, aes(x=log2_humanchondro_mean_FPKM, y=log2_ratchondro_mean_FPKM)) +
        geom_text(aes(label = GeneSymbol)) +
        stat_cor(method = "spearman", cor.coef.name = "rho", alternative = "two.sided") +
        geom_smooth(method = 'lm', se = TRUE) +
        theme_classic() +
        theme(aspect.ratio = 1) +
        ggsave(paste0("project_results/humanchondro_bulk/correlations/REACTOME/", COR_humanrat, "_Rat_human_correlation_", unique(tmp$gs_name), ".pdf"))
    } else {
      print("Low correlation rat")
    }
    
    
  } else {
    print("Too few genes")
  }
}

# Make labels for key terms
key_terms <- c("REACTOME_ECM_PROTEOGLYCANS", "REACTOME_COLLAGEN_BIOSYNTHESIS_AND_MODIFYING_ENZYMES", "REACTOME_SIGNALING_BY_FGFR")
cor_df_REACTOME$KEY <- ""
cor_df_REACTOME$key_terms <- cor_df_REACTOME$Term %in% key_terms
ix_label <- match(key_terms, cor_df_REACTOME$Term)
cor_df_REACTOME$KEY[ix_label] <- cor_df_REACTOME$Term[ix_label]

# Plot correlation distribution labelling key terms
cor_df_REACTOME %>%
  # arrange(rho) %>%
  # mutate(Term = factor(KEY, levels = KEY)) %>%
  ggplot(aes(x = rho_humanmouse, y = rho_humanrat, label = KEY)) +
  geom_point(aes(color = key_terms), alpha = 0.2) +
  scale_colour_manual(values = c("Black", "Red")) +
  geom_text_repel() +
  theme(aspect.ratio = 1) +
  theme_classic() +
  geom_hline(yintercept = c(0.5), linetype = 'dashed') +
  geom_vline(xintercept = c(0.5), linetype = 'dashed') +
  theme(axis.line=element_blank(), aspect.ratio = 1) +
  ggsave("project_results/humanchondro_bulk/correlations/Correlation_distribution_REACTOME_rat_mouse.pdf")

