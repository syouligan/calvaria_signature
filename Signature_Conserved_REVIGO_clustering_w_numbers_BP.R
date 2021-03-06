#!/usr/bin/Rscript

# --------------------------------------------------------------------------
#! Make I-graph from cytoscape (REVIGO - 0.7 allowed similarity)
# --------------------------------------------------------------------------

# Set adaptive working directory
if(dir.exists("/Users/mac/cloudstor/")) {
  setwd("/Users/mac/cloudstor/scott_projects/osteochondro_signature/") # Uses practice data (5% of cells from each sample) if running locally
  place <- "local"
} else {
  setwd("/share/ScratchGeneral/scoyou/sarah_projects/osteochondro_signature/")
  place <- "wolfpack"
}

library("xml2")
library("stringr")
library("tidyverse")
library("igraph")
library('network')
library("ggnet")
library("ggplot2")
library("intergraph")
library("RColorBrewer")
library('pals')

for(set in c("BP")) {
  # Load significantly enriched GO terms
  GOsig <- read.csv(paste0("project_results/humanchondro_bulk/Enriched_GO", set, "_ChondroConserved.csv"), header = TRUE)

  # Make iGraph in R from xgmml cytoscape file (from revigo). From https://gist.github.com/bkutlu/5a6f3b144d88169916f586cd2080d106
  # --------------------------------------------------------------------------
  
  # Import Revigo.xgmml
  x1 <- read_xml(paste0("project_results/humanchondro_bulk/Revigo_", set, ".xgmml"))
  
  # Get the nodes with xpath expression
  xpath_nodes <- str_replace_all("/graph/node", '/', '/d1:')
  nodes <- data.frame(xml_find_all(x1, xpath_nodes, xml_ns(x1)) %>%
                        xml_attrs() %>%
                        do.call(rbind, .))
  nodes$GOID <- nodes$id
  nodes$id <- paste(nodes$GOID, nodes$label, sep = "-")
  nodes <- nodes[!grepl("obsolete", nodes$id),]
  
  # Standardises GO names
  nodes <- data.frame(nodes %>% 
                        mutate_all(str_replace_all, "P-type ", ""))
  
  # Get the edges with xpath expression
  xpath_edges <- str_replace_all("/graph/edge",'/','/d1:')
  edges <- xml_find_all(x1, xpath_edges, xml_ns(x1)) %>% 
    xml_attrs() %>%
    do.call(rbind, .)
  
  # Curate into edges dataframe
  edges <- data.frame(edges[,c(2:3)])
  idx <- match(edges$source, nodes$GOID)
  edges$from <- nodes$id [idx]
  idx <- match(edges$target, nodes$GOID)
  edges$to <- nodes$id [idx]
  edges$connection <- paste(edges$from, edges$to, sep = ":")
  edges <- edges[,c(3:5, 1:2)]
  
  # Standardises GO names
  edges <- data.frame(edges %>% 
                        mutate_all(str_replace_all, "P-type ", ""))
  
  # Create iGraph object from dataframe
  ig_go <- graph_from_data_frame(d = edges, vertices = nodes, directed = F)
  
  # Plot top terms in each cluster
  # --------------------------------------------------------------------------
  GOsig$Name <- paste(GOsig$ID, GOsig$Description, sep = "-")
  idx <- match(GOsig$Name, names(igraph::components(ig_go)$membership))
  GOsig$membership <- igraph::components(ig_go)$membership [idx]
  GOsig <- GOsig[!is.na(GOsig$membership),]
  
  # Set singletons to the same group
  idx <- match(GOsig$membership, 1:length(igraph::components(ig_go)$csize))
  GOsig$csize <- igraph::components(ig_go)$csize [idx]
  GOsig$membership[GOsig$csize == 1] <- 'Singletons'
  
  # Annotate with membership information
  net <- asNetwork(ig_go)
  idx <- match(network.vertex.names(net), names(igraph::components(ig_go)$membership))
  net %v% "membership" <- igraph::components(ig_go)$membership [idx]
  
  # Set singletons to the same group
  idx <- match(get.vertex.attribute(net, "membership"), 1:length(igraph::components(ig_go)$csize))
  net %v% "csize" <- igraph::components(ig_go)$csize [idx]
  sing <- net %v% "csize" == 1
  membership <- net %v% "membership"
  membership[sing] <- 'Singletons'
  net %v% "membership" <- membership
  
  # Annotate with GO information
  idx <- match(network.vertex.names(net), GOsig$Name)
  net %v% "Count" <- GOsig$Count [idx]
  net %v% "p.adjust" <- -log10(GOsig$p.adjust) [idx]
  net %v% "Description" <- GOsig$Description [idx]
  net %v% "Name" <- GOsig$Name [idx]
  
  # Label with top terms per group and top 6 genes (by count among terms)
  rm(Top_terms)
  rm(Top_genes)
  for(i in unique(GOsig$membership)){
    # Get top terms
    tmp <- GOsig[GOsig$membership == i, ]
    num <- ifelse(nrow(tmp) > 4 | i == 'Singletons', 3, 2)
    tmp_vec <- tmp[order(tmp$p.adjust), "Name"][1:num]
    
    if(exists("Top_terms")){
      Top_terms <- c(Top_terms, tmp_vec)
    } else {
      Top_terms <- tmp_vec
    }
    
    # Get top genes
    genes_tmp <- data.frame(table(unlist(strsplit(paste(tmp$geneID, sep="", collapse="/"), "/"))))
    genes_tmp <- as.character(genes_tmp[order(genes_tmp$Freq, decreasing = TRUE), "Var1"][1:3])
    
    if(exists("Top_genes")){
      Top_genes <- rbind(Top_genes, data.frame('Cluster' = i, 'Genes' = paste(genes_tmp, sep="", collapse=" ")))
    } else {
      Top_genes <- data.frame('Cluster' = i, 'Genes' = paste(genes_tmp, sep="", collapse=" "))
    }
  }
  
  # Label clusters with unique genes
  idx <- match(net %v% "membership", Top_genes$Cluster)
  net %v% "membership" <- Top_genes$Genes [idx]
  
  idx <- match(GOsig$membership, Top_genes$Cluster)
  GOsig$Top_genes <- Top_genes$Genes [idx]
  
  # Plot the GOBP network
  getPalette <- cols25(n = length(unique(net %v% "membership")))
  names(getPalette) <- unique(net %v% "membership")
  ggnet2(net, color = "membership", size = "p.adjust", alpha = 0.5, size.cut = TRUE, color.legend = "Cluster", palette = getPalette, label = Top_terms, label.size = 2, label.color = "grey20") +
    ggsave(paste0("project_results/humanchondro_bulk/Revigo_GO", set, "_clustering_cytoscape_w_numbers.pdf"))
  
  write.csv(GOsig, paste0("project_results/humanchondro_bulk/Enriched_GO", set, "_conserved_signature_network_annotated.csv"))
  
}

