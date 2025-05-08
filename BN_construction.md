## This code is aimed at constructing a Bayesian Network (BN) to identify key regulator genes in Crohn's Disease
It uses data files generated previously with the code shown in DGE_WCGNA.md in this repo.

```
setwd("C:/Users/Елена/Downloads/Diploma_2025")
library(bnlearn)
library(Rgraphviz)
library(igraph)
library(scales)
library(dplyr)


# Read the hub gene list
hub_genes <- read.csv("WGCNA_hub_genes_for_Bayesian_network_42_samples.csv")

datExpr <- readRDS("datExpr.rds")

# Extract expression data for hub genes only
hub_expr <- datExpr[, hub_genes$gene]

# Subset only the hub genes (columns) and keep samples as rows
bn_data <- as.data.frame(datExpr[, colnames(datExpr) %in% hub_genes$gene])

# Learn structure with hill-climbing
bn_model <- hc(bn_data, score = "bic-g")  # Gaussian for continuous data

# Bootstrapped structure learning
boot_strength <- boot.strength(bn_data, algorithm = "hc", R = 100)
avg_bn <- averaged.network(boot_strength)

# Save bootstrap results
saveRDS(boot_strength, file = "boot_strength_results_primary.rds")

# Save averaged network
saveRDS(avg_bn, file = "averaged_network_primary.rds")

# Plot the network
plot(bn_model, main = "Bayesian Network Structure")

boot_strength <- readRDS("boot_strength_results_primary.rds")
avg_bn <- readRDS("averaged_network_primary.rds")

# Plot with arc strengths
par(mar = c(1, 1, 1, 1))  # Reduce margins
dev.new(width = 10, height = 10)  # Open a larger graphics device (adjust as needed)
strength.plot(avg_bn, boot_strength, shape = "ellipse")

# Create a blacklist (prevent arcs between genes in same module)
blacklist <- data.frame(
  from = hub_genes$gene,
  to = hub_genes$gene
)

# Assess arc confidence (from bootstrapping)
print(boot_strength[boot_strength$strength > 0.7 & boot_strength$direction > 0.7, ])

# Plot only high-confidence arcs (e.g., strength > 0.7)
strength.plot(avg_bn, boot_strength, threshold = 0.7, shape = "ellipse")

# Identify top influencers (nodes with highest out-degree)
library(igraph)
bn_igraph <- as.igraph(avg_bn)  # Convert to igraph for centrality
degree_df <- data.frame(
  gene = names(degree(bn_igraph, mode = "out")),
  out_degree = degree(bn_igraph, mode = "out")
)
degree_df <- degree_df[order(-degree_df$out_degree), ]
print(head(degree_df))

library(igraph)
bn_igraph <- as.igraph(avg_bn)

# Set vertex sizes
V(bn_igraph)$size <- 10  # Adjust size as needed

# Plot
pdf("network_igraph_42_samples.pdf", width = 40, height = 40)
plot(bn_igraph, vertex.label.cex = 1.2)
dev.off()


# Plot only high-confidence arcs (e.g., strength > 0.7)
strength.plot(avg_bn, boot_strength, threshold = 0.7, shape = "ellipse") 

pdf("network_simple_42_samples.pdf", width = 30, height = 30)
strength.plot(avg_bn, boot_strength, threshold = 0.7, 
              shape = "ellipse", fontsize = 24)
dev.off()

#################################
# Convert to igraph and calculate centrality
library(igraph)

# Convert to igraph and calculate centrality
bn_igraph <- as.igraph(avg_bn)

# Calculate out-degree centrality
degree_df <- data.frame(
  gene = names(degree(bn_igraph, mode = "out")),
  out_degree = degree(bn_igraph, mode = "out")
)
degree_df <- degree_df[order(-degree_df$out_degree), ]

# Identify top influencers (top 5%)
top_influencers <- head(degree_df, n = max(6, round(nrow(degree_df)*0.05)))$gene

# Size nodes by influence (scale degree to 8-15 range)
V(bn_igraph)$size <- scales::rescale(degree(bn_igraph, mode = "out"), to = c(8, 15))

# Label all nodes
V(bn_igraph)$label <- names(V(bn_igraph))

# Color labels: red for top influencers, black for others
label_colors <- ifelse(names(V(bn_igraph)) %in% top_influencers, "red", "black")

# Create the plot
pdf("fully_labeled_network_42_samples.pdf", width = 40, height = 40)
plot(bn_igraph,
     edge.arrow.size = 2,  # Controls arrowhead size (default: 1)
     edge.arrow.width = 1.2, # Controls arrowhead width 
     edge.arrow.mode = 2,    # 2 = directed arrows, 1 = backward, 0 = none
     vertex.label = V(bn_igraph)$label,  # Show all labels
     vertex.label.color = label_colors,  # Color by influence
     vertex.label.cex = 3,  # Slightly smaller for all labels
     vertex.label.font = 15,   # Bold font
     vertex.shape = "circle", # All nodes as circles
     edge.arrow.size = 0.3,
     layout = layout_with_fr(bn_igraph, weights = E(bn_igraph)$weight))  # Weighted layout

# Add enhanced legend
legend("topleft",
       legend = c("Blue module", "Turquoise module", "Other", "Top influencer (label)"),
       pt.bg = c("lightblue", "lightgreen", "gray", NA),
       col = c(NA, NA, NA, "red"),  # Only show color for influencer text
       pch = 21,  # All circles
       pt.cex = 12,
       text.col = c("black", "black", "black", "red"),  # Match label colors
       bty = "n")

# Add title with network metrics
title(paste("Bayesian Network (", length(V(bn_igraph)), "genes,", 
            length(E(bn_igraph)), "connections)"), cex.main = 5)

dev.off()

# Print top influencers with their degrees
cat("\nTop influential genes (red labels):\n")
print(head(degree_df, 10))

# Get edge strengths from your bootstrapping results
edge_strengths <- boot_strength[boot_strength$strength >= 0.7, ]  # Filter for significant edges
print(edge_strengths[order(-edge_strengths$strength), ]) 
```

# Apply bootstrapping results to the initial network
```
# Choose edges evealed after bootstrapping
edges_df <- boot_strength %>%
  mutate(edge = paste(from, to, sep = "|")) %>%
  group_by(from, to) %>%
  summarize(mean_strength = mean(strength), .groups = "drop")

# Choose edges with strength >= 0.7
strong_edges_df <- edges_df %>% filter(mean_strength >= 0.7)

# Create edge matrix with igraph
edge_list <- as.matrix(strong_edges_df[, c("from", "to")])

# Create graph
filtered_g <- graph_from_edgelist(edge_list, directed = TRUE)

# Add weight
E(filtered_g)$weight <- strong_edges_df$mean_strength

# Visualize ONLY the strong connections

# palette for modules
module_palette <- c(
  "yellow" = "#FFB3BA",
  "blue" = "#BAE1FF",
  "turquoise" = "#FFFFBA",
  "brown" = "#B5EAD7",
  "green" = "grey"
  )

# node names are genes
node_names <- V(bn_igraph)$name
modules_for_nodes <- setNames(rep(NA, length(node_names)), node_names)

modules_for_nodes[names(modules_for_nodes) %in% hub_genes$gene] <- 
  hub_genes$module[match(names(modules_for_nodes)[names(modules_for_nodes) %in% 
                                                      hub_genes$gene], hub_genes$gene)]
# pply palette based on node color
V(bn_igraph)$color <- module_palette[modules_for_nodes[V(bn_igraph)$name]]

# Save the graph
pdf("network_colored_by_module_final.pdf", width = 40, height = 40)
plot(bn_igraph,
     vertex.size = rescale(degree(filtered_g), to = c(10, 20)),
     vertex.label.color = label_colors,  # Color by influence
     vertex.size = scales::rescale(degree(bn_igraph), to = c(8, 15)),
     vertex.label = V(bn_igraph)$name,
     vertex.label.cex = 3.0, 
     vertex.color = V(bn_igraph)$color,
     edge.arrow.size = 2,
     edge.width = 1.5,  # ширина ребер
     layout = layout_with_kk)

legend("bottomleft",
       legend = paste("Strong edges (n =", ecount(filtered_g), ")"),
       col = "red", lwd = 3, bty = "n")

dev.off()

# Filter for statistically significant edges (strength >= 0.7)
significant_edges <- boot_strength[boot_strength$strength >= 0.7, ]

# Sort by strength (most to least significant)
significant_edges <- significant_edges[order(-significant_edges$strength), ]

# Print results
print(significant_edges)
```

## GO analysis of the resulting network
```
# Select genes with string connections
genes <- unique(c(significant_edges$from, significant_edges$to))

# Convert gene symbols to Entrez IDs (required by clusterProfiler)
entrez_ids <- bitr(genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
print(genes)

# Biological Process (BP) - MF" (Molecular Function) or "CC" (Cellular Component)
go_bp <- enrichGO(
  gene          = entrez_ids$ENTREZID,
  OrgDb         = org.Hs.eg.db,
  ont           = "CC",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.2,
  readable      = TRUE
)

# View top terms
head(go_bp)

# Plot
dotplot_go <- dotplot(go_bp, showCategory = 20, title = "GO Cellular Component")
print(dotplot_go)

pdf("GO_BN_over_0.7.pdf", width = 7, height = 7)
dotplot_go
dev.off()

# Save GO results
write.csv(go_bp@result, "GO_CC_enrichment.csv")
```

## Find the longest connection within the network
```
# Convert Bayesian network to igraph
bn_igraph <- as.igraph(avg_bn)

# Calculate diameter (longest shortest path)
diameter_info <- diameter(bn_igraph, directed = TRUE, unconnected = TRUE, weights = NA)
cat("Network diameter (longest shortest path):", diameter_info, "\n")

# Get the actual longest path
longest_path <- get_diameter(bn_igraph, directed = TRUE)
cat("Longest path:", paste(names(longest_path), collapse = " → "), "\n")

##### Longest connection for edges over 0.7:

# Calculate the longest path using THE SAME GRAPH OBJECT that I plotted
longest_path <- get_diameter(filtered_g, directed = TRUE)

# Now highlight using vertex indices from it
V(filtered_g)$color <- "lightgray"
V(filtered_g)[longest_path]$color <- "red"  # Uses indices from filtered_g

# Highlight edges along the path
E(filtered_g)$color <- "gray"
for(i in 1:(length(longest_path)-1)) {
  E(filtered_g)[longest_path[i] %->% longest_path[i+1]]$color <- "red"
  E(filtered_g)[longest_path[i] %->% longest_path[i+1]]$width <- 3
}

# Visualize
pdf("longest_path_highlighted.pdf", width = 15, height = 15)
plot(filtered_g, 
     layout = layout_with_fr,
     main = paste("Longest Path (Length =", length(longest_path)-1, ")"),
     edge.arrow.size = 0.5)
legend("topleft", 
       legend = c("Longest Path", "Other Nodes"),
       col = c("red", "lightgray"), 
       pch = 19, cex = 1.2)
dev.off()
```
