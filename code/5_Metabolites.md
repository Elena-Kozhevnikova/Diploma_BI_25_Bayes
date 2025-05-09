## This code analyses metabolites from the stool of CD patients and alignes them within gene expression BN
It uses the files generated at the previous steps
:grey_exclamation: NOTE: Metabolites made no meaningfull contributions, so they were not used further on.:grey_exclamation:  

## Prepare data
```
library(edgeR)
library(limma)
library(dplyr)
library(ggplot2)
library(tidyr)
library(bnlearn)
library(igraph)
library(scales)

# Load metadata
metadata <- read.csv("Israel_metadata.csv", na.strings = c("na", "NA", ""))
metadata <- metadata[complete.cases(metadata$pn_ID), ]  # Remove rows with missing IDs
metadata <- metadata[, !(colnames(metadata) %in% c("X", "X.1"))]

# Load metabolites
metabolites <- read.csv("Metabolites_Israel.csv")
rownames(metabolites) <- metabolites[,1]
metabolites <- metabolites[,-1]
metabolites <- t(metabolites)  # Transpose so samples are rows and metabolites are columns

# Ensure metabolites is a numeric matrix
metabolites <- apply(metabolites, 2, as.numeric)
rownames(metabolites) <- colnames(read.csv("Metabolites_Israel.csv"))[-1]

# Handle missing values by column mean imputation
metabolites <- apply(metabolites, 2, function(x) {
  x[is.na(x)] <- mean(x, na.rm = TRUE)
  return(x)
})

# Process metadata
metadata$Patient_group <- factor(metadata$Patient_group)
metadata$Gender <- factor(metadata$Gender)

# Impute missing values for Calprotectin and CRP
metadata <- metadata %>%
  group_by(Patient_group) %>%
  mutate(
    Calprotectin_ug_g = ifelse(
      is.na(Calprotectin_ug_g),
      mean(Calprotectin_ug_g, na.rm = TRUE),
      Calprotectin_ug_g
    ),
    CRP_mg_L = ifelse(
      is.na(CRP_mg_L),
      mean(CRP_mg_L, na.rm = TRUE),
      CRP_mg_L
    )
  ) %>%
  ungroup()

# Convert Age to factor and create dummy variables
metadata$Age <- factor(metadata$Age, levels = c("20_29", "30_39", "40_49", "50_59", "60_69", "70_79"))
metadata <- na.omit(metadata)

# Find common samples
common_samples <- intersect(rownames(metabolites), metadata$pn_ID)
metabolites <- metabolites[rownames(metabolites) %in% common_samples, ]
metadata_age_common <- metadata[metadata$pn_ID %in% common_samples, ]


metadata_age_common$Patient_group <- factor(metadata_age_common$Patient_group,
                                            levels = c("Israel_control", "Israel_CD"))

age_dummies <- model.matrix(~ Age - 1, data = metadata_age_common)
colnames(age_dummies) <- gsub("Age", "Age_", colnames(age_dummies))

metadata_age_common <- cbind(metadata_age_common, age_dummies)
```

## Differential expression analysis

```
# Create design matrix
design <- model.matrix(~ Patient_group + Gender + Age_20_29 +
                         Age_30_39 + Age_40_49 + Age_50_59 + 
                         Age_60_69 + Age_70_79,
                       data = metadata_age_common)

design <- design[, !colnames(design) %in% "Age_70_79"]

# Create DGEList
dge <- DGEList(counts = t(metabolites))

# Estimate dispersions
dge <- estimateDisp(dge, design)

# Create DGEList object
dge <- DGEList(counts = t(metabolites), 
               group = metadata_age_common$Patient_group,
               genes = colnames(metabolites))

# Filter low-abundance metabolites (keep features with at least 5 counts in at least 20% of samples)
keep <- filterByExpr(dge, design = design)
dge <- dge[keep, , keep.lib.sizes=FALSE]

# Normalization (TMM for metabolomics data)
dge <- calcNormFactors(dge, method = "TMM")

# Data Exploration (QC)
plotMDS(dge, col = as.numeric(dge$samples$group), 
        main = "MDS Plot of Metabolite Profiles")
legend("topright", legend = levels(dge$samples$group), 
       col = 1:nlevels(dge$samples$group), pch = 16)

# Differential Analysis
# Estimate dispersions
dge <- estimateDisp(dge, design)

# Fit the model
fit <- glmQLFit(dge, design, robust = TRUE)

colnames(design)
table(metadata_age_common$Patient_group)

# Make contrasts: Comparing all groups to control
contrasts <- makeContrasts(
  CD_vs_Control = Patient_groupIsrael_CD,  # This single term gives CD vs Control
  levels = design
)

# Test for differential metabolites
results <- list()
for (i in 1:ncol(contrasts)) {
  res <- glmQLFTest(fit, contrast = contrasts[,i])
  results[[colnames(contrasts)[i]]] <- topTags(res, n = Inf)$table
}

# Results Interpretation
CD_vs_control <- results$CD_vs_Control

# Add metabolite names if they were lost
CD_vs_control$Metabolite <- rownames(CD_vs_control)

# Filter significant results (FDR < 0.05 and absolute logFC > 1)
significant_metabolites <- CD_vs_control %>%
  filter(FDR < 0.05, abs(logFC) > 1) %>%
  arrange(FDR)

# View top results
print(significant_metabolites)

# Volcano plot
volcano_metab <- ggplot(significant_metabolites, aes(x = logFC, y = -log10(FDR))) +
  geom_point(aes(color = FDR < 0.05 & abs(logFC) > 1), alpha = 0.6) +
  scale_color_manual(values = c("grey", "red")) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  labs(title = "Volcano Plot: CD vs Control",
       x = "Log2 Fold Change",
       y = "-Log10 FDR") +
  theme_minimal(base_size = 14)

print(volcano_metab)
```

## Add significant metabolites into the BN (with genes)

```
# Add expression data
datExpr <- readRDS("datExpr.rds")
hub_gene_df <- read.csv("WGCNA_hub_genes_for_Bayesian_network_42_samples.csv")

sig_metab_names <- significant_metabolites$Metabolite

# Extract significant metabolites
metab_data <- as.data.frame(metabolites[, sig_metab_names])

# Ensure sample names match between datasets
rownames(metab_data) <- rownames(metabolites)

common_samples <- intersect(rownames(datExpr), rownames(metab_data))

# Remove immunoglubulin genes
filtered_hub_genes <- hub_gene_df %>%
  filter(!grepl("^(IG|IGH|IGK|IGL|SLC)", gene))

# Combine the datasets
combined_data <- cbind(
  datExpr[common_samples, filtered_hub_genes$gene],
  metab_data[common_samples, ]
)

# Verify dimensions
dim(combined_data)

# Learn structure with hill-climbing on combined data
bn_model_combined <- hc(combined_data, score = "bic-g")

# Bootstrapped structure learning with metabolites
boot_strength_combined <- boot.strength(combined_data, algorithm = "hc", R = 100)
avg_bn_combined <- averaged.network(boot_strength_combined)

# Create a blacklist (prevent arcs between genes in same module)
blacklist_combined <- data.frame(
  from = c(filtered_hub_genes$gene, sig_metab_names),
  to = c(filtered_hub_genes$gene, sig_metab_names)
)

# Run HC with blacklist
bn_combined <- hc(combined_data, 
                  blacklist = blacklist_combined,
                  score = "bic-g")

# Convert to igraph and calculate centrality
bn_metab_igraph <- as.igraph(bn_combined)

# Calculate out-degree centrality
degree_metab_df <- data.frame(
  gene = names(degree(bn_metab_igraph, mode = "out")),
  out_degree = degree(bn_metab_igraph, mode = "out")
)

degree_metab_df <- degree_metab_df[order(-degree_metab_df$out_degree), ]

# Identify top influencers (top 5%)
top_metab_influencers <- head(degree_metab_df, n = max(6, round(nrow(degree_metab_df)*0.05)))$gene
print(top_metab_influencers)

# Size nodes by influence (scale degree to 8-15 range)
V(bn_metab_igraph)$size <- scales::rescale(degree(bn_metab_igraph, mode = "out"), to = c(8, 15))

# Label all nodes
V(bn_metab_igraph)$label <- names(V(bn_metab_igraph))

# Color labels: red for top influencers, black for others
label_colors <- ifelse(names(V(bn_metab_igraph)) %in% top_metab_influencers, "red", "black")

# Print top influencers with their degrees
cat("\nTop influential genes (red labels):\n")
print(head(degree_metab_df, 10))

# Get edge strengths from your bootstrapping results
edge_strengths <- boot_strength_combined[boot_strength_combined$strength >= 0.7, ]
print(edge_strengths[order(-edge_strengths$strength), ])
```

## Build the graph

```
edges_metab_df <- boot_strength_combined %>%
  mutate(edge = paste(from, to, sep = "|")) %>%
  group_by(from, to) %>%
  summarize(mean_strength = mean(strength), .groups = "drop")

# Choose edges with strength >= 0.7
strong_edges_metab_df <- edges_metab_df %>% filter(mean_strength >= 0.7)

# Create te edge matrix for igraph
edge_metab_list <- as.matrix(strong_edges_metab_df[, c("from", "to")])

# Create graph
filtered_metab_g <- graph_from_edgelist(edge_metab_list, directed = TRUE)
min_edges <- 2
filtered_metab_g <- delete_vertices(
  filtered_metab_g,
  degree(filtered_metab_g) < min_edges
)
# Add weight
E(filtered_metab_g)$weight <- strong_edges_metab_df$mean_strength

# Create plot
module_palette <- c(
  "yellow" = "#FFB3BA",
  "blue" = "#BAE1FF",
  "turquoise" = "#FFFFBA",
  "brown" = "#B5EAD7",
  "green" = "grey"
  )

node_names <- V(bn_metab_igraph)$name
modules_for_nodes <- setNames(rep(NA, length(node_names)), node_names)

modules_for_nodes[names(modules_for_nodes) %in% filtered_hub_genes$gene] <- 
  filtered_hub_genes$module[match(names(modules_for_nodes)[names(modules_for_nodes) %in% 
                                                      filtered_hub_genes$gene], filtered_hub_genes$gene)]

V(bn_metab_igraph)$color <- module_palette[modules_for_nodes[V(bn_metab_igraph)$name]]

pdf("network_colored_by_module_metab2.pdf", width = 40, height = 40)
plot(bn_metab_igraph,
     vertex.size = rescale(degree(filtered_g), to = c(10, 20)),
     vertex.label.color = label_colors,  # Color by influence
     vertex.size = scales::rescale(degree(bn_metab_igraph), to = c(8, 15)),
     vertex.label = V(bn_metab_igraph)$name,
     vertex.label.cex = 3.0, 
     vertex.color = V(bn_metab_igraph)$color,
     edge.arrow.size = 2,
     edge.width = 1.5,
     layout = layout_with_kk)

legend("bottomleft",
       legend = paste("Strong edges (n =", ecount(filtered_metab_g), ")"),
       col = "red", lwd = 3, bty = "n")

dev.off()
```
## Analyze graph

```
# Filter edges by strength
strong_edges_metab_df <- edges_metab_df %>% filter(mean_strength > 0.5)

# Create filtered graph
filtered_metab_g <- graph_from_data_frame(
  d = strong_edges_metab_df[, c("from", "to", "mean_strength")],
  directed = TRUE
)

# Apply minimum edge degree filter (optional)
min_edges <- 1
filtered_metab_g <- delete_vertices(
  filtered_metab_g,
  degree(filtered_metab_g) < min_edges
)

# Apply custom module palette
module_palette <- c(
  "yellow" = "#FFB3BA",
  "blue" = "#BAE1FF", 
  "turquoise" = "#FFFFBA",
  "brown" = "#B5EAD7",
  "green" = "grey"
)

# Get module assignments for nodes in the FILTERED graph
filtered_nodes <- V(filtered_metab_g)$name
node_colors <- rep("grey", length(filtered_nodes)) # Default color

# Match nodes to their modules
matched_modules <- filtered_hub_genes$module[
  match(filtered_nodes, filtered_hub_genes$gene)
]

# Apply palette only to nodes with known modules
known_modules <- !is.na(matched_modules)
node_colors[known_modules] <- module_palette[matched_modules[known_modules]]

V(filtered_metab_g)$color <- node_colors

# Enhanced visualization
pdf("network_filtered_with_palette.pdf", width = 40, height = 40)

plot(filtered_metab_g,
     vertex.size = scales::rescale(degree(filtered_metab_g), to = c(10, 20)),
     vertex.label.cex = 3.0,
     vertex.frame.color = NA,
     vertex.label.color = "black",
     edge.arrow.size = 2,
     edge.width = scales::rescale(E(filtered_metab_g)$mean_strength, to = c(1, 4)),
     edge.color = adjustcolor("darkgrey", alpha.f = 0.7),
     layout = layout_with_lgl,
     main = "Filtered Metabolic Network (Strength ≥ 0.7)")

legend("bottomleft",
       legend = c(names(module_palette), "Other"),
       fill = c(module_palette, "grey"),
       title = "Modules",
       bty = "n",
       cex = 2)

legend("bottomright",
       legend = c(paste("Nodes:", vcount(filtered_metab_g)),
                  paste("Edges:", ecount(filtered_metab_g)),
                  "Edge width = Connection strength",
                  "Node size = Degree"),
       bty = "n",
       cex = 2)

dev.off()
```
