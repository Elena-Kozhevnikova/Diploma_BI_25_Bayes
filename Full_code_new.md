## 1. Prepare data
```
setwd("C:/Users/Елена/Downloads/Diploma_2025")

# Load packages
library(edgeR)
library(tidyverse)
library(WGCNA)
library(clusterProfiler)
library(ggplot2)
library(org.Hs.eg.db)
library(dplyr)


# Read data with duplicate handling
count_data <- read_csv("fc_data.csv", show_col_types = FALSE) %>%
  # Aggregate duplicate genes by summing TPM values
  group_by(gene_name, sample_id) %>%
  summarise(tpm = sum(tpm), .groups = 'drop') %>%
  # Convert to wide format
  pivot_wider(
    names_from = sample_id,
    values_from = tpm,
    values_fill = 0
  ) %>%
  # Now safe to set row names
  column_to_rownames("gene_name") %>%
  as.matrix()

# For count data (samples = columns)
count_samples <- colnames(count_data)

# For metadata
metadata <- read.csv("Israel_metadata.csv", row.names = "pn_ID")

# Keep only metadata for samples present in count data
matched_metadata <- metadata[count_samples, , drop = FALSE]

# Find missing samples
missing_in_metadata <- setdiff(count_samples, rownames(metadata))
missing_in_counts <- setdiff(rownames(metadata), count_samples)

# Remove empty columns
matched_metadata <- matched_metadata[, colSums(is.na(matched_metadata)) < nrow(matched_metadata)]

# Convert NAs
matched_metadata[is.na(matched_metadata)] <- "NA"

# View structure
str(matched_metadata)

# Find intersecting samples
common_samples <- intersect(colnames(count_data), rownames(metadata))

# Subset both objects to matching samples
count_data <- count_data[, common_samples]
metadata <- metadata[common_samples, ]

# Now create DGEList
y <- DGEList(counts = count_data, group = metadata$Patient_group)
```

## Perform DGE analysis
```
# Filter genes with low expression
keep <- filterByExpr(y, group = y$samples$group)
y <- y[keep, , keep.lib.sizes = FALSE]

y <- calcNormFactors(y)

design <- model.matrix(~0 + Patient_group + Age + Gender, data = metadata)
colnames(design) <- make.names(colnames(design))

y <- estimateDisp(y, design, robust = TRUE)

fit <- glmQLFit(y, design)

# Define contrast (CD vs Control)
contrast <- makeContrasts(
  "Patient_groupIsrael_CD - Patient_groupIsrael_control",
  levels = design
)

qlf <- glmQLFTest(fit, contrast = contrast)

results <- topTags(qlf, n = Inf, sort.by = "PValue")$table

# Filter significant genes (FDR < 0.05, |logFC| > 1)
sig_genes <- results %>%
  filter(FDR < 0.05, abs(logFC) > 1) %>%
  arrange(FDR)

# Save results
write_csv(results, "full_dge_results.csv")
write_csv(sig_genes, "significant_genes.csv")

pdf("dge_qc_plots.pdf")
plotMDS(y, col = as.numeric(factor(metadata$Patient_group)))
plotBCV(y)
plotMD(qlf)
dev.off()
```
## WCGNA

```
options(stringsAsFactors = FALSE)

# Start with DEGs
deg_counts <- count_data[rownames(sig_genes), ]

# Proceed with WGCNA steps from earlier, using:
datExpr <- t(deg_counts)  # Genes as columns

# Filter lowly expressed genes
keep_genes <- rowSums(datExpr > 1) >= 0.5*ncol(datExpr)  # Expressed in >50% samples
datExpr <- datExpr[keep_genes, ]

# Check for outliers
# sampleTree <- hclust(dist(datExpr), method = "average")
# plot(sampleTree, main = "Sample clustering")

########### BIG DEAL OUTLIERS #############
# Remove with PCA - removes 3 CD samples
# Proper PCA on samples (rows) × genes (columns)
# pca <- prcomp(datExpr, scale. = TRUE)  # Remove t() transpose

# Calculate outlier scores (Mahalanobis distance)
# library(mvoutlier)
# outlier_scores <- pca$x[,1]^2 + pca$x[,2]^2  # Sum of squared PC1 and PC2 scores

# Identify outlier SAMPLES (not genes)
# bad_samples <- rownames(datExpr)[outlier_scores > quantile(outlier_scores, 0.95)]  # Top 5% outliers
# print(bad_samples)

# Remove with dynamic clusters - leaves only 21 samples, strongly unbalanced

###### FINALLY I KEPT ALL SAMPLES! ##############################

# Choose power based on scale-free topology fit
powers <- c(seq(4, 10, by = 1), seq(12, 16, by = 2))
sft <- pickSoftThreshold(datExpr, 
                         powerVector = powers, 
                         networkType = "signed",
                         verbose = 5)

# Plot results
par(mfrow = c(1, 2))
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab = "Soft Threshold Power", 
     ylab = "Scale Free Topology Model Fit (signed R^2)",
     main = "Scale Independence")
abline(h = 0.80, col = "red")  # Adjusted threshold for DEGs

plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab = "Soft Threshold Power", 
     ylab = "Mean Connectivity",
     main = "Mean Connectivity")
dev.off()

# Select power (choose lowest power where R^2 > 0.80)
softPower <- 5  # Typically 6-8 for DEGs

# One-step network construction (auto-detects blocks)
net <- blockwiseModules(
  datExpr,
  power = softPower,
  networkType = "signed",
  TOMType = "signed",
  minModuleSize = 20,
  mergeCutHeight = 0.20,
  deepSplit = 2,
  pamRespectsDendro = FALSE,
  verbose = 3
)

# Module colors and labels
moduleColors <- net$colors
table(moduleColors)

# Plot dendrogram
plotDendroAndColors(
  net$dendrograms[[1]],
  moduleColors[net$blockGenes[[1]]],
  "Module colors",
  dendroLabels = FALSE,
  hang = 0.03,
  addGuide = TRUE,
  guideHang = 0.05
)

pdf("dendrogram_42_samples.pdf", width = 10, height = 6)
par(mar = c(6, 8.5, 3, 3))
plotDendroAndColors(
  net$dendrograms[[1]],
  moduleColors[net$blockGenes[[1]]],
  "Module colors",
  dendroLabels = FALSE,
  hang = 0.03,
  addGuide = TRUE,
  guideHang = 0.05
)
dev.off()
```
## Module-trait associations

```
# Prepare traits with proper type conversion and NA handling
traitData <- metadata[rownames(datExpr), ] %>%
  mutate(
    CD_status = as.numeric(Patient_group == "Israel_CD"),
    Male = as.numeric(Gender == "male"),
    # Convert to numeric, handling "na" strings
    CRP_mg_L = as.numeric(ifelse(CRP_mg_L == "na", NA, CRP_mg_L)),
    Calprotectin_ug_g = as.numeric(ifelse(Calprotectin_ug_g == "na", NA, Calprotectin_ug_g)),
    # Convert age to factor
    Age = factor(Age, levels = c("20_29", "30_39", "40_49", "50_59", "60_69", "70_79"))
  )

# Create dummy variables for age
age_dummies <- model.matrix(~ Age - 1, data = traitData)
colnames(age_dummies) <- gsub("Age", "Age_", colnames(age_dummies))

# Prepare numeric traits and impute NAs with column means
numeric_traits <- traitData %>%
  select(CD_status, Male, CRP_mg_L, Calprotectin_ug_g) %>%
  mutate(across(where(is.numeric), ~ ifelse(is.na(.x), mean(.x, na.rm = TRUE), .x))) %>%
  as.matrix()

# Combine all traits
traitMatrix <- cbind(numeric_traits, age_dummies)

# Final numeric conversion (safety check)
traitMatrix <- apply(traitMatrix, 2, as.numeric)
rownames(traitMatrix) <- rownames(numeric_traits)

# Align with module eigengenes
MEs <- net$MEs
common_samples <- intersect(rownames(MEs), rownames(traitMatrix))
traitMatrix <- traitMatrix[common_samples, ]
MEs <- MEs[common_samples, ]

# Calculate correlations
moduleTraitCor <- cor(MEs, traitMatrix, use = "complete.obs")
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nSamples = nrow(datExpr))

# Plot with improved formatting
pdf("module_trait_relationships_final_42_samples.pdf", width = 12, height = 8)
age_cols <- grep("^Age_", colnames(traitMatrix))

labeledHeatmap(
  Matrix = moduleTraitCor,
  xLabels = colnames(traitMatrix),
  yLabels = names(MEs),
  ySymbols = names(MEs),
  colorLabels = FALSE,
  colors = blueWhiteRed(50),
  textMatrix = ifelse(moduleTraitPvalue < 0.05, 
                      paste0(round(moduleTraitCor, 2), "*"),
                      round(moduleTraitCor, 2)),
  setStdMargins = FALSE,
  cex.text = 0.7,
  zlim = c(-1, 1),
  main = "Module-Trait Relationships\n(* = p < 0.05)",
  xLabelsAngle = 45,
  xColorLabels = ifelse(colnames(traitMatrix) %in% colnames(traitMatrix)[age_cols],
                        "blue", "black")
)
dev.off()

# Save module assignments and colors
module_labels <- data.frame(
  gene_id = colnames(datExpr),
  module_color = moduleColors
)
write_csv(module_labels, "gene_module_assignments_42_samples.csv")

# Save module eigengenes (MEs)
write_csv(as.data.frame(MEs) %>% rownames_to_column("sample_id"), "module_eigengenes_42_samples.csv")

# Calculate intramodular connectivity
kME <- signedKME(datExpr, MEs)  # Gene-module eigengene correlation
colnames(kME) <- gsub("kME", "", colnames(kME))

# For each module, get top 20 hub genes
hub_genes <- lapply(unique(moduleColors), function(mod) {
  kME_sub <- kME[moduleColors == mod, mod]
  names(sort(kME_sub, decreasing = TRUE))[1:20]
})
names(hub_genes) <- unique(moduleColors)
write_rds(hub_genes, "hub_genes_per_module_42_samples.rds")
```

# GO and KEGG analysis of module genes

```
analyze_module <- function(module_color, gene_universe) {
  # Get module genes
  module_genes <- colnames(datExpr)[moduleColors == module_color]
  
  # Convert symbols to ENTREZID (with error handling)
  gene_map <- tryCatch({
    bitr(module_genes, 
         fromType = "SYMBOL", 
         toType = "ENTREZID", 
         OrgDb = org.Hs.eg.db)
  }, error = function(e) {
    message("ID conversion failed for ", module_color, ": ", e$message)
    return(NULL)
  })
  
  if (is.null(gene_map) || nrow(gene_map) < 5) {
    message("Skipping ", module_color, " (too few genes or mapping failed)")
    return(NULL)
  }
  
  # Prepare universe genes (if provided)
  if (!missing(gene_universe)) {
    universe_map <- bitr(gene_universe, 
                         fromType = "SYMBOL", 
                         toType = "ENTREZID", 
                         OrgDb = org.Hs.eg.db)
    universe_entrez <- na.omit(universe_map$ENTREZID)
  } else {
    universe_entrez <- NULL
  }
  
  # Run GO and KEGG enrichment
  results <- list()
  
  # GO enrichment
  results$GO <- tryCatch({
    enrichGO(
      gene = gene_map$ENTREZID,
      universe = universe_entrez,
      OrgDb = org.Hs.eg.db,
      ont = "BP",
      pAdjustMethod = "BH",
      minGSSize = 5
    )
  }, error = function(e) {
    message("GO enrichment failed for ", module_color, ": ", e$message)
    NULL
  })
  
  # KEGG enrichment
  results$KEGG <- tryCatch({
    enrichKEGG(
      gene = gene_map$ENTREZID,
      organism = "hsa",
      pAdjustMethod = "BH"
    )
  }, error = function(e) {
    message("KEGG enrichment failed for ", module_color, ": ", e$message)
    NULL
  })
  
  return(results)
}

###################################
# Define gene universe (all expressed genes)
gene_universe <- colnames(datExpr)

# Get all modules (excluding grey)
modules_to_analyze <- setdiff(unique(moduleColors), "grey")

# Run enrichment for each module
enrichment_results <- lapply(modules_to_analyze, function(mod) {
  cat("\n=== Analyzing", mod, "module ===\n")
  analyze_module(mod, gene_universe)
})
names(enrichment_results) <- modules_to_analyze

# Save results
for (mod in modules_to_analyze) {
  if (!is.null(enrichment_results[[mod]])) {
    # Save GO results
    if (!is.null(enrichment_results[[mod]]$GO)) {
      write.csv(as.data.frame(enrichment_results[[mod]]$GO),
                file = paste0(mod, "_GO_enrichment_42_samples.csv"))
    }
    # Save KEGG results
    if (!is.null(enrichment_results[[mod]]$KEGG)) {
      write.csv(as.data.frame(enrichment_results[[mod]]$KEGG),
                file = paste0(mod, "_KEGG_enrichment_42_samples.csv"))
    }
  }
}

# Plot results

# Custom safe plotting function
safe_dotplot <- function(enrich_result, title = "") {
  if (is.null(enrich_result) || nrow(enrich_result) == 0) {
    message("No enrichment results to plot")
    return(ggplot() + 
             ggtitle(paste(title, "(No significant enrichment)")) + 
             theme_void())
  }
  
  tryCatch({
    dotplot(enrich_result, showCategory = 15) + 
      ggtitle(title) +
      theme(plot.title = element_text(hjust = 0.5))
  }, error = function(e) {
    message("Plotting failed: ", e$message)
    ggplot() + 
      ggtitle(paste(title, "(Plotting error)")) + 
      theme_void()
  })
}

# Generate all plots safely
pdf("module_enrichment_plots_robust_42_samples.pdf", width = 12, height = 8)

for (mod in names(enrichment_results)) {
  cat("\nGenerating plots for", mod, "module...")
  
  # GO Plot
  if (!is.null(enrichment_results[[mod]]$GO)) {
    p <- safe_dotplot(enrichment_results[[mod]]$GO, 
                      paste(mod, "Module: GO Enrichment"))
    print(p)
  } else {
    message("No GO results for ", mod)
  }
  
  # KEGG Plot
  if (!is.null(enrichment_results[[mod]]$KEGG)) {
    p <- safe_dotplot(enrichment_results[[mod]]$KEGG, 
                      paste(mod, "Module: KEGG Pathways"))
    print(p)
  } else {
    message("No KEGG results for ", mod)
  }
}

dev.off()
```
# Double check the modules

```
# Check if modules were properly assigned
table(moduleColors)

# Calculate signed adjacency matrix (using biweight midcorrelation for robustness)
library(WGCNA)
adjacency <- adjacency(datExpr, 
                       power = softPower,
                       type = "signed",
                       corFnc = "bicor")  # More robust than Pearson

# Calculate TOM and intramodular connectivity
TOM <- TOMsimilarity(adjacency)
kIN <- intramodularConnectivity(adjacency, moduleColors, scaleByMax = TRUE)

get_hub_genes <- function(module) {
  mod_genes <- names(moduleColors)[moduleColors == module]
  
  # Calculate three hubness metrics:
  # 1. Intra-modular connectivity
  mod_connectivity <- kIN[mod_genes, "kWithin"]
  
  # 2. Module membership (kME)
  mod_kME <- cor(datExpr[, mod_genes], MEs[, paste0("ME", module)])
  
  # 3. Gene significance (correlation with traits)
  if(exists("traitMatrix")) {
    gs <- as.data.frame(cor(datExpr[, mod_genes], traitMatrix, use = "p"))
  } else {
    gs <- NULL
  }
  
  # Combine metrics (weighted sum)
  hub_score <- scale(mod_connectivity) + 
    scale(abs(mod_kME)) + 
    if(!is.null(gs)) scale(rowMeans(abs(gs))) else 0
  
  # Get top 20 hub genes
  hub_genes <- mod_genes[order(hub_score, decreasing = TRUE)[1:20]]
  return(hub_genes)
}

# Apply to all non-grey modules
hub_genes_list <- lapply(c("blue", "brown", "green", "turquoise", "yellow"), get_hub_genes)
names(hub_genes_list) <- c("blue", "brown", "green", "turquoise", "yellow")

# Plot connectivity vs kME for each module
plot_hub_validation <- function(module) {
  mod_genes <- names(moduleColors)[moduleColors == module]
  
  # Calculate required metrics
  connectivity <- kIN[mod_genes, "kWithin"]
  kME_values <- cor(datExpr[, mod_genes], MEs[, paste0("ME", module)])
  
  # Create data frame for plotting
  plot_data <- data.frame(
    Gene = mod_genes,
    Connectivity = connectivity,
    ModuleMembership = as.vector(kME_values),  # Remove matrix attributes
    IsHub = mod_genes %in% hub_genes_list[[module]]
  )
  
  # Generate the plot
    ggplot(plot_data, aes(x = Connectivity, y = abs(ModuleMembership), color = IsHub)) +
    geom_point(aes(color = IsHub), alpha = 0.6) +
    scale_color_manual(values = c("grey", "red")) +
    ggtitle(paste(module, "Module Validation")) +
    labs(x = "Intramodular Connectivity", y = "Module Membership (|kME|)") +
    theme_minimal() +
    theme(legend.position = "none")
}

# Generate all plots
pdf("hub_gene_validation_plots_42_samples.pdf", width = 10, height = 8)
for (mod in names(hub_genes_list)) {
  if (length(hub_genes_list[[mod]]) > 0) {  # Only plot modules with hubs
    print(plot_hub_validation(mod))
  } else {
    message("Skipping ", mod, " - no hub genes")
  }
}
dev.off()

# Create comprehensive hub gene data frame
hub_gene_df <- do.call(rbind, lapply(names(hub_genes_list), function(mod) {
  data.frame(
    gene = hub_genes_list[[mod]],
    module = mod,
    kWithin = kIN[hub_genes_list[[mod]], "kWithin"],
    kME = cor(datExpr[, hub_genes_list[[mod]]], MEs[, paste0("ME", mod)]),
    GS = if(exists("traitMatrix")) rowMeans(abs(cor(datExpr[, hub_genes_list[[mod]]], traitMatrix))) else NA
  )
}))
```
## Remove green module
# Here I removed the green module as it was 99% immunoglobulin genes, which are surely no key regulators
```
# Filter out immunoglobulin genes and solute carriers
filtered_hub_genes <- hub_gene_df %>%
  filter(!grepl("^(IG|IGH|IGK|IGL|SLC)", gene))
write.csv(filtered_hub_genes, "WGCNA_hub_genes_for_Bayesian_network_42_samples.csv", row.names = FALSE)
saveRDS(datExpr, file = "datExpr.rds")
```
