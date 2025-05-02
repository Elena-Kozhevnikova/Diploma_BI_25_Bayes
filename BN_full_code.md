# Data aquisition and pre-processing

The Isratel patient sequence data was downloaded from here: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE199906
Metabolomic and metadata were downloaded from the aticle supplementary file: https://www.nature.com/articles/s41467-024-48106-6

File processing was performed in command line as follows:

# Combine the data from all archive files
```
tar -xvf GSE199906_RAW.tar

zcat GSM*.txt.gz | awk '
   BEGIN { FS="\t"; OFS="," }      # Set input/output field separators
  {
    split($1, parts, "|");        # Split target_id by "|"
    print parts[1], parts[6], $2  # Print transcript_id, gene_name, TPM
  }' > extracted_data.csv

  cat extracted_data.csv --head
  cat extracted_data.csv

  echo "sample_id,gene_name,tpm,transcript_id" > final_output.csv
  find . -name "*.txt.gz" -print0 | parallel -0 -j $(nproc) "
    zcat {} | awk -v file=\"{}\" '
    BEGIN {FS=\"\t\"; OFS=\",\"}
      {
       # Use full filename (without path) as sample_id
       sample_id = file
        sub(/^.*\//, \"\", sample_id)  # Remove directory path
        sub(/\.txt\.gz$/, \"\", sample_id)  # Remove extension
        split(\$1, parts, \"|\")
        if (length(parts) >= 6 && parts[6] != \"\") {
          print sample_id, parts[6], \$2, parts[1]
        }
      }
    '
  " >> final_output.csv


  echo "Total .txt.gz files found:"
  find . -name "*.txt.gz" | wc -l

  echo "sample_id,gene_name,tpm,transcript_id" > final_output1.csv
  find . -name "*.txt.gz" -print0 | while IFS= read -r -d $'\0' file; do   echo "Processing: $file" >&2;   zcat "$file" | awk -v file="$file" '
      BEGIN {FS="\t"; OFS=","; processed=0}
      {
        # Extract base filename without path/extension
        sample_id = file
        sub(/^.*\//, "", sample_id)
        sub(/\.txt\.gz$/, "", sample_id)

        # Count processed records
        split($1, parts, "|")
        if (length(parts) >= 6 && parts[6] != "") {
          print sample_id, parts[6], $2, parts[1]
          processed++
        }
      }
      END {
        if (processed == 0) {
          print "DEBUG: No valid records in " file > "/dev/stderr"
        }
      }
    ' >> final_output1.csv; qz
```
# Change sanple names
```
  echo "sample_id,gene_name,tpm,transcript_id" > final_output_combined.csv
  find . -name "*.txt.gz" -print0 | while IFS= read -r -d $'\0' file; do   zcat "$file" | awk -v file="$file" '
      BEGIN {FS="\t"; OFS=","}
      NR > 1 {
        sample_id = file
        sub(/^.*\//, "", sample_id)
        sub(/\.txt\.gz$/, "", sample_id)

        # Detect format
        if ($1 ~ /\|/) {
          # Pipe-delimited format (GSM691*)
          split($1, parts, "|")
          print sample_id, parts[6], $2, parts[1]
        } else if ($1 ~ /_/) {
          # Underscore-delimited format (GSM599*)
          split($1, parts, "_")
          print sample_id, parts[1], $2, parts[2]
        }
      }
    ' >> final_output_combined.csv; done

  awk -F, 'NR>1 {print $1}' final_output.csv | sort | uniq | wc -l
  awk -F, 'NR>1 {print $1}' final_output_combined.csv | sort | uniq | wc -l
  head -n 5 final_output_combined.csv && echo "..." && tail -n 5 final_output_combined.csv

  awk 'BEGIN {FS=OFS="\t"}
       NR==FNR {a[$1]=$2 OFS $3 OFS $4; next}
       $1 in a {print $0, a[$1]}' metadata.tsv expression.tsv > merged_data.tsv


  cat merged_data.tsv | head
output
  A025    DDX11L1 0.0     ENST00000456328.2       Israel_CD       30_39   male
  A025    DDX11L1 0.0     ENST00000450305.2       Israel_CD       30_39   male
  A025    WASH7P  1.85587 ENST00000488147.1       Israel_CD       30_39   male
  A025    MIR6859-1       0.0     ENST00000619216.1       Israel_CD       30_39   male
  A025    RP11-34P13.3    0.0     ENST00000473358.1       Israel_CD       30_39   male
  A025    RP11-34P13.3    0.0     ENST00000469289.1       Israel_CD       30_39   male
  A025    MIR1302-2       0.0     ENST00000607096.1       Israel_CD       30_39   male
  A025    FAM138A 0.0     ENST00000417324.1       Israel_CD       30_39   male
  A025    FAM138A 0.0209104       ENST00000461467.1       Israel_CD       30_39   male

  echo -e "sample_id\tgene_name\ttpm\ttranscript_id\tPatient_group\tAge\tGender" | cat - merged_data.tsv > headed_data.tsv
```
# Perform DGE in command line:
  awk 'BEGIN {
      FS=OFS="\t";
      print "gene_name\tmean_CD\tmean_Control\tlog2FC\tCD_samples\tControl_samples";
  }
  NR==1 {next} # Skip header
  {
      genes[$2][$5] += $3; # Sum TPM by gene and group
      counts[$2][$5]++;    # Count samples per group
  }
  END {
      for (gene in genes) {
          cd_mean = (counts[gene]["Israel_CD"] > 0) ? genes[gene]["Israel_CD"]/counts[gene]["Israel_CD"] : 0;
          ctrl_mean = (counts[gene]["Israel_control"] > 0) ? genes[gene]["Israel_control"]/counts[gene]["Israel_control"] : 0;

          # Avoid division by zero in fold change calculation
          log2fc = (ctrl_mean > 0) ? log(cd_mean/ctrl_mean)/log(2) : "NA";

          print gene, cd_mean, ctrl_mean, log2fc,
                counts[gene]["Israel_CD"]+0, counts[gene]["Israel_control"]+0;
      }
  }' headed_data.tsv > dge_results.tsv

### The next steps are performed in R

# Prepare data:
```
library(edgeR)
library(tidyverse)

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

count_samples <- colnames(count_data)
n_samples <- length(count_samples)

metadata <- read.csv("Israel_metadata.csv", row.names = "pn_ID")

matched_metadata <- metadata[count_samples, , drop = FALSE]

# Find missing samples
missing_in_metadata <- setdiff(count_samples, rownames(metadata))
missing_in_counts <- setdiff(rownames(metadata), count_samples)

if(length(missing_in_metadata) > 0) {
  warning("Missing metadata for: ", paste(missing_in_metadata, collapse = ", "))
}

if(length(missing_in_counts) > 0) {
  warning("Missing counts for: ", paste(missing_in_counts, collapse = ", "))
}
# Remove empty columns
matched_metadata <- matched_metadata[, colSums(is.na(matched_metadata)) < nrow(matched_metadata)]

# Convert NAs
matched_metadata[is.na(matched_metadata)] <- "NA"

# Find intersecting samples
common_samples <- intersect(colnames(count_data), rownames(metadata))
cat("Common samples:", length(common_samples), "\n")

# 3. Subset both objects to matching samples
count_data <- count_data[, common_samples]
metadata <- metadata[common_samples, ]

# 5. Now create DGEList
y <- DGEList(counts = count_data, group = metadata$Patient_group)
```

# Perfrom DGE using EdgeR: 
```
keep <- filterByExpr(y, group = y$samples$group)
y <- y[keep, , keep.lib.sizes = FALSE]

y <- calcNormFactors(y)

# Accound for factors Age and Gender
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
### WCGNA
# 
```
#################################################
library(WGCNA)
options(stringsAsFactors = FALSE)

# Start with 1000+ DEGs (from edgeR)
deg_counts <- count_data[rownames(sig_genes), ]  # Subset DEGs

# Proceed with WGCNA steps from earlier, using:
datExpr <- t(deg_counts)  # Genes as columns

# 1.2 Filter lowly expressed genes
keep_genes <- rowSums(datExpr > 1) >= 0.5*ncol(datExpr)  # Expressed in >50% samples
datExpr <- datExpr[keep_genes, ]

# 1.3 Check for outliers
sampleTree <- hclust(dist(datExpr), method = "average")
plot(sampleTree, main = "Sample clustering")

# Remove outliers (if any, adjust height cutoff)
dynamic_clusters <- cutreeDynamic(sampleTree, method = "tree", minClusterSize = 10)
table(dynamic_clusters)  # Check distribution of samples in clusters

keepSamples <- (dynamic_clusters == 1)
datExpr <- datExpr[keepSamples, ]

# ----------------------------
# Soft Thresholding (Adjusted for DEGs)
# ----------------------------

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

# ----------------------------
#  Network Construction
# ----------------------------

# One-step network construction (auto-detects blocks)
net <- blockwiseModules(
  datExpr,
  power = softPower,
  networkType = "signed",  # Emphasize positive correlations
  TOMType = "signed",
  minModuleSize = 20,      # Smaller modules for DEGs
  mergeCutHeight = 0.20,   # More conservative merging
  deepSplit = 2,           # Balanced splitting
  pamRespectsDendro = FALSE,
  verbose = 3
)

# Module colors and labels
moduleColors <- net$colors
table(moduleColors)  # Check module sizes

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

pdf("dendrogram.pdf", width = 10, height = 6)
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

# ----------------------------
# Module-Trait Associations
# ----------------------------

# Prepare factors
traitData <- metadata[rownames(datExpr), ] %>%
  mutate(
    CD_status = as.numeric(Patient_group == "Israel_CD"),
    Male = as.numeric(Gender == "male"),
    # Convert to numeric, handling "na" strings
    CRP_mg_L = as.numeric(ifelse(CRP_mg_L == "na", NA, CRP_mg_L)),
    Calprotectin_ug_g = as.numeric(ifelse(Calprotectin_ug_g == "na", NA, Calprotectin_ug_g)),
    # Convert age to factor
    Age = factor(Age, levels = c("20_29", "30_39", "40_49", "50_59", "60_69"))
  )

# Create dummy variables for age
age_dummies <- model.matrix(~ Age - 1, data = traitData)
colnames(age_dummies) <- gsub("Age", "Age_", colnames(age_dummies))

# Prepare numeric factors and impute NAs with column means
numeric_traits <- traitData %>%
  select(CD_status, Male, CRP_mg_L, Calprotectin_ug_g) %>%
  mutate(across(where(is.numeric), ~ ifelse(is.na(.x), mean(.x, na.rm = TRUE), .x))) %>%
  as.matrix()

# Combine all factors
traitMatrix <- cbind(numeric_traits, age_dummies)

# Verify all numeric
print("Column types after conversion:")
print(sapply(as.data.frame(traitMatrix), class))

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

# Plot
pdf("module_trait_relationships_final.pdf", width = 12, height = 8)
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
write_csv(module_labels, "gene_module_assignments.csv")

# Save module eigengenes (MEs)
write_csv(as.data.frame(MEs) %>% rownames_to_column("sample_id"), "module_eigengenes.csv")

# Calculate intramodular connectivity
kME <- signedKME(datExpr, MEs)  # Gene-module eigengene correlation
colnames(kME) <- gsub("kME", "", colnames(kME))

# For each module, get top 20 hub genes
hub_genes <- lapply(unique(moduleColors), function(mod) {
  kME_sub <- kME[moduleColors == mod, mod]
  names(sort(kME_sub, decreasing = TRUE))[1:20]
})
names(hub_genes) <- unique(moduleColors)
write_rds(hub_genes, "hub_genes_per_module.rds")

library(clusterProfiler)
library(org.Hs.eg.db)
```

# GO analysis of the module genes:
```
library(clusterProfiler)
library(org.Hs.eg.db)
library(dplyr)

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
                file = paste0(mod, "_GO_enrichment.csv"))
    }
    # Save KEGG results
    if (!is.null(enrichment_results[[mod]]$KEGG)) {
      write.csv(as.data.frame(enrichment_results[[mod]]$KEGG),
                file = paste0(mod, "_KEGG_enrichment.csv"))
    }
  }
}
```

# Build plots
```
library(clusterProfiler)
library(ggplot2)

# Custom plotting function
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
pdf("module_enrichment_plots_robust.pdf", width = 12, height = 8)

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

### Recalculate hub-genes
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
  # Intra-modular connectivity
  mod_connectivity <- kIN[mod_genes, "kWithin"]
  
  # Module membership (kME)
  mod_kME <- cor(datExpr[, mod_genes], MEs[, paste0("ME", module)])
  
  # Gene significance (correlation with traits)
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
# First, ensure you have the plottable data
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
pdf("hub_gene_validation_plots.pdf", width = 10, height = 8)
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
### Remove immunoglobulin genes to prepare data for a BN
```
## Remove green module
library(dplyr)

# Remove green module genes
hub_gene_df_filtered <- hub_gene_df %>%
  filter(module != "green")

# Verify removal
table(hub_gene_df_filtered$module)

# Filter out immunoglobulin genes and solute carriers
filtered_hub_genes <- hub_gene_df_filtered %>%
  filter(!grepl("^(IG|IGH|IGK|IGL|SLC)", gene))
write.csv(filtered_hub_genes, "WGCNA_hub_genes_for_Bayesian_network.csv", row.names = FALSE)
```
