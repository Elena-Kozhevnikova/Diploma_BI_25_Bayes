## This code tests the Bayesian Network's predictive power
It uses the code the results of previous steps

## Prepare data
```
library(bnlearn)
library(Rgraphviz)
library(dplyr)
library(tidyr)
library(tibble)
library(igraph)
library(ggplot2)
library(ggpubr)
library(scales)
library(lmtest)
library(tidyverse)

# Read the hub gene list
hub_genes <- read.csv("WGCNA_hub_genes_for_Bayesian_network_42_samples.csv")

# Read expression data
datExpr <- readRDS("datExpr.rds")

# Extract expression data for hub genes only
hub_expr <- datExpr[, hub_genes$gene]

# Subset only the hub genes (columns) and keep samples as rows
bn_data <- as.data.frame(datExpr[, colnames(datExpr) %in% hub_genes$gene])

# Download metadata
metadata <- read.csv("Israel_metadata.csv", row.names = "pn_ID")
metadata$Calprotectin_ug_g <- as.numeric(metadata$Calprotectin_ug_g)
metadata <- metadata %>% rownames_to_column("sample_id")

# Calculate group-wise medians for each biomarker
group_medians <- metadata %>%
  group_by(Patient_group) %>%
  summarise(
    Calprotectin_median = median(Calprotectin_ug_g, na.rm = TRUE)
  )

# Replace NAs with group medians
metadata <- metadata %>%
  right_join(group_medians, by = "Patient_group") %>%
  mutate(
    Calprotectin_ug_g = ifelse(is.na(Calprotectin_ug_g), 
                               Calprotectin_median, Calprotectin_ug_g))

# Remove temp columns
metadata <- metadata[, !names(metadata) %in% "Calprotectin_median"]

#  Verify sample IDs exist in both datasets
shared_samples <- intersect(rownames(bn_data), metadata$sample_id)

# Subset both datasets to matching samples
bn_data_aligned <- bn_data[shared_samples, ]
metadata_aligned <- metadata[metadata$sample_id %in% shared_samples, ]

# Convert rownames to column for merging
temp_df <- data.frame(sample_id = rownames(bn_data_aligned))
metadata_aligned <- merge(temp_df, metadata_aligned, by = "sample_id", all.x = TRUE)

# Reset rownames
rownames(metadata_aligned) <- metadata_aligned$sample_id
metadata_aligned <- metadata_aligned[rownames(bn_data_aligned), ]  # Final ordering

# Merge
bn_data_with_calprotectin <- cbind(
  bn_data_aligned,
  Calprotectin_ug_g = metadata_aligned$Calprotectin_ug_g
)

# Save BN dataframe with calprotectin
write.csv(bn_data_with_disease_state, "Israel_training_data_cal.csv")

bn_data_with_disease_state <- cbind(
  bn_data_aligned,
  disease = metadata_aligned$Patient_group
)

# Save BN dataframe with disease (need it later)
write.csv(bn_data_with_disease_state, "Israel_training_data_dis.csv")
```
## Build BN with calprotectin data
```

# Define whitelist to force Calprotectin as terminal node
whitelist <- data.frame(
  from = setdiff(colnames(bn_data_with_calprotectin), "Calprotectin_ug_g"),
  to = "Calprotectin_ug_g"
)

# Learn initial structure with constraints
bn_model_calprotectin <- hc(
  bn_data_with_calprotectin,
  whitelist = whitelist,
  score = "bic-g"
)

# Bootstrapping with R replicates
boot_results_calprotectin <- boot.strength(
  data = bn_data_with_calprotectin,
  R = 100,
  algorithm = "hc",
  algorithm.args = list(
    #whitelist = whitelist,
    score = "bic-g"
  )
)

saveRDS(boot_results_calprotectin, file = "BN_calprotectin.rds")

# Filter for significant edges (threshold = 0.85)
significant_edges_calprotectin <- boot_results_calprotectin[
  boot_results_calprotectin$strength > 0.7 
  & boot_results_calprotectin$direction > 0.5, ]

print(significant_edges_calprotectin)

# Build consensus network
consensus_bn <- averaged.network(boot_results_calprotectin, threshold = 0.5)

# Verify calprotectin has incoming edges (parents)
consensus_bn$nodes$Calprotectin_ug_g$parents
```

## Identify calprotectin parents
```

parents <- consensus_bn$nodes$Calprotectin_ug_g$parents
print(parents)

library(dplyr)

# Filter boot_results for edges pointing to Calprotectin
calprotectin_edges <- boot_results_calprotectin %>%
  filter(to == "Calprotectin_ug_g") %>%
  # Keep only edges from the actual parents
  filter(from %in% parents) %>%
  # Sort by strength (descending)
  arrange(desc(strength))

# Print the results
print(calprotectin_edges)
```
## Build BN with calprotectin

```
bn_igraph_calprotectin <- as.igraph(consensus_bn)
darjeeling_colors <- wes_palette("Moonrise3")

module_palette <- c(
  "yellow" = darjeeling_colors[1],
  "blue" = darjeeling_colors[2],
  "turquoise" = darjeeling_colors[3],
  "brown" = darjeeling_colors[5],
  "green" = darjeeling_colors[4]
  
)

node_names <- V(bn_igraph_calprotectin)$name
modules_for_nodes <- setNames(rep(NA, length(node_names)), node_names)

modules_for_nodes[names(modules_for_nodes) %in% hub_genes$gene] <- 
  hub_genes$module[match(names(modules_for_nodes)[names(modules_for_nodes) %in% 
                                                      hub_genes$gene], hub_genes$gene)]

node_degrees <- degree(bn_igraph_calprotectin)
label_colors <- ifelse(node_degrees >= 0.7, "black", "grey")
V(bn_igraph_calprotectin)$color <- module_palette[modules_for_nodes[V(bn_igraph_calprotectin)$name]]

# Теперь визуализируем граф с цветами
pdf("network_colored_calprotectin_bn_0805a.pdf", width = 40, height = 40)
plot(bn_igraph_calprotectin,
     vertex.size = rescale(degree(bn_igraph_calprotectin), to = c(10, 20)),
     vertex.label.color = label_colors,  # Color by influence
     vertex.size = scales::rescale(degree(bn_igraph_calprotectin), to = c(8, 15)),
     vertex.label = V(bn_igraph_calprotectin)$name,
     vertex.label.cex = 3.0, 
     vertex.color = V(bn_igraph_calprotectin)$color,
     edge.arrow.size = 2,
     edge.width = 1.5,  # ширина ребер
     layout = layout_with_dh)

legend("bottomleft",
       legend = paste("Strong edges (n =", ecount(bn_igraph_calprotectin), ")"),
       col = "red", lwd = 3, bty = "n")

dev.off()
```

## Fit the model to predict calprotectin

```
# Check if the graph remained acyclic after bootstrapping
is_dag(bn_igraph_calprotectin)

# Remove the cycles
consensus_dag <- cextend(consensus_bn_file)

# Test again
# consensus_igraph_calprotectin <- as.igraph(consensus_dag)
# is_dag(consensus_igraph_calprotectin)

# fit the graph model
fitted_bn_calprotectin <- bn.fit(consensus_dag, data = bn_data_with_calprotectin)
print(fitted_bn_calprotectin)

# Remove calprotectin column to create "new" data
new_data_calprotectin <- bn_data_with_calprotectin[, !colnames(bn_data_with_calprotectin) %in% "Calprotectin_ug_g"]

# Predict using Bayesian weighting
predictions <- predict(
  fitted_bn_calprotectin,
  node = "Calprotectin_ug_g",
  data = new_data_calprotectin,
  method = "bayes-lw"
)

# Save the fit
saveRDS(fitted_bn_calprotectin, "fitted_bn_calprotectin.rds")
```

## Test fitted model
```

# Compare to observed values

plot_data <- data.frame(
  Observed = bn_data_with_calprotectin$Calprotectin_ug_g,
  Predicted = predictions
)

# Calculate correlation stats
cor_test <- cor.test(plot_data$Observed, plot_data$Predicted, method = "pearson")
r_value <- round(cor_test$estimate, 3)
p_value <- ifelse(cor_test$p.value < 0.001, "< 0.001", 
                  round(cor_test$p.value, 3))

# Create the plot
ggplot(plot_data, aes(x = Observed, y = Predicted)) +
  geom_point(aes(color = abs(Observed - Predicted)),
             size = 3, alpha = 0.8) +
  geom_smooth(method = "lm", se = TRUE, color = "darkred", fill = "lightpink") +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "grey40") +
  
  # Regression equation and stats
  stat_regline_equation(label.x = 0.1, label.y = 0.9 * max(plot_data$Predicted),
                        aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~~"))) +
  
  # Customize appearance
  scale_color_gradientn(name = "Error", 
                        colors = c("blue", "yellow", "red"),
                        values = rescale(c(0, median(abs(plot_data$Observed - plot_data$Predicted)), 
                                           max(abs(plot_data$Observed - plot_data$Predicted))))) +
  labs(title = "Observed vs Predicted Calprotectin Levels",
       subtitle = paste("Pearson r =", r_value, " (p ", p_value, ")", sep = ""),
       x = "Observed Calprotectin (μg/g)",
       y = "Predicted Calprotectin (μg/g)") +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "right",
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    panel.grid.minor = element_blank()
  ) +
  coord_equal()

# Save the graph
ggsave("calprotectin_prediction_correlation_0405.pdf", width = 8, height = 7, dpi = 300)
```

## Test linnear dependancies
```

# Fit linear and nonlinear models
linear <- lm(Predicted ~ Observed, data = plot_data)
nonlinear <- lm(Predicted ~ poly(Observed, 2), data = plot_data) # Quadratic

# Compare models
anova(linear, nonlinear)  # Significant p-value indicates nonlinearity

# Extract regression coefficients
coefs <- fitted_bn_calprotectin$Calprotectin_ug_g$coefficients
sd_effects <- coefs[-1] * sapply(bn_data_with_calprotectin[, parents, drop=FALSE], sd)
data.frame(
  Gene = parents,
  Coefficient = coefs[-1],  # Remove intercept
  SD_Effect = sd_effects,
  Relative_Influence = abs(sd_effects) / sum(abs(sd_effects))
) %>% arrange(desc(Relative_Influence))

coef_df <- data.frame(
  Gene = factor(parents, levels = parents[order(abs(sd_effects))]),
  Effect = sd_effects
)

effect_size <- ggplot(coef_df, aes(x = Gene, y = Effect, fill = Effect > 0)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("red", "blue")) +
  coord_flip() +
  labs(title = "Gene Effects on Calprotectin Prediction",
       y = "Effect Size (SD units change in calprotectin per SD increase in gene)")
print(effect_Size)

ggsave(
  filename = "calprotectin_parent_effect_size.pdf",
  plot = effect_size,
  width = 8,      # Width in inches (adjust as needed)
  height = 6,     # Height in inches
  dpi = 300,      # High resolution for publications
  bg = "white",    # Background color
  device = "pdf"
)

simple_lm <- lm(Calprotectin_ug_g ~ ISX + ACTN1 + CECR1 + IFI16 + TMEM120A , data = bn_data_with_calprotectin)
summary(simple_lm)
```

## Apply BN predictor to another set of CD data
```
#Load data
metadata_china <- read_csv("china_metadata.csv") %>% as.data.frame()
expression_data <- read_csv("modified_output.csv") %>% as.data.frame()

# Get intersecting IDs
common_ids <- intersect(metadata_china$pn_ID, expression_data$id)
merged_data <- metadata_china %>% 
  filter(pn_ID %in% common_ids) %>% 
  left_join(
    expression_data %>% 
      filter(id %in% common_ids) %>% 
      group_by(id, gene_name) %>% 
      summarise(tpm = sum(tpm), .groups = 'drop') %>% 
      pivot_wider(names_from = gene_name, values_from = tpm, values_fill = 0),
    by = c("pn_ID" = "id")
  ) %>% 
  as.data.frame()

# Identify genes in BN but missing in new data
required_genes <- nodes(fitted_bn_calprotectin)
missing_genes <- setdiff(required_genes, colnames(merged_data))

# Select and order columns to match BN
prediction_data <- merged_data[, nodes(fitted_bn_calprotectin), drop = FALSE]
prediction_data <- prediction_data[, setdiff(colnames(prediction_data), "Calprotectin_ug_g")]

# Make predictions
predictions <- predict(
  fitted_bn_calprotectin,
  node = "Calprotectin_ug_g",
  data = prediction_data,
  method = "bayes-lw"
)

# Combine with metadata
results <- merged_data %>%
  dplyr::select(pn_ID, Patient_group, Age, Gender, CRP_mg_L) %>%
  dplyr::mutate(Predicted_Calprotectin = predictions)

china_control_calprotectin <- results$Predicted_Calprotectin[results$Patient_group == "China_urban"]
china_CD_calprotectin <- results$Predicted_Calprotectin[results$Patient_group == "China_CD"]

wilcox_test_result <- wilcox.test(china_control_calprotectin, china_CD_calprotectin)

print(wilcox_test_result)

my_fill_colors <- c("#E69F00", "#56B4E9", "#009E73")

# Create the boxplot filled with custom colors and black-outlined jittered points
my_plot <- ggplot(results, aes(x = Patient_group, y = Predicted_Calprotectin)) +
  geom_boxplot(aes(fill = Patient_group), outlier.shape = NA) + # Fill boxplots with group colors
  geom_jitter(aes(fill = Patient_group), shape = 21, color = "black", width = 0.2, alpha = 0.7, size = 3) + # Fill points with group colors, outline black
  scale_fill_manual(values = my_fill_colors) + # Apply your custom fill colors
  labs(
    title = "Predicted Calprotectin by Patient Group",
    x = "Patient Group",
    y = "Predicted Calprotectin",
    fill = "Patient Group" # Label for the fill legend
  )

p_value <- wilcox_test_result$p.value
formatted_p_value <- sprintf("p = %.3f", p_value) # Format to 3 decimal places

annotation_x <- 1.5 # Default to between the first two positions if indexes aren't found
annotation_y <- max(results$Predicted_Calprotectin, na.rm = TRUE) * 0.95

final_plot <- my_plot +
  annotate("text",
           x = annotation_x,
           y = annotation_y,
           label = formatted_p_value,
           hjust = 0.5, # Center the text horizontally
           vjust = 1) # Align the bottom of the text with the y-coordinate

# Print the final plot
print(final_plot)

# Save high-quality version
ggsave("predicted_calprotectin_by_group.pdf", width = 10, height = 6, dpi = 300)
```

## Predict disease state
```
# Here work with this dataset: "bn_data_with_disease_state" 

# Convert Patient Group to a Discrete Node

bn_data_with_disease_state$disease <- as.factor(bn_data_with_disease_state$disease)

# Remove the gene not present in my test dataframe
bn_data_with_disease_state <- bn_data_with_disease_state %>%
  select(-"DHRS4.AS1")

# Identify continuous variables (assuming all except disease are continuous)
continuous_vars <- names(which(sapply(bn_data_with_disease_state, is.numeric)))

# Discretize continuous variables into 3 bins (Low/Medium/High)
bn_orig_data_discrete <- bn_data_with_disease_state %>%
  mutate(across(
    all_of(continuous_vars),
    ~ cut(.x, 
          breaks = 3, 
          labels = c("Low", "Medium", "High"),
          include.lowest = TRUE)
  ))

bn_orig_data_discrete <- bn_orig_data_discrete[, -1]

# Verify structure
str(bn_orig_data_discrete)

# Blacklist: Block all edges FROM disease
blacklist <- data.frame(
  from = "disease",
  to = setdiff(colnames(bn_orig_data_discrete), "disease")
)

# Learn structure with constraints
bn_model_orig_discrete <- hc(
  bn_orig_data_discrete,
  blacklist = blacklist,
  score = "aic"
)

# Verify disease is terminal (no children)
stopifnot(length(bn_model_orig_discrete$nodes$disease$children) == 0)


# Get edge strengths with 100 bootstrap replicates
boot_strengths <- boot.strength(
  bn_orig_data_discrete,
  algorithm = "hc",
  algorithm.args = list(
    blacklist = blacklist,
    score = "aic"
  ),
  R = 100
)

# Refine the network using averaged.network with a strength threshold of 0.7
refined_bn_structure <- averaged.network(boot_strengths, threshold = 0.7)

consensus_dis <- cextend(refined_bn_structure)

fitted_bn_dis <- bn.fit(consensus_dis, data = bn_orig_data_discrete)
print(fitted_bn_dis)

# Remove calprotectin column to simulate "new" data
new_data_dis <- bn_orig_data_discrete[, !colnames(bn_orig_data_discrete) %in% "disease"]

# Predict using Bayesian weighting
predictions <- predict(
  fitted_bn_dis,
  node = "disease",
  data = new_data_dis,
  method = "bayes-lw"
)

library(caret)

confusion_matrix <- confusionMatrix(bn_orig_data_discrete$disease, predictions)

# Print the confusion matrix and associated statistics
print(bn_orig_data_discrete$disease)
```

## Test the predictor BN using another set of data
```
# Filter edges pointing to disease (terminal node)
disease_edges <- boot_strengths %>%
  filter(to == "disease") %>%
  arrange(desc(strength))

print(disease_edges)

bn_model_orig_discrete$nodes$disease$parents

library(igraph)

bn_igraph_disease <- as.igraph(bn_model)

module_palette <- c(
  "yellow" = "#FFB3BA",
  "blue" = "#BAE1FF",
  "turquoise" = "#FFFFBA",
  "brown" = "#B5EAD7",
  "green" = "grey"
 )
node_names <- V(bn_igraph_disease)$name
modules_for_nodes <- setNames(rep(NA, length(node_names)), node_names)

modules_for_nodes[names(modules_for_nodes) %in% hub_genes$gene] <- 
  hub_genes$module[match(names(modules_for_nodes)[names(modules_for_nodes) %in% 
                                                    hub_genes$gene], hub_genes$gene)]

node_degrees <- degree(bn_igraph_disease)
label_colors <- ifelse(node_degrees >= 0.7, "black", "grey")
V(bn_igraph_disease)$color <- module_palette[modules_for_nodes[V(bn_igraph_disease)$name]]

# Build the graph
pdf("network_colored_calprotectin_bn_0505disc.pdf", width = 40, height = 40)
plot(bn_igraph_disease,
     vertex.size = rescale(degree(bn_igraph_disease), to = c(10, 20)),
     vertex.label.color = label_colors,  # Color by influence
     vertex.size = scales::rescale(degree(bn_igraph_disease), to = c(8, 15)),
     vertex.label = V(bn_igraph_disease)$name,
     vertex.label.cex = 3.0, 
     vertex.color = V(bn_igraph_disease)$color,
     edge.arrow.size = 2,
     edge.width = 1.5,  # ширина ребер
     layout = layout_with_fr)

legend("bottomleft",
       legend = paste("Strong edges (n =", ecount(bn_igraph_disease), ")"),
       col = "red", lwd = 3, bty = "n")

dev.off()
```

## Do the fit
```

# Refit the BN with Patient_group as a node
fitted_bn_disease_disc <- bn.fit(bn_model_orig_discrete, data = bn_orig_data_discrete, method = "bayes")  # For discrete nodes

metadata_china <- read_csv("china_metadata.csv") 
metadata_china <- metadata_china %>%
  rename(disease = Patient_group)

expression_data <- read_csv("modified_output.csv") %>% as.data.frame()

# Get intersecting IDs
common_ids <- intersect(metadata_china$pn_ID, expression_data$id)
merged_data <- metadata_china %>% 
  filter(pn_ID %in% common_ids) %>% 
  left_join(
    expression_data %>% 
      filter(id %in% common_ids) %>% 
      group_by(id, gene_name) %>% 
      summarise(tpm = sum(tpm), .groups = 'drop') %>% 
      pivot_wider(names_from = gene_name, values_from = tpm, values_fill = 0),
    by = c("pn_ID" = "id")
  ) %>% 
  as.data.frame()

# Use only genes from the network
gene_names <- nodes(bn_model_orig_discrete)
gene_names <- gene_names[gene_names != "disease"]

new_data_prepared <- merged_data %>%
  select(all_of(gene_names), disease)  # Keep only BN-relevant columns 

# Rename the patient group from Chinese to Israel
new_data_prepared <- new_data_prepared %>%
  mutate(disease = case_when(
    disease == "China_urban" ~ "Israel_control",
    disease == "China_rural" ~ "Israel_control",
    disease == "China_CD"    ~ "Israel_CD",
    TRUE                      ~ disease # Keep other values as they are
  ))

# Discretize continuous variables into 3 bins (Low/Medium/High)
new_data_prepared_discrete <- new_data_prepared %>%
  mutate(across(
    all_of(continuous_vars),
    ~ cut(.x, 
          breaks = 3, 
          labels = c("Low", "Medium", "High"),
          include.lowest = TRUE)
  ))

new_data_prepared_discrete$disease <- as.factor(new_data_prepared_discrete$disease)
new_data_prepared_discrete <- as.data.frame(new_data_prepared_discrete)

new_data_for_prediction <- new_data_prepared_discrete[, !(names(new_data_prepared_discrete) %in% "disease")]

predictions <- predict(
  fitted_bn_disease_disc,
  node = "disease",
  data = new_data_for_prediction,
  method = "bayes-lw"
)

print(predictions)

# Create the confusion matrix

confusion_matrix <- confusionMatrix(new_data_prepared_discrete$disease, predictions)

# Print the confusion matrix and associated statistics
print(confusion_matrix)

# Create the confusion matrix
conf_matrix <- matrix(c(12, 1, 8, 19), 
                      nrow = 2,
                      byrow = TRUE,
                      dimnames = list(Prediction = c("Israel_CD", "Israel_control"),
                                      Reference = c("Israel_CD", "Israel_control")))

# Melt the matrix for ggplot
melted_matrix <- melt(conf_matrix)
colnames(melted_matrix) <- c("Prediction", "Reference", "value")

# Reorder the factor levels to put Israel_CD first (top)
melted_matrix$Prediction <- factor(melted_matrix$Prediction, 
                                   levels = rev(c("Israel_CD", "Israel_control")))
melted_matrix$Reference <- factor(melted_matrix$Reference,
                                  levels = c("Israel_CD", "Israel_control"))

# Create the heatmap
my_plot <- ggplot(melted_matrix, aes(x = Reference, y = Prediction, fill = value)) +
  geom_tile(color = "white") +
  geom_text(aes(label = value), color = "black", size = 6) +
  scale_fill_gradient(low = "#00A08A", high = "#FF0000") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 0, vjust = 0.5),
        panel.grid = element_blank()) +
  labs(title = "Confusion Matrix Heatmap",
       x = "Reference", y = "Prediction") +
  coord_fixed()

ggsave("confusion_matrix_heatmap.pdf", my_plot, width = 5, height = 5)
```
