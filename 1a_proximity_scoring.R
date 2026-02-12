####################################################################################
## overview: script to perform CN annotation and proximity scoring
## created: 7/29/25

## script inputs:
## 1. processed data file: /nfs/turbo/umms-ukarvind/shared_data/mfIHC_collab_data/stanford_crc/results/g_cross_clustered_k20.csv

## script outputs:
## 1. 
## 2. 
####################################################################################

# Load necessary libraies
# Load necessary libraries
library(dplyr)
library(readr)
library(tidyr)
library(RANN)
library(rstatix)
library(dendextend)
library(pheatmap)

# Load the data
# Read the CSV file
df <- read_csv("/nfs/turbo/umms-ukarvind/shared_data/mfIHC_collab_data/stanford_crc/results/g_cross_clustered_k20.csv")

df <- df %>% filter(CELL_TYPE != "Other")

# Quick check
print(dim(df))
print(colnames(df))
head(df)

# Make sure cluster is numeric or factor
df <- df %>%
  mutate(
    cluster = as.integer(cluster_k20),
    X = as.numeric(X),
    Y = as.numeric(Y)
  )

# Extract Tumor cells
# Get coordinates of Individual Tumor cells
tumor_cells <- df %>%
  filter(CELL_TYPE == "Tumor") %>%
  select(X, Y)

cat("Number of tumor cells:", nrow(tumor_cells), "\n")
head(tumor_cells)

#----------------------------------------CN Annotation------------------------------
# 1. Identify AUC Columns (excluding 'Other')
auc_cols <- grep("→", colnames(df), value = TRUE)
auc_cols <- auc_cols[!grepl("Other", auc_cols)]
cat("Number of AUC columns:", length(auc_cols), "\n")

# 2. Define thresholds for enrichment
thresholds <- df %>%
  summarise(across(all_of(auc_cols), ~ quantile(.x, 0.75, na.rm = TRUE)))

# 3. Calculate interaction enrichment per CN
interaction_enrichment <- df %>%
  group_by(cluster) %>%
  summarise(across(all_of(auc_cols), 
                   ~ sum(.x > thresholds[[cur_column()]], na.rm = TRUE) / n(),
                   .names = "{.col}")) %>%
  ungroup()

# 4. Cell type proportions per CN
cluster_celltype_prop <- df %>%
  group_by(cluster, CELL_TYPE) %>%
  summarise(count = n(), .groups = "drop_last") %>%
  mutate(freq = count / sum(count)) %>%
  ungroup()

# 5. Top 3 cell types per CN with percentages
top_celltypes <- cluster_celltype_prop %>%
  group_by(cluster) %>%
  arrange(desc(freq)) %>%
  slice_head(n = 3) %>%
  summarise(top_cells = paste0(CELL_TYPE, "(", round(freq * 100, 1), "%)", collapse = ", "),
            .groups = "drop")

# 6. Top 5 enriched interactions per CN
top_interactions <- interaction_enrichment %>%
  pivot_longer(cols = all_of(auc_cols), names_to = "Interaction", values_to = "Enrichment") %>%
  group_by(cluster) %>%
  arrange(desc(Enrichment)) %>%
  slice_head(n = 5) %>%
  summarise(enriched_interactions = paste(Interaction, collapse = ", "), .groups = "drop")

# 7. Final CN annotation (merged only once)
cluster_labels <- top_celltypes %>%
  inner_join(top_interactions, by = "cluster") %>%
  mutate(CN_Label = paste0("CN (", top_cells, ") | Enriched: ", enriched_interactions))

# View result
head(cluster_labels)


#-----------------------------------------Tumor proximity scoring------------------------------------------------
# Compute raw proximity
tumor_cells <- df %>%
  filter(CELL_TYPE == "Tumor") %>%
  select(X, Y)

distances <- nn2(
  data = tumor_cells, 
  query = df %>% select(X, Y), 
  k = 1
)$nn.dists

df$dist_to_tumor <- distances

# Median proximity per patient x CN
patient_cluster_dist <- df %>%
  group_by(image, cluster_k20) %>%
  summarise(median_distance = median(dist_to_tumor, na.rm = TRUE), .groups = "drop")

# Tumor proportion per patient x CN
patient_cluster_tumor_prop <- df %>%
  group_by(image, cluster_k20) %>%
  summarise(
    tumor_count = sum(CELL_TYPE == "Tumor"),
    total_count = n(),
    tumor_prop = tumor_count / total_count,
    .groups = "drop"
  )

# Combine proximity and tumor proprotion
patient_cluster_metrics <- patient_cluster_dist %>%
  left_join(patient_cluster_tumor_prop %>% select(image, cluster_k20, tumor_prop),
            by = c("image", "cluster_k20"))

# Calculate weighted proximity score
patient_cluster_metrics <- patient_cluster_metrics %>%
  mutate(weighted_proximity = (1 - tumor_prop) * median_distance)

# Patient frequency scoring
# Calculate frequency of each CN per patient
patient_cluster_freq <- df %>%
  group_by(image, cluster_k20) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(image) %>%
  mutate(freq = count / sum(count)) %>%
  ungroup()

# --------------------------------------------Significance analysis---------------------------------------
# Merge metrics
metrics_with_outcome <- patient_cluster_metrics %>%
  left_join(patient_cluster_freq %>% select(image, cluster_k20, freq),
            by = c("image", "cluster_k20")) %>%
  left_join(
    df %>%
      select(image, primary_outcome) %>%
      distinct(),
    by = "image"
  )

# Statistical Test
results <- metrics_with_outcome %>%
  group_by(cluster_k20) %>%
  group_modify(~ {
    data <- as.data.frame(.x)  # ensure it's a data.frame
    
    # handle cases where all values are same or only one outcome group
    freq_p <- if(length(unique(data$primary_outcome)) > 1) {
      wilcox.test(data$freq[data$primary_outcome == 0],
                  data$freq[data$primary_outcome == 1])$p.value
    } else { NA }
    
    prox_p <- if(length(unique(data$primary_outcome)) > 1) {
      wilcox.test(data$weighted_proximity[data$primary_outcome == 0],
                  data$weighted_proximity[data$primary_outcome == 1])$p.value
    } else { NA }
    
    tibble(
      freq_pval = freq_p,
      proximity_pval = prox_p,
      tumor_prop_good = median(data$tumor_prop[data$primary_outcome == 0], na.rm = TRUE),
      tumor_prop_poor = median(data$tumor_prop[data$primary_outcome == 1], na.rm = TRUE)
    )
  }) %>%
  ungroup()

# CNs sorted by p-value for frequency
top_freq_cns <- results %>%
  arrange(freq_pval) %>%
  slice_head(n = 3) %>%
  pull(cluster_k20)

# CNs sorted by p-value for weighted proximity
top_prox_cns <- results %>%
  arrange(proximity_pval) %>%
  slice_head(n = 3) %>%
  pull(cluster_k20)

library(ggplot2)

metrics_with_outcome %>%
  filter(cluster_k20 %in% top_freq_cns) %>%
  ggplot(aes(x = factor(primary_outcome), y = freq, fill = factor(primary_outcome))) +
  geom_boxplot() +
  facet_wrap(~ cluster_k20, scales = "free_y") +
  labs(
    title = "Top CNs with Significant Frequency Differences",
    x = "Outcome (0 = good, 1 = poor)",
    y = "CN Frequency"
  ) +
  theme_minimal()


metrics_with_outcome %>%
  filter(cluster_k20 %in% top_prox_cns) %>%
  ggplot(aes(x = factor(primary_outcome), y = weighted_proximity, fill = factor(primary_outcome))) +
  geom_boxplot() +
  facet_wrap(~ cluster_k20, scales = "free_y") +
  labs(
    title = "Top CNs with Significant Weighted Proximity Differences",
    x = "Outcome (0 = good, 1 = poor)",
    y = "Weighted Proximity to Tumor"
  ) +
  theme_minimal()

# Identify significant CNs (from Wilcoxon tests)
sig_cns <- unique(c(
  results$cluster_k20[results$freq_pval < 0.05],
  results$cluster_k20[results$proximity_pval < 0.05]
))

#---------------------------Reducing redundancy---------------------------------------------------
# --- 1. Create CN × Cell Type composition matrix ---
cn_matrix <- df %>%
  group_by(cluster_k20, CELL_TYPE) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(cluster_k20) %>%
  mutate(freq = count / sum(count)) %>%
  select(cluster_k20, CELL_TYPE, freq) %>%
  pivot_wider(names_from = CELL_TYPE, values_from = freq, values_fill = 0)

# --- 2. Hierarchical clustering and cut tree at height 2 ---
dist_matrix <- dist(cn_matrix %>% select(-cluster_k20))
hc <- hclust(dist_matrix, method = "ward.D2")

plot(hc)  # visualize dendrogram
# Cut at height = 2
meta_cn_assignments <- cutree(hc, h = 2)

# Create mapping table
meta_cn_df <- data.frame(
  cluster_k20 = cn_matrix$cluster_k20,
  meta_cn = meta_cn_assignments
)

# --- 3. Add meta-CNs to main dataframe ---
df_merged <- df %>%
  left_join(meta_cn_df, by = "cluster_k20")

# --- 4. Visualize merged CN composition ---
cn_matrix$meta_cn <- meta_cn_assignments
cn_matrix_ordered <- cn_matrix %>%
  arrange(meta_cn)

rownames(cn_matrix_ordered) <- paste0(
  "CN", cn_matrix_ordered$cluster_k20,
  " (Meta ", cn_matrix_ordered$meta_cn, ")"
)

pheatmap(
  cn_matrix_ordered %>% select(-cluster_k20, -meta_cn),
  cluster_rows = FALSE,
  show_rownames = TRUE,
  main = "Cell Type Composition of Meta-CNs (Height 2)"
)

# --- 5. Recompute metrics for merged CNs ---
patient_meta_freq <- df_merged %>%
  group_by(image, meta_cn) %>%
  summarise(freq = n() / nrow(df_merged), .groups = "drop")

# Proximity calculation
library(RANN)
tumor_cells <- df_merged %>%
  filter(CELL_TYPE == "Tumor") %>%
  select(X, Y)

distances <- nn2(
  data = tumor_cells,
  query = df_merged %>% select(X, Y),
  k = 1
)$nn.dists

df_merged$dist_to_tumor <- distances

patient_meta_dist <- df_merged %>%
  group_by(image, meta_cn) %>%
  summarise(median_distance = median(dist_to_tumor, na.rm = TRUE), .groups = "drop")

patient_meta_tumor_prop <- df_merged %>%
  group_by(image, meta_cn) %>%
  summarise(
    tumor_count = sum(CELL_TYPE == "Tumor"),
    total_count = n(),
    tumor_prop = tumor_count / total_count,
    .groups = "drop"
  )

patient_meta_metrics <- patient_meta_dist %>%
  left_join(patient_meta_tumor_prop, by = c("image", "meta_cn")) %>%
  mutate(weighted_proximity = (1 - tumor_prop) * median_distance)

# Merge with outcome
metrics_with_outcome <- patient_meta_metrics %>%
  left_join(patient_meta_freq, by = c("image", "meta_cn")) %>%
  left_join(
    df_merged %>% select(image, primary_outcome) %>% distinct(),
    by = "image"
  )

# --- 6. Wilcoxon tests ---
results_meta <- metrics_with_outcome %>%
  group_by(meta_cn) %>%
  group_modify(~ {
    data <- as.data.frame(.x)
    freq_p <- if(length(unique(data$primary_outcome)) > 1) {
      wilcox.test(
        data$freq[data$primary_outcome == 0],
        data$freq[data$primary_outcome == 1]
      )$p.value
    } else { NA }
    
    prox_p <- if(length(unique(data$primary_outcome)) > 1) {
      wilcox.test(
        data$weighted_proximity[data$primary_outcome == 0],
        data$weighted_proximity[data$primary_outcome == 1]
      )$p.value
    } else { NA }
    
    tibble(
      freq_pval = freq_p,
      proximity_pval = prox_p,
      tumor_prop_good = median(data$tumor_prop[data$primary_outcome == 0], na.rm = TRUE),
      tumor_prop_poor = median(data$tumor_prop[data$primary_outcome == 1], na.rm = TRUE)
    )
  }) %>%
  ungroup()

print(results_meta)

sig_meta_cns <- results_meta %>%
  filter(freq_pval < 0.05 | proximity_pval < 0.05)

results_meta <- results_meta %>%
  mutate(freq_effect = tumor_prop_poor - tumor_prop_good)

# --- 7. Updated dendrogram ---
library(dendextend)
dend <- as.dendrogram(hc)
dend_colored <- color_branches(dend, k = length(unique(meta_cn_assignments)))
labels_colors(dend_colored) <- meta_cn_assignments[order.dendrogram(dend)]

plot(dend_colored, 
     main = "Updated Hierarchical Clustering with Meta-CNs (Height 2)", 
     ylab = "Height", 
     cex = 0.8)

abline(h = 2, col = "red", lty = 2)

length(unique(meta_cn_assignments))
table(meta_cn_assignments)
library(dplyr)
library(tidyr)

# --- Cell type proportions per meta-CN ---
meta_cn_composition <- df_merged %>%
  group_by(meta_cn, CELL_TYPE) %>%
  summarise(count = n(), .groups = "drop_last") %>%
  mutate(freq = count / sum(count)) %>%
  ungroup()

# --- Top 3 cell types per meta-CN ---
top_celltypes_meta <- meta_cn_composition %>%
  group_by(meta_cn) %>%
  arrange(desc(freq)) %>%
  slice_head(n = 3) %>%
  summarise(top_cells = paste0(CELL_TYPE, "(", round(freq * 100, 1), "%)", collapse = ", "),
            .groups = "drop")

# --- Interaction enrichment per meta-CN ---
auc_cols <- grep("→", colnames(df_merged), value = TRUE)
thresholds_meta <- df_merged %>%
  summarise(across(all_of(auc_cols), ~ quantile(.x, 0.75, na.rm = TRUE)))

interaction_enrichment_meta <- df_merged %>%
  group_by(meta_cn) %>%
  summarise(across(all_of(auc_cols), 
                   ~ sum(.x > thresholds_meta[[cur_column()]], na.rm = TRUE) / n(),
                   .names = "{.col}")) %>%
  ungroup()

top_interactions_meta <- interaction_enrichment_meta %>%
  pivot_longer(cols = all_of(auc_cols), names_to = "Interaction", values_to = "Enrichment") %>%
  group_by(meta_cn) %>%
  arrange(desc(Enrichment)) %>%
  slice_head(n = 5) %>%
  summarise(enriched_interactions = paste(Interaction, collapse = ", "), .groups = "drop")

# --- Final annotated labels ---
meta_cn_labels <- top_celltypes_meta %>%
  inner_join(top_interactions_meta, by = "meta_cn") %>%
  mutate(CN_Label = paste0("Meta-CN(", meta_cn, ") | Top Cells: ", top_cells,
                           " | Enriched: ", enriched_interactions))

print(meta_cn_labels)

library(ggplot2)
library(scales)

# --- Dynamic colors ---
n_meta_cn <- length(unique(df_merged$meta_cn))
n_cell_types <- length(unique(df_merged$CELL_TYPE))

color_meta_cn <- scale_color_manual(values = hue_pal()(n_meta_cn))
color_cell_type <- scale_color_manual(values = hue_pal()(n_cell_types))
color_outcome <- scale_fill_manual(values = c("0" = "#1f77b4", "1" = "#ff7f0e"))

# --- Custom theme ---
custom_theme <- theme_minimal(base_size = 14) +
  theme(
    panel.grid = element_blank(),
    strip.text = element_text(face = "bold"),
    axis.title = element_text(face = "bold"),
    plot.title = element_text(face = "bold", hjust = 0.5)
  )

# ---- 1. Meta-CN Frequency Differences ----
ggplot(metrics_with_outcome, aes(x = factor(primary_outcome), y = freq, fill = factor(primary_outcome))) +
  geom_boxplot(alpha = 0.8, outlier.color = "black") +
  facet_wrap(~ meta_cn, scales = "free_y") +
  labs(
    title = "Meta-CN Frequency Differences",
    x = "Outcome (0 = Good, 1 = Poor)",
    y = "Frequency of Meta-CN",
    fill = "Outcome"
  ) +
  color_outcome +
  custom_theme

# ---- 2. Meta-CN Weighted Proximity Differences ----
ggplot(metrics_with_outcome, aes(x = factor(primary_outcome), y = weighted_proximity, fill = factor(primary_outcome))) +
  geom_boxplot(alpha = 0.8, outlier.color = "black") +
  facet_wrap(~ meta_cn, scales = "free_y") +
  labs(
    title = "Meta-CN Weighted Proximity Differences",
    x = "Outcome (0 = Good, 1 = Poor)",
    y = "Weighted Proximity to Tumor",
    fill = "Outcome"
  ) +
  color_outcome +
  custom_theme

# ---- 3. Spatial Localization for Selected Image ----
selected_image <- "Charville_c002_v001_r001_reg056"  # change as needed

ggplot(df_merged %>% filter(image == selected_image), 
       aes(x = X, y = Y, color = factor(meta_cn))) +
  geom_point(alpha = 0.6, size = 0.8) +
  color_meta_cn +
  labs(
    title = paste("Spatial Localization of Meta-CNs for", selected_image),
    color = "Meta-CN"
  ) +
  custom_theme

# ---- 4. Meta-CNs with Tumor Cells Highlighted ----
ggplot(df_merged %>% filter(image == selected_image), aes(x = X, y = Y)) +
  geom_point(aes(color = factor(meta_cn)), alpha = 0.5, size = 0.6) +
  geom_point(
    data = df_merged %>% filter(image == selected_image, CELL_TYPE == "Tumor"),
    color = "black", shape = 4, size = 0.9
  ) +
  color_meta_cn +
  labs(
    title = paste("Meta-CNs with Tumor Cells Highlighted:", selected_image),
    color = "Meta-CN"
  ) +
  custom_theme

# ---- 5. Spatial Distribution of Cell Types ----
ggplot(df_merged %>% filter(image == selected_image), 
       aes(x = X, y = Y, color = CELL_TYPE)) +
  geom_point(alpha = 0.6, size = 0.8) +
  color_cell_type +
  labs(
    title = paste("Spatial Distribution of Cell Types:", selected_image),
    color = "Cell Type"
  ) +
  custom_theme

