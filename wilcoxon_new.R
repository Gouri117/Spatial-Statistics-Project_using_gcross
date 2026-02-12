# ------------------ Load Libraries ------------------
library(dplyr)
library(readr)
library(ggpubr)
library(rstatix)

# ------------------ Load Data ------------------
df <- read_csv("/nfs/turbo/umms-ukarvind/shared_data/mfIHC_collab_data/stanford_crc/results/cn_sample_matrix_k20_with_outcome.csv")

# ------------------ Prepare ------------------
df <- df %>% filter(!is.na(primary_outcome))
df$primary_outcome <- as.factor(df$primary_outcome)

# Gather CNs into long format
df_long <- df %>%
  pivot_longer(
    cols = -c(image, primary_outcome),
    names_to = "CN",
    values_to = "cn_proportion"
  )

# ------------------ Wilcoxon Rank-Sum Test ------------------
stat.test <- df_long %>%
  group_by(CN) %>%
  wilcox_test(cn_proportion ~ primary_outcome, detailed = TRUE) %>%
  adjust_pvalue(method = "fdr") %>%
  add_significance("p.adj")

# ------------------ Significant CNs ------------------
signif_pairs_test <- stat.test %>% filter(p.adj.signif != "ns")
signif_cns <- signif_pairs_test$CN

cat(length(signif_cns), "significant CNs found\n")
results <- lapply(colnames(cn_matrix), function(cn) {
  fit <- glm(primary_outcome ~ cn_matrix[[cn]], family = "binomial")
  summary(fit)$coefficients[2, c("Estimate", "Pr(>|z|)")]
})
