################################################################################
# Title: Statistical Analysis of Genotype–Phenotype Associations in ID Genes
# Author: Hussein Kassir
# Description: This script performs statistical enrichment analysis to identify
#              significant associations between diagnosed genes and HPO terms
#              or top-level phenotypic categories, focusing on intellectual
#              disability (ID) vs. other disease phenotypes.
# Sections:
#   1. Heatmap Visualization of Gene ↔ HPO Term Frequencies
#   2. Fisher's Exact Test for HPO Terms per Gene
#   3. Fisher's Test for Top-Level HPO Ancestor Terms
#   4. Statistical Comparison: HPO Terms in ID vs. Non-ID Diagnoses
################################################################################

##############################
# SECTION 1: Heatmap - Diagnosed Gene vs HPO Terms
##############################

# Prepare HPO term percentage by gene
genopheno_all_separated_hm <- genopheno %>%
  filter(`Diagnosed Gene` %in% id_genes$`Diagnosed Gene`) %>%
  separate_rows(`HPO terms(s)`, sep = "\\|") %>%
  group_by(`Diagnosed Gene`, `HPO terms(s)`) %>%
  tally() %>%
  arrange(-n)

hpo_percentage_all_hm <- genopheno_all_separated_hm %>%
  left_join(df_disease, by = "Diagnosed Gene") %>%
  mutate(Percentage = (n.x / n.y) * 100) %>%
  select(1, 2, 3, 4, 6) %>%
  filter(Percentage > 15) %>%
  filter(`HPO terms(s)` %in% most_occuring_hpoterms$`HPO terms(s)`)

# Load libraries
library(ggplot2)
library(reshape2)

# Format for heatmap plotting
heatmap_data <- dcast(hpo_percentage_all_hm, `Diagnosed Gene` ~ `HPO terms(s)`, value.var = "Percentage")
heatmap_melted <- melt(heatmap_data, id.vars = "Diagnosed Gene")

# Plot heatmap
p <- ggplot(heatmap_melted, aes(y = variable, x = `Diagnosed Gene`, fill = value)) +
  geom_tile(color = "white") +
  scale_fill_gradient(low = "blue", high = "red", name = "Percentage") +
  labs(title = "Heatmap of HPO Terms per Diagnosed Gene",
       y = "HPO Term",
       x = "Diagnosed Gene") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

print(p)

##############################
# SECTION 2: Fisher's Exact Test - Diagnosed Gene vs. HPO Terms
##############################

# Prepare tallied data
genopheno_all_separated <- genopheno %>%
  filter(`Diagnosed Gene` %in% id_genes$`Diagnosed Gene`) %>%
  separate_rows(`HPO terms(s)`, sep = "\\|") %>%
  group_by(`Diagnosed Gene`, `HPO terms(s)`) %>%
  tally() %>%
  arrange(-n)

genes <- unique(genopheno_all_separated$`Diagnosed Gene`)
hpo_terms <- unique(genopheno_all_separated$`HPO terms(s)`)

results_fisher_hpoterm <- data.frame(
  Gene = character(), 
  HPO_term = character(), 
  p_value = numeric(), 
  odds_ratio = numeric(),
  stringsAsFactors = FALSE
)

# Fisher test for each gene-term pair
for (gene in genes) {
  for (hpo in hpo_terms) {
    
    A <- sum(genopheno_all_separated$n[genopheno_all_separated$`Diagnosed Gene` == gene & genopheno_all_separated$`HPO terms(s)` == hpo])
    B <- sum(genopheno_all_separated$n[genopheno_all_separated$`Diagnosed Gene` == gene & genopheno_all_separated$`HPO terms(s)` != hpo])
    C <- sum(genopheno_all_separated$n[genopheno_all_separated$`Diagnosed Gene` != gene & genopheno_all_separated$`HPO terms(s)` == hpo])
    D <- sum(genopheno_all_separated$n[genopheno_all_separated$`Diagnosed Gene` != gene & genopheno_all_separated$`HPO terms(s)` != hpo])
    
    contingency_table <- matrix(c(A, B, C, D), nrow = 2) + 0.05
    fisher_test <- fisher.test(contingency_table)
    
    odds_ratio <- ifelse(!is.null(fisher_test$estimate), fisher_test$estimate, NA)
    confidence_interval <- fisher_test$conf.int
    
    results_fisher_hpoterm <- rbind(results_fisher_hpoterm, data.frame(
      `Diagnosed Gene` = gene, 
      `HPO terms(s)` = hpo, 
      p_value = fisher_test$p.value, 
      confidence_interval = confidence_interval,
      odds_ratio = odds_ratio
    ))
  }
}

# Adjust and filter results
print(results_fisher_hpoterm)

results_fisher_hpoterm$p_adjusted <- p.adjust(results_fisher_hpoterm$p_value, method = "BH")
results_fisher_hpoterm <- results_fisher_hpoterm %>%
  filter(p_value < 0.05) %>%
  unique()

##############################
# SECTION 3: Fisher's Test - Diagnosed Gene vs. Top-Level HPO Categories
##############################

hpo_toplevel_df <- joined_hpo %>%
  group_by(`Diagnosed Gene`, hpo_ancestors_description) %>%
  tally() %>%
  arrange(-n)

genes <- unique(hpo_toplevel_df$`Diagnosed Gene`)
hpo_terms <- unique(hpo_toplevel_df$hpo_ancestors_description)

results_ancestor_fisher <- data.frame(
  Gene = character(), 
  HPO_term = character(), 
  p_value = numeric(), 
  odds_ratio = numeric(),
  stringsAsFactors = FALSE
)

# Loop through each gene/ancestor term
for (gene in genes) {
  for (hpo in hpo_terms) {
    A <- sum(hpo_toplevel_df$n[hpo_toplevel_df$`Diagnosed Gene` == gene & hpo_toplevel_df$hpo_ancestors_description == hpo])
    B <- sum(hpo_toplevel_df$n[hpo_toplevel_df$`Diagnosed Gene` == gene & hpo_toplevel_df$hpo_ancestors_description != hpo])
    C <- sum(hpo_toplevel_df$n[hpo_toplevel_df$`Diagnosed Gene` != gene & hpo_toplevel_df$hpo_ancestors_description == hpo])
    D <- sum(hpo_toplevel_df$n[hpo_toplevel_df$`Diagnosed Gene` != gene & hpo_toplevel_df$hpo_ancestors_description != hpo])
    
    contingency_table <- matrix(c(A, B, C, D), nrow = 2)
    fisher_test <- fisher.test(contingency_table)
    
    odds_ratio <- ifelse(!is.null(fisher_test$estimate), fisher_test$estimate, NA)
    confidence_interval <- fisher_test$conf.int
    
    results_ancestor_fisher <- rbind(results_ancestor_fisher, data.frame(
      `Diagnosed Gene` = gene, 
      `hpo_ancestors_description` = hpo, 
      p_value = fisher_test$p.value, 
      confidence_interval = confidence_interval,
      odds_ratio = odds_ratio
    ))
  }
}

# Adjust p-values
print(results_ancestor_fisher)

results_ancestor_fisher$p_adjusted <- p.adjust(results_ancestor_fisher$p_value, method = "BH")

##############################
# SECTION 4: ID vs. Other Disease – Term Enrichment Comparison
##############################

# Compare HPO terms in ID vs other diseases
compare_df_id_od <- genopheno_id_separated %>%
  left_join(genopheno_no_id_separated, by = c("Diagnosed Gene", "HPO terms(s)")) %>%
  mutate(n_id = n.x, n_od = n.y) %>%
  select(1, 2, 5, 6) %>%
  replace(is.na(.), 0)

results_id_od <- data.frame(
  Diagnosed_Gene = character(),
  HPO_Term = character(),
  P_Value = numeric(),
  stringsAsFactors = FALSE
)

# Loop for Chi-Square or Fisher test based on cell counts
for (i in 1:nrow(compare_df_id_od)) {
  row <- compare_df_id_od[i, ]
  
  contingency_table <- matrix(
    c(row$n_id, row$n_od,
      sum(compare_df_id_od$n_id) - row$n_id,
      sum(compare_df_id_od$n_od) - row$n_od),
    nrow = 2, byrow = TRUE,
    dimnames = list(
      Category = c("Intellectual_Disability", "Other_Disease"),
      HPO_Term_Presence = c("Present", "Absent")
    )
  )
  
  test_result <- if (any(contingency_table < 5)) {
    fisher.test(contingency_table)
  } else {
    chisq.test(contingency_table)
  }
  
  results_id_od <- rbind(results_id_od, data.frame(
    Diagnosed_Gene = row$`Diagnosed Gene`,
    HPO_Term = row$`HPO terms(s)`,
    P_Value = test_result$p.value,
    stringsAsFactors = FALSE
  ))
}

results_id_od$Adjusted_P_Value <- p.adjust(results_id_od$P_Value, method = "BH")

# Final results
print(results_id_od)
