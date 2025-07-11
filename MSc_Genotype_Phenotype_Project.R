###########################
# GenoPheno-ID Analysis
# Author: Hussein Kassir
# Description: This script analyzes genotype-phenotype data focusing on Intellectual Disability (ID)
###########################

# Load required libraries and set library paths
.libPaths(c(.libPaths(), "/tools/aws-workspace-apps/ce/R/4.0.2"))
library(tidyverse)

###########################
# SECTION 1: Data Import
###########################

# Set working directory to the GenoPheno data location
setwd("/nas/weka.gel.zone/re_gecip/shared_allGeCIPs/dsmedley/FUSIL/GenoPheno")

# Import HPO ontology data
hpo1 = read_delim("hpo_ontology.txt")

# Import genopheno dataset
genopheno = read_delim("extract_allelic_series_30_08_24_plus_denovo.tsv") %>%
  unique()

# Harmonize disease names
genopheno <- genopheno %>%
  mutate(Disease = ifelse(Disease == "intellectual disability", "Intellectual disability", Disease))
genopheno <- genopheno %>%
  mutate(Disease = ifelse(Disease == "intellectual-disability", "Intellectual disability", Disease))

# Import phenotype annotations
phenotype_300824 = read_delim("phenotype_300824.hpoa")

# Switch to annotations directory
setwd("/nas/weka.gel.zone/re_gecip/shared_allGeCIPs/dsmedley/FUSIL/GenoPheno/annotations")

# Load gene panel and phenotype annotation files
panelapp_genes = read_delim("panelapp_by_gene_panel_aus_data_accessed_190824.txt")
gene_to_phenotype = read_delim("genes_to_phenotype_300824.txt")

###########################
# SECTION 2: Basic Exploration
###########################

dim(genopheno)
length(unique(genopheno$`GEL family ID`))

# Check for duplicate family IDs
genopheno_dups = genopheno$`GEL family ID`[duplicated(genopheno$`GEL family ID`)]
genopheno_dups_check = genopheno %>%
  filter(`GEL family ID` %in% unique(genopheno_dups))

# Count duplicate genes per family
gene_dups <- genopheno_dups_check %>%
  group_by(`Diagnosed Gene`) %>%
  tally()

# Consistency check: families with only one diagnosed gene
gene_consistency <- genopheno %>%
  group_by(`GEL family ID`) %>%
  summarise(Same_Gene = n_distinct(`Diagnosed Gene`) == 1)

# Summarize diagnosed genes per family ID
gene_occurrence <- genopheno %>%
  group_by(`GEL family ID`) %>%
  summarise(`Diagnosed Gene` = paste(unique(`Diagnosed Gene`), collapse = ", "))

# Frequency of ID diagnosis
genopheno %>%
  count(Disease == "Intellectual disability")

# Tally all diseases
all_disease <- genopheno %>%
  group_by(Disease) %>%
  tally

# Tally all genes
all_genes <- genopheno %>%
  group_by(`Diagnosed Gene`) %>%
  tally() %>%
  arrange(-n)

###########################
# SECTION 3: Focus on ID-Associated Genes
###########################

# Identify genes with ID and ≥10 cases
id_genes <- genopheno %>%
  group_by(Disease, `Diagnosed Gene`) %>%
  tally() %>%
  arrange(-n) %>%
  filter(Disease == "Intellectual disability") %>%
  filter(n >= 10)

# Filter genopheno dataset by ID genes
genopheno_genes <- genopheno %>%
  filter(`Diagnosed Gene` %in% id_genes$`Diagnosed Gene`)

# Count unique GEL family IDs for ID genes
genopheno_familyid <- genopheno_genes %>%
  group_by(`GEL family ID`) %>%
  tally() %>%
  select(1)

# Subset into ID vs non-ID cases
genopheno_id <- genopheno %>%
  filter(`Diagnosed Gene` %in% id_genes$`Diagnosed Gene`) %>%
  filter(Disease == "Intellectual disability")

genopheno_no_id <- genopheno %>%
  filter(`Diagnosed Gene` %in% id_genes$`Diagnosed Gene`) %>%
  filter(!Disease == "Intellectual disability")

# Frequency table of genes with ID
df_disease <- genopheno %>%
  group_by(`Diagnosed Gene`) %>%
  tally() %>%
  arrange(-n) %>%
  filter(`Diagnosed Gene` %in% id_genes$`Diagnosed Gene`) %>%
  mutate(percentage = (n / sum(n) * 100))

# Plot: Percentage per ID gene
ggplot(data = df_disease, aes(x = percentage, y = reorder(`Diagnosed Gene`, percentage), fill = `Diagnosed Gene`)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "Percentage of ID genes", x = "Percentage", y = "Diagnosed Gene") +
  theme(legend.position = "none")

# Plot: Count per ID gene
ggplot(data = df_disease, aes(x = n, y = reorder(`Diagnosed Gene`, percentage), fill = `Diagnosed Gene`)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "Number of diagnoses per gene (genes > 10 ID diagnoses)", x = "Number of diagnoses", y = "Diagnosed Gene") +
  theme(legend.position = "none")

###########################
# SECTION 4: Disease Comparisons by Gene
###########################

# Subset data to ID cases
df_disease_id <- genopheno %>%
  group_by(Disease, `Diagnosed Gene`) %>%
  filter(`Diagnosed Gene` %in% id_genes$`Diagnosed Gene`) %>%
  filter(Disease == "Intellectual disability") %>%
  tally() %>%
  arrange(-n) %>%
  mutate(n_id = n) %>%
  select(1, 2, 4)

# Subset data to non-ID cases
df_disease_od <- genopheno %>%
  filter(`Diagnosed Gene` %in% id_genes$`Diagnosed Gene`) %>%
  filter(Disease != "Intellectual disability") %>%
  group_by(`Diagnosed Gene`) %>%
  tally() %>%
  arrange(-n) %>%
  mutate(n_other_disease = n) %>%
  select(1, 3)

# Combine and reshape
compare_table_genes <- df_disease_id %>%
  left_join(df_disease_od, by = c("Diagnosed Gene" = "Diagnosed Gene"))

data_long_genes <- compare_table_genes %>%
  pivot_longer(cols = c(n_id, n_other_disease), names_to = "Type", values_to = "Count")

# Plot: Compare ID vs non-ID per gene
ggplot(data = data_long_genes, aes(y = `Diagnosed Gene`, x = Count, fill = Type)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_grid(~Type)

###########################
# SECTION 5: MOI and Variant-Level Summaries
###########################

# Mode of inheritance
genopheno_moi <- genopheno %>%
  group_by(MOI) %>%
  tally()

genopheno_moi <- genopheno_genes %>%
  group_by(MOI) %>%
  tally()

# MOI per gene
genopheno_genes_moi <- genopheno_genes %>%
  group_by(`Diagnosed Gene`, MOI) %>%
  tally() %>%
  arrange(-n) %>%
  inner_join(df_disease, by = c("Diagnosed Gene" = "Diagnosed Gene")) %>%
  mutate(percentage_zygosity = (n.x / n.y) * 100) %>%
  select(1, 2, 6)

# Plot: MOI per gene
ggplot(data = genopheno_genes_moi) +
  geom_bar(mapping = aes(y = `Diagnosed Gene`, x = percentage_zygosity, fill = MOI), stat = "identity", position = "dodge")

# Variant summary
variant_all <- genopheno %>%
  group_by(`Diagnosed Variant(s)`) %>%
  tally()

variant_all <- genopheno_genes %>%
  group_by(`Diagnosed Variant(s)`) %>%
  tally()

variant_consequence <- genopheno %>%
  group_by(`Consequence(s)`) %>%
  tally()

# De novo mutation summary
denovo <- genopheno %>%
  group_by(`De novo`) %>%
  tally()

denovo <- genopheno_genes %>%
  group_by(`De novo`) %>%
  tally()

###########################
# SECTION 6: Disease Breakdown per Gene
###########################

# Breakdown of diseases per gene
genopheno_genes_disease <- genopheno_genes %>%
  group_by(`Diagnosed Gene`, Disease) %>%
  tally() %>%
  arrange(-n) %>%
  inner_join(df_disease, by = c("Diagnosed Gene" = "Diagnosed Gene")) %>%
  mutate(percentage_disease = (n.x / n.y) * 100) %>%
  select(1, 2, 6)

# Exported summary data
data <- data.frame(genopheno_genes_disease)

# Function to plot disease distribution per gene
plot_diseases_by_gene <- function(gene_data) {
  ggplot(gene_data, aes(y = reorder(Disease, -percentage_disease), x = percentage_disease, fill = Disease)) +
    geom_bar(stat = "identity", width = 0.7) +
    labs(title = paste("Disease Occurrences for Gene:", gene_data$Diagnosed.Gene[1]),
         x = "Disease",
         y = "Count (n)") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "none") +
    scale_fill_viridis_d(option = "plasma")
}

# Loop to plot diseases for each gene
unique_genes <- unique(data$Diagnosed.Gene)
for (gene in unique_genes) {
  gene_data <- filter(data, Diagnosed.Gene == gene)
  print(plot_diseases_by_gene(gene_data))
}

# Genes with only ID
genes_only_id <- genopheno_genes_disease %>%
  group_by(`Diagnosed Gene`) %>%
  summarise(all_diseases = paste(unique(Disease), collapse = ", ")) %>%
  filter(all_diseases == "Intellectual disability")

###########################
# SECTION 7: HPO Term Analysis
###########################

# Extract and count HPO terms by gene
most_occuring_hpoterms <- genopheno %>%
  filter(`Diagnosed Gene` %in% id_genes$`Diagnosed Gene`) %>%
  separate_rows(`HPO terms(s)`, sep = "\\|") %>%
  group_by(`HPO terms(s)`, `Diagnosed Gene`) %>%
  tally()

# Aggregate to find most common HPO terms overall
most_occuring_hpoterms <- most_occuring_hpoterms %>%
  group_by(`HPO terms(s)`) %>%
  tally() %>%
  mutate(percentage_top_hpo = (n / 22) * 100) %>%
  filter(percentage_top_hpo > 50)

##############################
# SECTION 8: Gene ↔ HPO Term Mapping
##############################

# Breakdown of HPO term occurrences per gene (separated out by pipe-delimited terms)
genopheno_all_separated <- genopheno %>%
  filter(`Diagnosed Gene` %in% id_genes$`Diagnosed Gene`) %>%
  separate_rows(`HPO terms(s)`, sep = "\\|") %>%
  group_by(`Diagnosed Gene`, `HPO terms(s)`) %>%
  tally() %>%
  arrange(-n)

# Get total number of diagnoses per gene (used for normalization)
df_disease <- genopheno %>%
  group_by(`Diagnosed Gene`) %>%
  tally() %>%
  arrange(-n) %>%
  filter(`Diagnosed Gene` %in% id_genes$`Diagnosed Gene`) %>%
  mutate(percentage = (n / sum(n) * 100))

# Calculate % of cases with each HPO term in each gene
hpo_percentage_all <- genopheno_all_separated %>%
  left_join(df_disease, by = "Diagnosed Gene") %>%
  mutate(percentage_hpo_all = (n.x / n.y) * 100) %>%
  select(1, 2, 6)

# Filter for HPO terms that occur in >30% of diagnoses in each gene
hpo_percentage_all <- hpo_percentage_all %>%
  filter(percentage_hpo_all > 30)

# Tally overall HPO term usage by gene (another breakdown)
hpoterm_gene_all <- genopheno %>%
  filter(`Diagnosed Gene` %in% id_genes$`Diagnosed Gene`) %>%
  separate_rows(`HPO terms(s)`, sep = "\\|") %>%
  group_by(`HPO terms(s)`, `Diagnosed Gene`) %>%
  tally() %>%
  arrange(-n)

# Add percentage per gene
hpoterm_gene_all <- hpoterm_gene_all %>%
  left_join(df_disease, by = "Diagnosed Gene") %>%
  mutate(Percentage = (n.x / n.y) * 100) %>%
  select(1, 2, 6) %>%
  filter(Percentage > 20)

# Tally number of genes in which each HPO term occurs >20%
hpoterm_gene_all_2 <- hpoterm_gene_all %>%
  group_by(`HPO terms(s)`) %>%
  tally()

##############################
# SECTION 9: Visualizing HPO ↔ Gene Relationships
##############################

# Barplot: % of diagnoses per gene with each HPO term
data <- data.frame(hpo_percentage_all) %>%
  arrange(Diagnosed.Gene)

plot_hpo_by_gene <- function(gene_data) {
  ggplot(gene_data, aes(y = reorder(HPO.terms.s., percentage_hpo_all), x = percentage_hpo_all, fill = HPO.terms.s.)) +
    geom_bar(stat = "identity", width = 0.7) +
    labs(title = paste("HPO Term Occurrences for:", gene_data$Diagnosed.Gene[1]),
         x = "Frequency (%)",
         y = "HPO Term") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1),
          legend.position = "none") +
    scale_fill_viridis_d(option = "plasma")
}

# Plot for each gene
unique_genes <- unique(data$Diagnosed.Gene)
for (gene in unique_genes) {
  gene_data <- filter(data, Diagnosed.Gene == gene)
  print(plot_hpo_by_gene(gene_data))
}

# NOTES:
# - These plots show percentage of patients with each HPO term
# - Percent is relative to total number of diagnoses for that gene

##############################
# SECTION 10: Visualize HPO Terms Across Genes
##############################

data <- data.frame(hpoterm_gene_all) %>%
  arrange(Diagnosed.Gene)

plot_gene_by_hpo <- function(gene_data) {
  ggplot(gene_data, aes(y = reorder(Diagnosed.Gene, -Percentage), x = Percentage, fill = Diagnosed.Gene)) +
    geom_bar(stat = "identity", width = 0.7) +
    labs(title = paste("Gene Observed in HPO Term:", gene_data$HPO.terms.s.[1]),
         x = "Percentage of Diagnoses",
         y = "Diagnosed Gene") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1),
          legend.position = "none") +
    scale_fill_viridis_d(option = "plasma")
}

# Plot for each HPO term
unique_hpoterms <- unique(data$HPO.terms.s.)
for (hpoterm in unique_hpoterms) {
  gene_data <- filter(data, HPO.terms.s. == hpoterm)
  print(plot_gene_by_hpo(gene_data))
}

##############################
# SECTION 11: Mapping to HPO Top-Level Categories
##############################

# Separate HPO terms and join with ontology file to get top-level categories
genopheno_all_separated2 <- genopheno %>%
  filter(`Diagnosed Gene` %in% id_genes$`Diagnosed Gene`) %>%
  separate_rows(`HPO terms(s)`, sep = "\\|")

joined_hpo <- genopheno_all_separated2 %>%
  left_join(hpo1, by = c("HPO terms(s)" = "hpo_term_description")) %>%
  filter(top_level_phenotypic_abnormality != "-") %>%
  select(1, 2, 17) %>%
  unique()

# Count top-level terms per gene
hpo_toplevel_df <- joined_hpo %>%
  group_by(`Diagnosed Gene`, hpo_ancestors_description) %>%
  tally() %>%
  arrange(-n)

# Get unique HPO terms per gene (denominator for percentage)
df_disease2 <- genopheno %>%
  separate_rows(`HPO terms(s)`, sep = "\\|") %>%
  filter(`Diagnosed Gene` %in% id_genes$`Diagnosed Gene`) %>%
  left_join(hpo1, by = c("HPO terms(s)" = "hpo_term_description")) %>%
  filter(top_level_phenotypic_abnormality != "-") %>%
  select(1, 2) %>%
  unique() %>%
  group_by(`Diagnosed Gene`) %>%
  tally()

# Compute percentage of top-level terms per gene
hpo_toplevel_df <- hpo_toplevel_df %>%
  left_join(df_disease2, by = "Diagnosed Gene") %>%
  mutate(Percentage = (n.x / n.y) * 100) %>%
  filter(Percentage > 20)

##############################
# SECTION 12: Visualize HPO Top-Level Terms by Gene
##############################

data3 <- data.frame(hpo_toplevel_df)

plot_hpotoplevel_by_gene <- function(gene_data) {
  ggplot(gene_data, aes(y = reorder(hpo_ancestors_description, Percentage), x = Percentage, fill = hpo_ancestors_description)) +
    geom_bar(stat = "identity", width = 0.7) +
    labs(title = paste("Top-Level HPO Terms for:", gene_data$Diagnosed.Gene[1]),
         y = "Top-Level HPO Term",
         x = "Frequency (%)") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1),
          legend.position = "none") +
    scale_fill_viridis_d(option = "plasma")
}

# Plot for each gene
unique_genes3 <- unique(data3$Diagnosed.Gene)
for (gene in unique_genes3) {
  gene_data <- filter(data3, Diagnosed.Gene == gene)
  print(plot_hpotoplevel_by_gene(gene_data))
}

##############################
# SECTION 13: Top-Level HPO ↔ Gene Heatmap Preparation
##############################

# Aggregate top-level terms across genes
example_ancestor <- joined_hpo %>%
  group_by(hpo_ancestors_description, `Diagnosed Gene`) %>%
  tally() %>%
  left_join(df_disease2, by = "Diagnosed Gene") %>%
  mutate(Percentage = (n.x / n.y) * 100) %>%
  filter(Percentage > 20)

# Plot: Gene-level presence in top-level HPO categories
data <- data.frame(example_ancestor) %>%
  arrange(Diagnosed.Gene)

plot_gene_by_hpo <- function(gene_data) {
  ggplot(gene_data, aes(y = reorder(Diagnosed.Gene, -Percentage), x = Percentage, fill = Diagnosed.Gene)) +
    geom_bar(stat = "identity", width = 0.7) +
    labs(title = paste("Genes Observed in HPO Ancestor Term:", gene_data$hpo_ancestors_description[1]),
         x = "Percentage of Diagnoses",
         y = "Diagnosed Gene") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1),
          legend.position = "none") +
    scale_fill_viridis_d(option = "plasma")
}

# Plot for each top-level HPO category
unique_hpoterms <- unique(data$hpo_ancestors_description)
for (hpoterm in unique_hpoterms) {
  gene_data <- filter(data, hpo_ancestors_description == hpoterm)
  print(plot_gene_by_hpo(gene_data))
}

# Combine all genes per ancestor term
combined_df <- example_ancestor %>%
  group_by(hpo_ancestors_description) %>%
  summarize(`Diagnosed Gene` = paste(`Diagnosed Gene`, collapse = ", "))

# Count how many genes per ancestor term
nb_df <- example_ancestor %>%
  group_by(hpo_ancestors_description) %>%
  tally()


######################################
# SECTION 14: HPO Term Enrichment in ID vs Non-ID Cases
######################################

# -- HPO terms in ID cases only --

genopheno_id_separated <- genopheno %>%
  filter(`Diagnosed Gene` %in% id_genes$`Diagnosed Gene`) %>%
  filter(Disease == "Intellectual disability") %>%
  separate_rows(`HPO terms(s)`, sep = "\\|") %>%
  group_by(`Diagnosed Gene`, `HPO terms(s)`) %>%
  tally() %>%
  arrange(-n)

df_disease_id <- genopheno %>%
  group_by(Disease, `Diagnosed Gene`) %>%
  tally() %>%
  arrange(-n) %>%
  filter(`Diagnosed Gene` %in% id_genes$`Diagnosed Gene`) %>%
  filter(Disease == "Intellectual disability") %>%
  mutate(percentage_id = (n / sum(n) * 100))

hpo_percentage_id <- genopheno_id_separated %>%
  left_join(df_disease_id, by = "Diagnosed Gene") %>%
  mutate(Percentage_hpo_Id = (n.x / n.y) * 100) %>%
  select(1, 2, 4, 7) %>%
  filter(Percentage_hpo_Id > 10)

# export_2 (Optional export placeholder)

# -- HPO terms in non-ID cases --

genopheno_no_id_separated <- genopheno %>%
  filter(`Diagnosed Gene` %in% id_genes$`Diagnosed Gene`) %>%
  filter(Disease != "Intellectual disability") %>%
  separate_rows(`HPO terms(s)`, sep = "\\|") %>%
  group_by(`Diagnosed Gene`, `HPO terms(s)`) %>%
  tally() %>%
  arrange(-n)

df_disease_no_id <- genopheno %>%
  filter(Disease != "Intellectual disability") %>%
  group_by(`Diagnosed Gene`) %>%
  tally() %>%
  arrange(-n) %>%
  filter(`Diagnosed Gene` %in% id_genes$`Diagnosed Gene`) %>%
  mutate(percentage = (n / sum(n) * 100))

hpo_percentage_no_id <- genopheno_no_id_separated %>%
  left_join(df_disease_no_id, by = "Diagnosed Gene") %>%
  mutate(Percentage_hpo_Other_Diseases = (n.x / n.y) * 100) %>%
  select(1, 2, 6) %>%
  filter(Percentage_hpo_Other_Diseases > 10)

# export_3 (Optional export placeholder)

######################################
# SECTION 15: Compare HPO Terms in ID vs Other Diseases
######################################

# Merge and pivot data for visualization
compare_table_hpoterm <- full_join(hpo_percentage_id, hpo_percentage_no_id,
                                   by = c("Diagnosed Gene", "HPO terms(s)")) %>%
  select(1, 2, 4, 5) %>%
  arrange(`Diagnosed Gene`)

df_compare_table_hpoterm <- data.frame(compare_table_hpoterm)

data_long_genes_hpoterm <- df_compare_table_hpoterm %>%
  pivot_longer(cols = c(Percentage_hpo_Other_Diseases, Percentage_hpo_Id),
               names_to = "Type",
               values_to = "Percentage")

# Plot HPO term comparison by gene
genes <- unique(df_compare_table_hpoterm$Diagnosed.Gene)

for (gene in genes) {
  gene_data <- subset(data_long_genes_hpoterm, Diagnosed.Gene == gene)
  
  p <- ggplot(gene_data, aes(y = reorder(HPO.terms.s., -Percentage), x = Percentage, fill = Type)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.7), width = 0.7) +
    scale_fill_manual(values = c("Percentage_hpo_Other_Diseases" = "blue", "Percentage_hpo_Id" = "red")) +
    labs(title = paste("HPO Terms for", gene),
         y = "HPO Term",
         x = "Percentage",
         fill = "Percentage Type") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  print(p)
}

######################################
# SECTION 16: Top-Level HPO Terms in ID vs Non-ID
######################################

# -- Top-level HPO terms for ID cases --

genopheno_all_separated_id <- genopheno %>%
  filter(`Diagnosed Gene` %in% id_genes$`Diagnosed Gene`) %>%
  filter(Disease == "Intellectual disability") %>%
  separate_rows(`HPO terms(s)`, sep = "\\|")

joined_hpo_id <- genopheno_all_separated_id %>%
  left_join(hpo1, by = c("HPO terms(s)" = "hpo_term_description")) %>%
  filter(top_level_phenotypic_abnormality != "-") %>%
  select(1, 2, 17) %>%
  unique()

df_disease_id2 <- genopheno_all_separated_id %>%
  left_join(hpo1, by = c("HPO terms(s)" = "hpo_term_description")) %>%
  filter(top_level_phenotypic_abnormality != "-") %>%
  select(1, 2) %>%
  unique() %>%
  group_by(`Diagnosed Gene`) %>%
  tally()

hpo_toplevel_id <- joined_hpo_id %>%
  group_by(`Diagnosed Gene`, hpo_ancestors_description) %>%
  tally() %>%
  arrange(-n) %>%
  left_join(df_disease_id2, by = "Diagnosed Gene") %>%
  mutate(Percentage_hpo_ancetor_ID = (n.x / n.y) * 100) %>%
  filter(Percentage_hpo_ancetor_ID > 10)

# -- Top-level HPO terms for other diseases --

genopheno_all_separated_no_id <- genopheno %>%
  filter(`Diagnosed Gene` %in% id_genes$`Diagnosed Gene`) %>%
  filter(Disease != "Intellectual disability") %>%
  separate_rows(`HPO terms(s)`, sep = "\\|")

joined_hpo_no_id <- genopheno_all_separated_no_id %>%
  left_join(hpo1, by = c("HPO terms(s)" = "hpo_term_description")) %>%
  filter(top_level_phenotypic_abnormality != "-") %>%
  select(1, 2, 17) %>%
  unique()

df_disease_no_id2 <- genopheno_all_separated_no_id %>%
  left_join(hpo1, by = c("HPO terms(s)" = "hpo_term_description")) %>%
  filter(top_level_phenotypic_abnormality != "-") %>%
  select(1, 2) %>%
  unique() %>%
  group_by(`Diagnosed Gene`) %>%
  tally()

hpo_toplevel_no_id <- joined_hpo_no_id %>%
  group_by(`Diagnosed Gene`, hpo_ancestors_description) %>%
  tally() %>%
  left_join(df_disease_no_id2, by = "Diagnosed Gene") %>%
  mutate(Percentage_hpo_ancetor_Other_diseases = (n.x / n.y) * 100) %>%
  filter(Percentage_hpo_ancetor_Other_diseases > 10)

######################################
# SECTION 17: Visual Comparison of Top-Level Terms
######################################

# Combine ID and non-ID top-level percentages
compare_table_hpoterm_2 <- full_join(hpo_toplevel_id, hpo_toplevel_no_id,
                                     by = c("Diagnosed Gene", "hpo_ancestors_description")) %>%
  arrange(`Diagnosed Gene`)

df_compare_table_hpoterm_2 <- data.frame(compare_table_hpoterm_2)

data_long_genes_hpoterm_2 <- df_compare_table_hpoterm_2 %>%
  pivot_longer(cols = c(Percentage_hpo_ancetor_ID, Percentage_hpo_ancetor_Other_diseases),
               names_to = "Type",
               values_to = "Frequencies")

# Plot each gene's top-level term usage
genes <- unique(df_compare_table_hpoterm_2$Diagnosed.Gene)

for (gene in genes) {
  gene_data <- subset(data_long_genes_hpoterm_2, Diagnosed.Gene == gene)
  
  p <- ggplot(gene_data, aes(y = reorder(hpo_ancestors_description, -Frequencies), x = Frequencies, fill = Type)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.7), width = 0.7) +
    scale_fill_manual(values = c("Percentage_hpo_ancetor_ID" = "green4", "Percentage_hpo_ancetor_Other_diseases" = "orange")) +
    labs(title = paste("Top-Level HPO Terms for", gene),
         x = "Frequencies %",
         y = "HPO Ancestor Term",
         fill = "Frequencies Type") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  print(p)
}
