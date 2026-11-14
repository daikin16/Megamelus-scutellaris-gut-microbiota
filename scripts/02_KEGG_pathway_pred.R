################################################################################
# Predicted functional (KEGG pathway) analysis of the gut microbiome of
# Megamelus scutellaris
# Supplementary R script
#
# Sections
#   1. Load libraries
#   2. Data import
#   3. Differential abundance analysis (LinDA)
#   4. Select significant pathways per site comparison
#   5. Match sample IDs to groups
#   6. Waainek vs Middle Lake
#   7. Waainek vs Civic Dam
#   8. Export results
################################################################################


################################################################################
# 1. LOAD LIBRARIES
################################################################################

library(ggplot2)
library(readr)
library(ggpicrust2)
library(tibble)
library(tidyverse)
library(ggprism)
library(patchwork)
library(dplyr)


################################################################################
# 2. DATA IMPORT
################################################################################

# 2.1 Load metadata as a tibble
# data(metadata)
metadata <- read_delim("metadata_ko.txt", delim = "\t", escape_double = FALSE, trim_ws = TRUE)

# 2.2 Load KEGG pathway abundance
# data(kegg_abundance)
kegg_abundance <- ko2kegg_abundance("ko_unstrat.tsv")


################################################################################
# 3. DIFFERENTIAL ABUNDANCE ANALYSIS (LinDA)
################################################################################

set.seed(123)

# 3.1 Pathway differential abundance analysis (DAA) using the LinDA method
# Change group to "your_group_column" if you are not using the example dataset
daa_results_df <- pathway_daa(abundance = kegg_abundance, metadata = metadata, group = "Site",
                              daa_method = "LinDA", select = NULL, p_adjust_method = "BH", reference = "Waainek")

# 3.2 Annotate pathway results using KO to KEGG conversion
daa_annotated_sub_method_results_df <- pathway_annotation(pathway = "KO",
                                                          daa_results_df = daa_results_df, ko_to_kegg = TRUE)


################################################################################
# 4. SELECT SIGNIFICANT PATHWAYS PER SITE COMPARISON
################################################################################

# 4.1 Waainek vs Middle Lake (daa_one)
daa_one <- daa_annotated_sub_method_results_df %>%
  filter(group1 == "Waainek", group2 == "Middle Lake")

top_pathways <- daa_one %>%
  filter(p_adjust < 0.05) %>%
  arrange(p_adjust) %>%
  slice(1:30) %>%
  pull(feature)

# 4.2 Waainek vs Civic Dam (daa_two)
daa_two <- daa_annotated_sub_method_results_df %>%
  filter(group1 == "Waainek", group2 == "Civic Dam")

top_pathways_daa_two <- daa_two %>%
  filter(p_adjust < 0.05) %>%
  arrange(p_adjust) %>%
  slice(1:30) %>%
  pull(feature)


################################################################################
# 5. MATCH SAMPLE IDs TO GROUPS
################################################################################

# Replace "Sample_id" with the actual sample ID column in your metadata
sample_col <- "Sample_id"

identical(as.character(metadata[[sample_col]]), colnames(kegg_abundance))
# If this returns FALSE, then Group = metadata$Site is not safe for plotting.

sample_col <- "Sample_id"  # replace with your actual sample ID column
Group <- setNames(metadata$Site, metadata[[sample_col]])


################################################################################
# 6. WAAINEK VS MIDDLE LAKE (daa_one)
################################################################################

# 6.1 Subset metadata and abundance table
sample_col <- "Sample_id"  # replace with your actual sample ID column
keep <- metadata$Site %in% c("Waainek", "Middle Lake")
metadata_2 <- metadata[keep, ]

abundance_2 <- kegg_abundance[, metadata_2[[sample_col]], drop = FALSE]
Group_2 <- setNames(metadata_2$Site, metadata_2[[sample_col]])

# 6.2 Remove selected pathways from the plot
drop_B <- c("Adherens junction", "Endocytosis", "Sphingolipid signaling pathway",
            "cAMP signaling pathway", "Ras signaling pathway",
            "Neuroactive ligand-receptor interaction", "Basal transcription factors",
            "Staphylococcus aureus infection", "Arachidonic acid metabolism",
            "Primary bile acid biosynthesis", "Secondary bile acid biosynthesis",
            "Lipoarabinomannan (LAM) biosynthesis",
            "Biosynthesis of type II polyketide products",
            "Nonribosomal peptide structures", "Biosynthesis of various alkaloids",
            "Fatty acid elongation")
ids_B <- daa_one %>% filter(pathway_name %in% drop_B) %>% pull(feature)
top_pathways <- setdiff(top_pathways, ids_B)

# 6.3 Error bar plot (custom colours)
p_2 <- pathway_errorbar(
  abundance = abundance_2,
  daa_results_df = daa_one,
  Group = Group_2,
  p_values_threshold = 0.05,
  order = "pathway_class",
  select = top_pathways,
  ko_to_kegg = TRUE,
  colors = c("#3685bc", "#d93c3e"), # Waainek, then the other site
  legend_position = "top",
  legend_title = "Sample Groups",
  pvalue_format = "smart",
  pvalue_stars = TRUE
)

print(p_2)
ggsave("Waainek_Middle Lake.png", plot = p_2, width = 22, height = 8, dpi = 300)

# 6.4 Error bar plot (Nature-Springer colour theme)
p_nature_daa_one <- pathway_errorbar(
  abundance = abundance_2,
  daa_results_df = daa_one,
  Group = Group_2,
  p_values_threshold = 0.05,
  order = "pathway_class",
  select = top_pathways,
  ko_to_kegg = TRUE,
  color_theme = "nature",
  legend_position = "top",
  legend_title = "Sample Groups",
  pvalue_format = "smart",
  pvalue_stars = TRUE
)

print(p_nature_daa_one)

ggsave("Waainek_Middle Lake_11.png", plot = p_nature_daa_one, width = 22, height = 8, dpi = 300)


################################################################################
# 7. WAAINEK VS CIVIC DAM (daa_two)
################################################################################

# 7.1 Subset metadata and abundance table
sample_col_1 <- "Sample_id"  # replace with your actual sample ID column
keep_1 <- metadata$Site %in% c("Waainek", "Civic Dam")
metadata_3 <- metadata[keep_1, ]

abundance_3 <- kegg_abundance[, metadata_3[[sample_col_1]], drop = FALSE]
Group_3 <- setNames(metadata_3$Site, metadata_3[[sample_col]])

# 7.2 Remove selected pathways from the plot
drop_A <- c("Adherens junction", "Basal transcription factors",
            "Drug metabolism - cytochrome P450",
            "Metabolism of xenobiotics by cytochrome P450",
            "Lipoarabinomannan (LAM) biosynthesis",
            "Penicillin and cephalosporin biosynthesis",
            "Biosynthesis of type II polyketide products",
            "Tropane, piperidine and pyridine alkaloid biosynthesis",
            "Fatty acid elongation")
ids_A <- daa_two %>% filter(pathway_name %in% drop_A) %>% pull(feature)
top_pathways_daa_two <- setdiff(top_pathways_daa_two, ids_A)

# 7.3 Error bar plot (custom colours)
p_1 <- pathway_errorbar(
  abundance = abundance_3,
  daa_results_df = daa_two,
  Group = Group_3,
  p_values_threshold = 0.05,
  order = "pathway_class",
  select = top_pathways_daa_two,
  ko_to_kegg = TRUE,
  colors = c("#3685bc", "#6faa3e"),   # Waainek, then the other site
  legend_position = "top",
  legend_title = "Sample Groups",
  pvalue_format = "smart",
  pvalue_stars = TRUE
)

print(p_1)
ggsave("Waainek_Civic_DamFinal.png", plot = p_1, width = 22, height = 8, dpi = 300)

# 7.4 Error bar plot (Nature-Springer colour theme)
p_nature_daa_two <- pathway_errorbar(
  abundance = abundance_3,
  daa_results_df = daa_two,
  Group = Group_3,
  p_values_threshold = 0.05,
  order = "pathway_class",
  select = top_pathways_daa_two,
  ko_to_kegg = TRUE,
  color_theme = "nature",
  legend_position = "top",
  legend_title = "Sample Groups",
  pvalue_format = "smart",
  pvalue_stars = TRUE
)

print(p_nature_daa_two)

ggsave("Waainek_Civic_Dam_11.png", plot = p_nature_daa_two, width = 22, height = 8, dpi = 300)


################################################################################
# 8. EXPORT RESULTS
################################################################################

write.csv(daa_annotated_sub_method_results_df, file = "Waainek_Middle_Lake.csv", row.names = FALSE)
write.csv(top_pathways, file = "Top_30_Waainek_Middle_Lake.csv", row.names = FALSE)

#### END OF KEGG PATHWAY ANALYSES ####
