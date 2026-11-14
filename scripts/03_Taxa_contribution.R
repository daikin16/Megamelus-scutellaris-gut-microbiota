################################################################################
# Taxa contributions to predicted MetaCyc pathways in the gut microbiome of
# Megamelus scutellaris
# Supplementary R script
#
# Sections
#   1. Load libraries
#   2. Define pathways of interest
#   3. Import contributions and aggregate by genus
#   4. Rename pathway IDs
#   5. Plot taxa contributions
#   6. Save figure
################################################################################


################################################################################
# 1. LOAD LIBRARIES
################################################################################

library(ggpicrust2)
library(ggplot2)


################################################################################
# 2. DEFINE PATHWAYS OF INTEREST
################################################################################

# Pathways, in panel order
# (rows 1-2: amino acid / B vitamin biosynthesis, row 3: aromatic degradation + ppGpp)
labels <- c(
  "DAPLYSINESYN-PWY"                   = "L-lysine biosynthesis I",
  "COMPLETE-ARO-PWY"                   = "Aromatic amino acid biosynthesis",
  "ARGSYN-PWY"                         = "L-arginine biosynthesis I",
  "PWY-702"                            = "L-methionine biosynthesis II",
  "RIBOSYN2-PWY"                       = "Flavin (riboflavin) biosynthesis I",
  "PYRIDOXSYN-PWY"                     = "Pyridoxal 5'-phosphate biosynthesis I",
  "PROTOCATECHUATE-ORTHO-CLEAVAGE-PWY" = "Protocatechuate degradation II",
  "CATECHOL-ORTHO-CLEAVAGE-PWY"        = "Catechol degradation to beta-ketoadipate",
  "PPGPPMET-PWY"                       = "ppGpp biosynthesis"
)


################################################################################
# 3. IMPORT CONTRIBUTIONS AND AGGREGATE BY GENUS
################################################################################

# Rebuild from the raw contribution file with all 9 IDs
path_data <- read_pathway_contrib_file("path_abun_contrib.tsv.gz")

path_taxa_contrib <- aggregate_taxa_contributions(
  path_data, taxonomy = taxonomy, tax_level = "Genus", top_n = 10,
  pathway_ids = names(labels)
)
stopifnot(length(unique(path_taxa_contrib$function_id)) == 9)  # stops if any ID is missing


################################################################################
# 4. RENAME PATHWAY IDs
################################################################################

# Rename IDs to readable names, kept as text
path_taxa_contrib$function_id <- as.character(labels[as.character(path_taxa_contrib$function_id)])


################################################################################
# 5. PLOT TAXA CONTRIBUTIONS
################################################################################

# 5.1 Plot all 9 pathways (function_ids overrides the default of 6)
p <- taxa_contribution_bar(
  contrib_agg     = path_taxa_contrib,
  metadata        = contrib_metadata,   # needs sample_name and Site columns
  group           = "Site",
  function_ids    = unname(labels),
  n_functions     = 9,
  facet_by        = "function",
  legend_position = "bottom"
)

# 5.2 Force 3 x 3 layout in the order of 'labels'
p$data$function_id <- factor(p$data$function_id, levels = unname(labels))

p


################################################################################
# 6. SAVE FIGURE
################################################################################

ggsave("taxa_contribution_9pathways.png", p,
       width = 10, height = 14, units = "in", dpi = 300, bg = "white")

#### END OF TAXA CONTRIBUTION ANALYSES ####
