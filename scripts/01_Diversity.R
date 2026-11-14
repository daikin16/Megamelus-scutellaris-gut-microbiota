################################################################################
# Gut microbiome of Megamelus scutellaris: composition and diversity analyses
# Supplementary R script
#
# Sections
#   1. Software installation
#   2. Load libraries
#   3. Data import (QIIME2 to phyloseq)
#   4. Alpha diversity
#   5. Beta diversity (NMDS ordination and PERMANOVA)
#   6. Microbiota composition
#   7. Core microbiota
################################################################################


################################################################################
# 1. SOFTWARE INSTALLATION
################################################################################

# 1.1 qiime2R
if (!requireNamespace("devtools", quietly = TRUE)){install.packages("devtools")}
devtools::install_github("jbisanz/qiime2R")

# 1.2 BiocManager
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
BiocManager::install(version = "3.20")

# 1.3 phyloseq (Bioconductor)
BiocManager::install("phyloseq")

# 1.4 MicrobiotaProcess
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("MicrobiotaProcess")

# 1.5 pairwiseAdonis
install.packages("devtools")
if (!require("pairwiseAdonis")) devtools::install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis"); library(pairwiseAdonis)


################################################################################
# 2. LOAD LIBRARIES
################################################################################

# 2.1 Core libraries
library(qiime2R)
library(phyloseq)
library(ggplot2)
library(dplyr)

# 2.2 Additional libraries (load as needed, some affect the functioning of others)
library(lubridate)
library(ggplot2)
library(doBy)
library(reshape2)
library(vegan)
library(SPECIES)
library(phyloseq)
library(phangorn)
library(GUniFrac)

# 2.3 Further libraries (from Franco)
library("xfun", quietly = TRUE)
library("qiime2R", quietly = TRUE)
library("phyloseq", quietly = TRUE)
library("ggplot2", quietly = TRUE)
library("readxl", quietly = TRUE)
library("dplyr", quietly = TRUE)
library("tibble", quietly = TRUE)
library("utils", quietly = TRUE)
library("writexl", quietly = TRUE)
library("vegan", quietly = TRUE)
library("ggpubr", quietly = TRUE)
library("gplots", quietly = TRUE)
library("svglite", quietly = TRUE)
library("MicrobiotaProcess", quietly = TRUE)
library("plotly", quietly = TRUE)
library("pairwiseAdonis", quietly = TRUE)
library("microbial", quietly = TRUE)
library("tidyverse", quietly = TRUE)
library("microeco", quietly = TRUE)
library("GUniFrac", quietly = TRUE)
library(RColorBrewer)


################################################################################
# 3. DATA IMPORT (QIIME2 TO PHYLOSEQ)
################################################################################

# 3.1 Required QIIME2 outputs
#   filtered-table.qza       Feature (OTU/ASV) table
#   taxonomy.qza             Taxonomic assignments
#   rooted-tree.qza          Phylogenetic tree (for UniFrac)
#   sample-metadata.tsv      Metadata (plain text)

library(qiime2R)
library(phyloseq)
library(ggplot2)

# 3.2 Inspect individual components
ASVs <- read_qza("filtered-table-opt250.qza")
names(ASVs)
ASVs$data[1:5, 1:5]

metadata <- read_q2metadata("sample-metadata.tsv")
head(metadata)

taxonomy <- read_qza("taxonomy-opt250.qza")
head(taxonomy)

tree <- read_qza("rooted-tree-opt250.qza")
head(tree)

# 3.3 Import all components into one phyloseq object
physeq <- qza_to_phyloseq(
  features = "filtered-table-opt250.qza",
  taxonomy = "taxonomy-opt250.qza",
  tree = "rooted-tree-opt250.qza",
  metadata = "sample-metadata.tsv"
)

physeq


################################################################################
# 4. ALPHA DIVERSITY
################################################################################

library(ggpubr)

# 4.1 Rarefaction (subsampling to 40000 reads)
sample_sums(physeq)
datasubsampled <- rarefy_even_depth(physeq, rngseed = 123)

sample_sums(datasubsampled)

# 4.2 Alpha diversity plots (default colours)
plot_richness(datasubsampled, x = "Site", color = 'Site', measures = c("Chao1","ACE","Shannon","Simpson"), title = "Alpha diversity indices of samples") +
  geom_boxplot(alpha = 0.6) +
  theme(legend.position = "right", axis.text.x=element_text(angle=45, hjust=1, vjust=1, size=10))
ggsave("alphadiversitydatasampled.png", height = 6, width = 9, dpi = 400, device = "png")

# 4.3 Alpha diversity plots (custom colours)
plot_richness(datasubsampled, x = "Site", color = "Site",
              measures = c("Chao1", "ACE", "Shannon", "Simpson"),
              title = "Alpha diversity indices of samples") +
  geom_boxplot(alpha = 0.6) +
  scale_color_manual(values = c("Waainek" = "#F8766D",        # salmon red
                                "Benoni Lake" = "#619CFF",    # cornflower blue
                                "Civic Lake" = "#00BA38")) +  # green
  scale_fill_manual(values = c("Waainek" = "#F8766D",
                               "Benoni Lake" = "#619CFF",
                               "Civic Lake" = "#00BA38")) +
  theme(legend.position = "right",
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 10))

ggsave("alphadiversitydatasampled2.png", height = 6, width = 9, dpi = 400, device = "png")

# 4.4 Alpha diversity plots with pairwise statistics (Wilcoxon)
a_my_comparisons_group <- list(
  c("Civic Lake", "Waainek"),
  c("Benoni Lake", "Waainek"),
  c("Civic Lake", "Benoni Lake")
)

map_signif_level <- function(p) {
  ifelse(p > 0.05, "ns",
         ifelse(p <= 0.0001, "****",
                ifelse(p <= 0.001, "***",
                       ifelse(p <= 0.01, "**", "*"))))
}

plot_richness(datasubsampled, x = "Site", color = 'Site',
              measures = c("Chao1", "ACE", "Shannon", "Simpson"),
              title = "") +
  geom_boxplot(alpha = 0.6) +
  stat_compare_means(
    method = "wilcox.test",
    comparisons = a_my_comparisons_group,
    label = "p.signif",
    map_signif_level = map_signif_level,
    size = 3.5,             # Smaller label size
    fontface = "plain",     # No bold
    vjust = -0.2
  ) +
  theme(
    legend.position = "right",
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 10),  # Y-axis title
    axis.text.y = element_text(size = 10)    # Y-axis tick numbers
  )

ggsave("alphadiversitydatasubsampled_stats2.tiff", height = 6, width = 10, dpi = 320, device = "tiff")

#### END OF ALPHA DIVERSITY ANALYSES ####


################################################################################
# 5. BETA DIVERSITY
################################################################################

library(qiime2R)
library(phyloseq)
library(ggplot2)
library(dplyr)

# 5.1 Rarefied dataset (subsampled at 40000 reads)

# Rarefy once (run only the first time, then save)
#set.seed(123)
#datasubsampled <- rarefy_even_depth(physeq, rngseed = 123, verbose = FALSE)
#saveRDS(datasubsampled, "datasubsampled_fixed.rds")

# This fixed dataset was used for all analyses
datasubsampled <- readRDS("datasubsampled_fixed.rds")
#set.seed(123)


#------------------------------------------------------------------------------#
# 5.2 Bray-Curtis NMDS
#------------------------------------------------------------------------------#

Bray <- ordinate(datasubsampled, method = "NMDS", distance = "bray")

# Extract ordination points (Bray-Curtis NMDS)
ord_bray <- as.data.frame(Bray$points)
ord_bray$SampleID <- rownames(ord_bray)

# Get metadata and merge with ordination coordinates
metadata <- as(sample_data(datasubsampled), "data.frame")
metadata$SampleID <- rownames(metadata)
ord_bray <- merge(ord_bray, metadata, by = "SampleID")

# Create NMDS plot for Bray-Curtis
ordination_plot_bray <- ggplot(ord_bray, aes(x = MDS1, y = MDS2, color = Site, fill = Site)) +
  geom_point(size = 3) +
  stat_ellipse(geom = "polygon", alpha = 0.1, level = 0.95) +
  scale_color_manual(values = c("Waainek" = "#F8766D",        # salmon red
                                "Benoni Lake" = "#619CFF",    # cornflower blue
                                "Civic Lake" = "#00BA38")) +  # green
  scale_fill_manual(values = c("Waainek" = "#F8766D",
                               "Benoni Lake" = "#619CFF",
                               "Civic Lake" = "#00BA38")) +
  theme_bw() +
  labs(
    title = "",
    x = "NMDS1",
    y = "NMDS2",
    color = "Site",
    fill = "Site"
  ) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

# Save the Bray-Curtis NMDS plot
ggsave("BrayCurtis_RJW_reco.png", plot = ordination_plot_bray,
       dpi = 320, width = 8, height = 6, units = "in")


#------------------------------------------------------------------------------#
# 5.3 Weighted UniFrac NMDS
#------------------------------------------------------------------------------#

W_Unifrac <- ordinate(datasubsampled, method = "NMDS", distance = "unifrac", weighted = TRUE)

# Extract ordination points (Weighted UniFrac NMDS)
ord_W_Unifrac <- as.data.frame(W_Unifrac$points)
ord_W_Unifrac$SampleID <- rownames(ord_W_Unifrac)

# Get metadata and merge with ordination coordinates
metadata <- as(sample_data(datasubsampled), "data.frame")
metadata$SampleID <- rownames(metadata)
ord_W_Unifrac <- merge(ord_W_Unifrac, metadata, by = "SampleID")

# Create NMDS plot for Weighted UniFrac
ordination_plot_W_Unifrac <- ggplot(ord_W_Unifrac, aes(x = MDS1, y = MDS2, color = Site, fill = Site)) +
  geom_point(size = 3) +
  stat_ellipse(geom = "polygon", alpha = 0.1, level = 0.95) +
  scale_color_manual(values = c("Waainek" = "#F8766D",        # salmon red
                                "Benoni Lake" = "#619CFF",    # cornflower blue
                                "Civic Lake" = "#00BA38")) +  # green
  scale_fill_manual(values = c("Waainek" = "#F8766D",
                               "Benoni Lake" = "#619CFF",
                               "Civic Lake" = "#00BA38")) +
  theme_bw() +
  labs(
    title = "",
    x = "NMDS1",
    y = "NMDS2",
    color = "Site",
    fill = "Site"
  ) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

# Save the Weighted UniFrac NMDS plot
ggsave("ord_W_Unifrac_RJW_reco.png", plot = ordination_plot_W_Unifrac,
       dpi = 320, width = 8, height = 6, units = "in")


#------------------------------------------------------------------------------#
# 5.4 PERMANOVA, dispersion and pairwise tests
#------------------------------------------------------------------------------#

library(phyloseq)
library(vegan)

meta <- as(sample_data(datasubsampled), "data.frame")
table(meta$Site)   # should show 4 per site

# Distance matrices
dists <- list(
  BrayCurtis        = phyloseq::distance(datasubsampled, "bray"),
  WeightedUniFrac   = phyloseq::distance(datasubsampled, "wunifrac"),
  UnweightedUniFrac = phyloseq::distance(datasubsampled, "unifrac")
)

# Global PERMANOVA, dispersion test and pairwise PERMANOVA (BH-adjusted)
set.seed(123)
global_tab <- list(); disp_tab <- list(); pair_tab <- list()

for (nm in names(dists)) {
  d  <- dists[[nm]]
  dm <- as.matrix(d)
  stopifnot(all(rownames(dm) == rownames(meta)))

  # global PERMANOVA
  a <- adonis2(d ~ Site, data = meta, permutations = 999)
  global_tab[[nm]] <- data.frame(Metric = nm, R2 = round(a$R2[1], 3),
                                 F = round(a$F[1], 3), P = a$`Pr(>F)`[1])

  # dispersion test
  bd <- permutest(betadisper(d, meta$Site), permutations = 999)
  disp_tab[[nm]] <- data.frame(Metric = nm, F = round(bd$tab$F[1], 3),
                               P = bd$tab$`Pr(>F)`[1])

  # pairwise PERMANOVA
  pairs <- combn(unique(as.character(meta$Site)), 2, simplify = FALSE)
  pr <- do.call(rbind, lapply(pairs, function(p) {
    idx <- meta$Site %in% p
    x <- adonis2(as.dist(dm[idx, idx]) ~ Site, data = meta[idx, , drop = FALSE],
                 permutations = 999)
    data.frame(Metric = nm, Comparison = paste(p, collapse = " vs "),
               R2 = round(x$R2[1], 3), F = round(x$F[1], 2), P = x$`Pr(>F)`[1])
  }))
  pr$P_adj_BH <- round(p.adjust(pr$P, method = "BH"), 3)
  pair_tab[[nm]] <- pr
}

# Results tables
do.call(rbind, global_tab)
do.call(rbind, disp_tab)
do.call(rbind, pair_tab)

# Full dispersion test output per distance metric
for (nm in names(dists)) {
  cat("\n", nm, "\n")
  print(permutest(betadisper(dists[[nm]], meta$Site), permutations = 999)$tab)
}

#### END OF BETA DIVERSITY ANALYSES ####


################################################################################
# 6. MICROBIOTA COMPOSITION
################################################################################

library(phyloseq)
library(vegan)

# 6.1 Aggregate to genus level
ps_genus <- tax_glom(physeq, taxrank = "Genus")

# 6.2 Replace taxon IDs with genus names
tax_table_genus <- tax_table(ps_genus)
genus_names <- as.character(tax_table_genus[, "Genus"])

# Handle NAs or blank genus names
genus_names[is.na(genus_names) | genus_names == ""] <- "Unclassified"

# Update taxa names in the phyloseq object
taxa_names(ps_genus) <- genus_names

# 6.3 Keep the most abundant genera
top_genera <- names(sort(taxa_sums(ps_genus), decreasing = TRUE))[1:10]
ps_genus_top20 <- prune_taxa(top_genera, ps_genus)

# 6.4 Extract and normalise abundance data
abundance_matrix <- as.data.frame(otu_table(ps_genus_top20))
if (taxa_are_rows(ps_genus_top20)) {
  abundance_matrix <- t(abundance_matrix)
}

# Convert to relative abundance per sample
rel_abundance <- sweep(abundance_matrix, 1, rowSums(abundance_matrix), FUN = "/")

# 6.5 Principal component analysis
pca_res <- prcomp(rel_abundance, scale. = TRUE)

# 6.6 Compositional biplot
biplot(pca_res, main = "Compositional Biplot - Top 20 Genera")

#### END OF MICROBIOTA COMPOSITION ####


################################################################################
# 7. CORE MICROBIOTA
################################################################################

library(phyloseq)
library(dplyr)
library(ggplot2)
library(tibble)

# 7.1 Transform to relative abundance
ps.rel <- transform_sample_counts(physeq, function(x) x / sum(x))

# 7.2 Calculate prevalence and mean abundance per OTU
otu_mat <- as(otu_table(ps.rel), "matrix")

# Ensure taxa are rows
if (!taxa_are_rows(ps.rel)) {
  otu_mat <- t(otu_mat)
}

otu_stats <- tibble(
  OTU = rownames(otu_mat),
  Prevalence = rowSums(otu_mat > 0) / ncol(otu_mat),
  MeanAbundance = rowMeans(otu_mat)
)

# 7.3 Add taxonomy labels
tax_df <- as.data.frame(tax_table(ps.rel)) %>%
  rownames_to_column("OTU")

otu_stats <- otu_stats %>%
  left_join(tax_df, by = "OTU") %>%
  mutate(Label = ifelse(!is.na(Genus), Genus, OTU))

# 7.4 Define core taxa thresholds
otu_stats <- otu_stats %>%
  mutate(
    Core = ifelse(
      Prevalence >= 0.80 & MeanAbundance >= 1e-4,
      "CORE",
      "NON CORE"
    )
  )

# 7.5 Prevalence-abundance plot
p <- ggplot(otu_stats, aes(x = Prevalence, y = MeanAbundance)) +
  geom_point(aes(color = Core), size = 2, alpha = 0.8) +

  # Core thresholds
  geom_vline(xintercept = 0.8, linetype = "dashed", color = "red") +
  geom_hline(yintercept = 1e-4, linetype = "dashed", color = "red") +

  # Log scale for abundance
  scale_y_log10() +

  # Colors
  scale_color_manual(values = c("CORE" = "red", "NON CORE" = "black")) +

  labs(
    title = "OTU Prevalence vs. Mean Abundance",
    x = "Prevalence (Proportion of Samples)",
    y = "Mean Relative Abundance",
    color = ""
  ) +

  theme_bw(base_size = 14) +
  theme(
    legend.position = "right",
    panel.grid.minor = element_blank()
  )

p

# 7.6 Label core taxa
p +
  geom_text(
    data = subset(otu_stats, Core == "CORE"),
    aes(label = Label),
    size = 3,
    hjust = -0.1,
    vjust = 0.5,
    color = "red"
  ) +
  coord_cartesian(xlim = c(0, 1.05))

# 7.7 Label core taxa with ggrepel (avoids overlapping labels)
library(ggrepel)

p +
  geom_text_repel(
    data = subset(otu_stats, Core == "CORE"),
    aes(label = Label),
    color = "red",
    size = 3,
    max.overlaps = 30
  )

# 7.8 Save final plot (explicit resolution)
p_final <- p +
  geom_text_repel(
    data = subset(otu_stats, Core == "CORE"),
    aes(label = Label),
    color = "red",
    size = 3,
    max.overlaps = 30
  )

ggsave(
  filename = "prevalence_abundance_core_taxa.png",
  plot = p_final,
  width = 8,       # inches
  height = 6,      # inches
  dpi = 300,
  units = "in",
  bg = "white"
)

#### END OF CORE MICROBIOTA ANALYSES ####
