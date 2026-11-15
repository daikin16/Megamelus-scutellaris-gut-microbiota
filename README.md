# Supplementary code and data for "Comparison of the gut microbiota and predicted functional potential of the water hyacinth planthopper *Megamelus scutellaris* (Hemiptera: Delphacidae) between field and mass-reared populations"

[![DOI](https://zenodo.org/badge/1458047945.svg)](https://doi.org/10.5281/zenodo.24228800)

R scripts, command-line scripts and input data supporting the manuscript above.
Study populations: Waainek (mass-rearing station) and the field sites Middle Lake (also labelled Benoni Lake in the QIIME2 metadata) and Civic Dam (also labelled Civic Lake), East Rand, Johannesburg, South Africa.

## Contents

| File | Purpose |
|---|---|
| `scripts/QIIME2_pipeline.txt` | QIIME2 processing of the 16S rRNA gene reads (primer removal, DADA2, taxonomy, filtering, tree, core diversity metrics) |
| `scripts/01_Diversity.R` | QIIME2 import, alpha diversity, beta diversity (NMDS, PERMANOVA), genus composition, core microbiota |
| `scripts/02_KEGG_pathway_pred.R` | PICRUSt2 KO to KEGG pathway differential abundance (LinDA, `ggpicrust2`) |
| `scripts/03_Taxa_contribution.R` | Genus contributions to selected MetaCyc pathways |
| `data/qiime2/` | QIIME2 inputs for script 01: feature table, taxonomy, rooted tree, sample metadata |
| `data/metadata_ko.txt` | Sample metadata (sample ID and site) for the functional analyses |
| `data/ko_unstrat.tsv` | PICRUSt2 unstratified KO abundance table |
| `data/path_abun_contrib.tsv.gz` | PICRUSt2 stratified pathway contributions |
| `data/picrust2/` | PICRUSt2 nearest sequenced taxon index (NSTI) outputs: per-sample weighted NSTI (`weighted_nsti.tsv`) and per-ASV NSTI (`combined_marker_predicted_and_nsti.tsv`), summarised in Supplementary Table 10 |

## Run order

1. `QIIME2_pipeline.txt` (command-line, QIIME2)
2. `01_Diversity.R`
3. `02_KEGG_pathway_pred.R`
4. `03_Taxa_contribution.R`

Set the working directory to `data/` for scripts 02 and 03 (they read files by name).

## Data notes

Raw 16S rRNA gene sequence data are available in the NCBI Sequence Read Archive under BioProject PRJNA1454826 (https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1454826).

`QIIME2_pipeline.txt` reads a manifest file (`Manifestfrank.txt`) listing the paths to the raw FASTQ files. The manifest is not included; recreate it from the SRA data.

`01_Diversity.R` reads the QIIME2 files in `data/qiime2/` (set the working directory there). Data are rarefied to 40,000 reads per sample (`rarefy_even_depth`, seed 123); the rarefied object `datasubsampled_fixed.rds` is created by the commented lines in section 5.1 of that script.

`03_Taxa_contribution.R` builds its own `taxonomy` and `contrib_metadata` objects from `data/qiime2/taxonomy-opt250.qza` and `data/metadata_ko.txt`.

## Software

Analyses were run in R 4.4.2 (2024-10-31 ucrt) on Windows 11 x64. PICRUSt2 v2.6.2 was used for functional prediction. QIIME2 version: ADD VERSION.

## Citation

ADD MANUSCRIPT CITATION AFTER PUBLICATION. Code and data archive: https://doi.org/10.5281/zenodo.24228800
