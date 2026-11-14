# Supplementary code and data for "Characterisation of the gut microbiota and predicted functional potential of the water hyacinth planthopper *Megamelus scutellaris* (Hemiptera: Delphacidae) between field and mass-reared populations"

R scripts and input data supporting the manuscript above. Field sites: Waainek, Middle Lake (also labelled Benoni Lake in the QIIME2 metadata) and Civic Dam (also labelled Civic Lake), East Rand, Johannesburg, South Africa.

## Contents

| File | Purpose |
|---|---|
| `scripts/01_Diversity.R` | QIIME2 import, alpha diversity, beta diversity (NMDS, PERMANOVA), genus composition, core microbiota |
| `scripts/02_KEGG_pathway_pred.R` | PICRUSt2 KO to KEGG pathway differential abundance (LinDA, `ggpicrust2`) |
| `scripts/03_Taxa_contribution.R` | Genus contributions to selected MetaCyc pathways |
| `data/metadata_ko.txt` | Sample metadata (sample ID and site) for the functional analyses |
| `data/qiime2/` | QIIME2 inputs for script 01: feature table, taxonomy, rooted tree, sample metadata |
| `data/ko_unstrat.tsv` | PICRUSt2 unstratified KO abundance table |
| `data/path_abun_contrib.tsv.gz` | PICRUSt2 stratified pathway contributions |

## Run order

1. `01_Diversity.R`
2. `02_KEGG_pathway_pred.R`
3. `03_Taxa_contribution.R`

Set the working directory to `data/` for scripts 02 and 03 (they read files by name).

## Data notes

`01_Diversity.R` reads the QIIME2 files in `data/qiime2/` (set the working directory there).
The rarefied object `datasubsampled_fixed.rds` is created by the commented lines in section 5.1
(`rarefy_even_depth`, seed 123). 

`03_Taxa_contribution.R` uses two objects created elsewhere in the session: `taxonomy` (ASV taxonomy table) and
`contrib_metadata` (data frame with `sample_name` and `Site` columns).

## Software

Analyses were run in R 4.4.2 (2024-10-31 ucrt) on Windows 11 x64.

