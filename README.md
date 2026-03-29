# Pan-PTM Proteomics — *Drosophila* Aging & Longevity

Reproducible R-based visualization pipeline for a pan-PTM (post-translational modification) proteomics study of skeletal muscle aging in *Drosophila melanogaster*. PTMs were identified using the JUMPptm computational pipeline applied to TMT mass spectrometry data from long-lived fly strains and a parental control.

Published in: *npj Aging*, 2025;11:23.
https://doi.org/10.1038/s41514-025-00215-2

---

## Project Overview

This repository contains the R visualization scripts used to generate figures for a study comparing 8 post-translational modifications (PTMs) across Drosophila strains with exceptional longevity and those with normal aging trajectories.

**Biological system:**
- **B3** — parental strain (normal lifespan, age-associated muscle decline)
- **O1, O3** — lab-evolved long-lived strains (negligible functional senescence)
- **Ages:** Young (1-week-old) and Old (6-week-old)
- **Tissue:** Thoracic skeletal muscle

**PTMs analyzed (8 total):**
Acetylation, Carboxylation, Deamidation, Dihydroxylation, Mono-methylation, Oxidation, Phosphorylation, Ubiquitination

**Key finding:** PTMs that differentiate O vs. B strains are largely constitutive (present in both young and old age) rather than purely age-responsive, suggesting they represent an additional layer of phenotypic evolution underlying the longevity of O strains. The most prevalent significantly regulated PTMs are deamidation (n=319) and oxidation (n=115), occurring in evolutionarily conserved muscle contractile proteins and metabolic enzymes.

---

## Author Contributions

PTM identification and statistical analysis were performed by Suresh Poudel and Him K. Shrestha (St. Jude). **R-based graphing and visualization (this repository) were performed by Chia-Lung Chuang.**

---

## Repository Structure

```
PTM-proteomics-fly-aging/
│
├── README.md
└── PTM_fly_aging_visualization.R    # Complete visualization pipeline
```

---

## Analytical Workflow

The script processes JUMPptm output (Supplementary Table 4 from the published dataset) and generates the following visualizations:

### 1. Data Import and Reshaping
Parses the wide-format JUMPptm output into long format, extracting Time1, Time2, Gene1, Gene2, and Statistic (logFC / adj.P.Val) from encoded column names.

### 2. Volcano Plots (7 comparisons × 2 variants = 14 PDFs)

**Genotype comparisons (long-lived vs. parental):**
- O1 vs. B3 at 1 week
- O3 vs. B3 at 1 week
- O1 vs. B3 at 6 weeks
- O3 vs. B3 at 6 weeks

**Aging comparisons within strain (6 weeks vs. 1 week):**
- B3 old vs. young
- O1 old vs. young
- O3 old vs. young

Each comparison is generated in two variants: with gene + modification site labels on significant points, and without labels (clean version).

Significance threshold: adj.P.Val < 0.05 and |log2FC| > 1. Points are color-coded by PTM type.

### 3. PTM Regulation Bar Plots
Bidirectional bar plots showing the count of up- vs. down-regulated PTM sites (|logFC| > 1) per PTM type at 1-week and 6-week time points, faceted by genotype comparison.

### 4. Aging-Specific PTM Change Patterns
Grouped bar plot showing how many PTM sites changed with aging in each strain exclusively (B3-specific, O1-specific, O3-specific) vs. those shared across strains.

### 5. Global PTM Pattern Heatmap
Heatmap showing the percentage of significantly changed sites (adj.P.Val < 0.05) for each PTM type across all 7 genotype/age comparisons, with count annotations.

### 6. Global PTM Pattern Bar Plot
Same data as the heatmap displayed as bar plots faceted by PTM type, with bars colored by mean fold change.

### 7. Per-PTM Heatmaps
Log2FC heatmaps for each of the 8 PTM types showing regulated genes across all comparisons.

### 8. UpSet Plots (ComplexHeatmap)
Two sets of UpSet plots showing the intersection of significantly regulated PTM sites across the 7 comparisons, generated at two thresholds:
- **|logFC| > 1** — fold change threshold (excludes ubiquitination due to low n)
- **adj.P.Val < 0.05** — significance threshold (all 8 PTM types)

Each PTM type generates a separate PDF with bar annotations showing intersection sizes and right-side bars showing set sizes. Summary CSVs listing the Gene_Modsite identifiers per intersection column are also exported.

---

## Statistical Approach

PTM quantification and statistical testing were performed upstream (JUMPptm + limma) by Poudel et al. This script performs downstream visualization only. The input data contains pre-computed log2FC and BH-corrected adjusted P-values from limma's Empirical Bayes method.

**Significance thresholds used in visualization:**
- Volcano plots: adj.P.Val < 0.05 **and** |log2FC| > 1
- Heatmaps and bar plots: adj.P.Val < 0.05
- UpSet plots: both thresholds (separate outputs)

---

## Requirements

**R version:** ≥ 4.0

**CRAN packages:**
```r
install.packages(c(
  "tidyverse", "ggplot2", "ggrepel", "readxl",
  "broom", "pheatmap", "gridExtra", "RColorBrewer",
  "UpSetR", "scales"
))
```

**Bioconductor packages** (install once):
```r
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(c("ComplexHeatmap", "enrichplot"))
```

---

## Input Data

The raw input file (`O & B strains.xlsx`) corresponds to **Supplementary Table 4** of the published paper and is available through the ProteomeXchange Consortium (PRIDE accession: PXD014223).

The file contains pre-normalized, pre-tested JUMPptm output with columns encoding TMT comparisons in the format `<Time1><Gene1>.<Time2><Gene2>_<Statistic>` (e.g., `OnewkO1.OnewkB3_logFC`).

---

## Usage

```r
setwd("path/to/your/data")   # folder containing O & B strains.xlsx
source("PTM_fly_aging_visualization.R")
```

---

## Citation

Poudel S\*, Chuang C-L\*, Shrestha HK\*, Demontis F. Pan-PTM profiling identifies post-translational modifications associated with exceptional longevity and preservation of skeletal muscle function in *Drosophila*. *npj Aging*. 2025;11:23.
https://doi.org/10.1038/s41514-025-00215-2

\* These authors contributed equally.
