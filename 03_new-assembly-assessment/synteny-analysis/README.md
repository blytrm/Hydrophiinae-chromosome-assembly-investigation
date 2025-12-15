# Whole-Genome Synteny Visualization

Circular synteny plots comparing chromosome-level assemblies of *Hydrophis major* with V2R gene annotations.

## Overview

This repository contains scripts to generate circular synteny plots visualising whole-genome alignments between old and new chromosome-level assemblies of *Hydrophis major*. The plots highlight V2R (vomeronasal type-2 receptor) gene positions and their syntenic relationships across assemblies.

## Preprocessing

Whole-genome alignment and preprocessing is performed separately using MUMmer (nucmer) on Google Cloud Platform. See the [data processing](https://github.com/blytrm/Hydrophiinae-chromosome-assembly-investigation/blob/main/03_new-assembly-assessment/synteny-analysis/01_data-processing_genome_assembly-synteny-analysis.md) for:
- MUMmer alignment setup and execution
- Alignment filtering (identity ≥ 90%, length ≥ 1000bp)
- Reference file preparation
- Data export to Google Cloud Storage

Required input files from preprocessing:
*see* (*see* [data processing](https://github.com/blytrm/Hydrophiinae-chromosome-assembly-investigation/blob/main/03_new-assembly-assessment/synteny-analysis/01_data-processing_genome_assembly-synteny-analysis.md))
- `combined_ref.csv` - Chromosome sizes for both assemblies
- `filtered.csv` - Filtered synteny alignments
- `bed-v2r/new-final.bed` - V2R positions (new assembly)
- `bed-v2r/old-final.bed` - V2R positions (old assembly)


## Visualization

### Main Plot: All Chromosomes
*Complete synteny plot showing all chromosomes with rainbow-coloured links and orange V2R annotations*

![Full Synteny Plot v2](wg-synteny.png)
*Alternative view of full synteny visualisation*

![Full Synteny Plot with V2R](wg-synteny-v2r.png)
*Full synteny plot highlighting V2R synteny ribbons in orange*

![Full Synteny Plot with V2R v2](wg-synteny-v2r-v2.png)
*Refined V2R synteny visualisation*

### Subset Plot: Chromosomes 2 and Z

![Subset Plot](synteny_chr2_Z.png)
*Focused view of chromosomes 2 and Z showing detailed synteny relationships*

![Subset Plot with V2R](synteny-chr2_Z%20_v2r.png)
*Chromosomes 2 and Z with V2R gene annotations and synteny links*

### Version History

![Version 1](synteny_plot_v1.png)
*Initial version of the synteny plot*

## Features
- **Rainbow synteny links**: Colour-coded by chromosome, connecting homologous regions
- **V2R gene tracks**: Orange rectangular tracks showing V2R gene positions
- **V2R synteny ribbons**: Orange ribbons highlighting synteny links overlapping V2R regions
- **Dual plot outputs**: Full genome view and focused chr2/chrZ subset

## Requirements
- R (>= 3.6)
- R packages: `circlize`, `dplyr`

```r
install.packages(c("circlize", "dplyr"))
```


## Parameters

- Link filtering: identity ≥ 90%, length ≥ 1000bp
- Maximum links: 200,000 (sampled by identity if exceeded)
- V2R track height: 0.05
- Link transparency: 0.6 (regular), 0.85 (V2R)

