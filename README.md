# Integrative Multi-Omics Analysis Identifies JUN and HSPB1 as Therapeutic Targets of Huangqi Guizhi Wuwu Decoction in Diabetic Peripheral Neuropathy

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.XXXXXXX.svg)](https://doi.org/10.5281/zenodo.XXXXXXX)

**⚠️ Manuscript Status: Under Review**

The corresponding manuscript is currently under review at *Journal of Diabetes Research*.  
This repository contains the analysis code supporting the findings. Once the paper is published, please cite the final article.  
For now, you may cite this repository as:  
`[Authors]. (2026). Code for "Integrative Multi-Omics Analysis Identifies JUN and HSPB1 as Therapeutic Targets of HGWD in DPN" (v1.0.0). Zenodo. https://doi.org/10.5281/zenodo.XXXXXXX`

## 1. Project Overview
This repository provides the complete analytical code used in the study **"Integrative Multi-Omics and Experimental Validation Identify JUN and HSPB1 as Direct Therapeutic Targets of Huangqi Guizhi Wuwu Decoction in Diabetic Peripheral Neuropathy"**.  

The study integrates:
- **Mendelian randomization (MR)** to identify causal plasma protein targets for DPN
- **Bulk RNA-seq** differential expression analysis
- **Single-cell transcriptomics** (scRNA-seq) of human tibial nerve

Key findings: HSPB1 and JUN are causal protein targets for DPN, and HGWD may exert therapeutic effects through these targets.

## 2. Repository Structure
.
├── README.md
├── LICENSE
├── .gitignore
├── code/
│ ├── 01_mr_analysis.R # Mendelian randomization (TwoSampleMR, MR-PRESSO)
│ ├── 02_bulk_deg.R # Bulk RNA-seq differential expression
│ └── 03_sc_rnaseq.R # Single-cell RNA-seq analysis (Seurat, CellChat, Monocle)
├── data/
│ └── README.md # Data sources and accession information
└── results/
├── figures/ # Key output plots
└── tables/ # Key result tables

text

## 3. System Requirements
- **R** version ≥ 4.1.0
- **R packages**: `TwoSampleMR`, `MRPRESSO`, `DESeq2`, `limma`, `clusterProfiler`, `Seurat`, `CellChat`, `Monocle3`, `Harmony`, `data.table`, `dplyr`, `ggplot2`
- **External software**: PLINK 1.9 (for LD clumping in MR)
- **Reference panel**: 1000 Genomes Phase 3 European subset

## 4. Installation Guide
```r
install.packages(c("data.table", "dplyr", "ggplot2", "BiocManager"))

# Bioconductor packages
BiocManager::install(c("DESeq2", "limma", "clusterProfiler", "org.Hs.eg.db"))

# TwoSampleMR
install.packages("TwoSampleMR", repos = c("https://mrcieu.r-universe.dev", "https://cloud.r-project.org"))

# MR-PRESSO
install.packages("MRPRESSO")

# Seurat and single-cell packages
install.packages("Seurat")
install.packages(c("CellChat", "Harmony", "Monocle3"))
PLINK can be downloaded from https://www.cog-genomics.org/plink/.
``` 
## 5. Usage
Run the scripts in the following order:

01_mr_analysis.R

Input: GWAS summary statistics (pQTL, mQTL, DPN outcome)

Output: MR estimates, heterogeneity/pleiotropy tests, scatter/forest/funnel/leave-one-out plots

Modify: File paths and column mappings at the top of the script.

02_bulk_deg.R

Input: Bulk RNA-seq count matrix

Output: DEG lists, volcano plots, heatmaps, pathway enrichment.

03_sc_rnaseq.R

Input: Single-cell expression matrix (e.g., GSE266026)

Output: Clustering, annotation, differential expression, cell communication, pseudotime.

All output files are automatically saved under the results/ directory.

## 6. Data Availability
All input data are publicly available. See data/README.md for detailed accessions and links.

Dataset	Source	Accession
DPN GWAS	GWAS Catalog	GCST90435713
pQTL	deCODE study	See reference [24]
mQTL	CLSA study	See references [26,27]
Single-cell RNA-seq (tibial nerve)	GEO	GSE266026
1000 Genomes LD reference	MRC IEU	Download
## 7. License
This project is licensed under the MIT License – see the LICENSE file.

## 8. Citation
For the code only: [Shuxian Lu]. (2026). Code for "Integrative Multi-Omics Analysis Identifies JUN and HSPB1 as Therapeutic Targets of HGWD in DPN" (v1.0.0). Zenodo. DOI: 10.5281/zenodo.XXXXXXX

For the paper: Citation will be updated upon publication.

## 9. Contact
For questions, please contact [Shuxian Lu] at [shuxianlu1216@163.com] or open an issue in this repository.
