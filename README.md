# Zebrafish neuroblastoma single cell transcriptomic and lineage tracing analysis

**WORK IN PROGRESS**  
  
Scripts and notebooks used for analyses in manuscript 'Lineage origin and microenvironment shape neuroblastoma transcriptional state and plasticity' ([https://doi.org/10.1101/2025.10.13.682025](https://www.biorxiv.org/content/10.1101/2025.10.13.682025v1)).  
  
**Sequencing data** used for this project is deposited in GEO under GSE301974 (scRNA-seq data), GSE301860 (lineage tracing data from scRNA-seq), GSE301869 (spatial transcriptomics data).  

## Contents:

### Zebrafish neuroblastoma single cell transcriptomic analysis:
- Filtering, processing and cell type assignment of zebrafish NB primary tumor dataset and of full dataset containing data from primary tumors and allograft samples. Detailed information is provided [here](./transcriptome_analysis/GE_analysis_README.md).
- Analysis of neuroblastoma gene modules is described in more detail [here](./gene_modules/modules_README.md) and covers the following steps:
	- Gene module detection
	- Gene module scoring in scRNA-seq data

### CRISPR/Cas9-based high-throughput lineage tracing data analysis:
- All steps from sequencing fastq-files to clone calling and assigning grafted cells to primary tumor clones are described in more detail [here](./lineage/lineage_README.md)

### Gene expression analysis in clones and groups of cells across transplantation time points:
- Analysis of gene module expression differences between groups of cells and variance within groups of cells is described in more detail [here](./lineage_aware_expression_analysis/lineage_plus_transcriptome_README.md):
	- Differential module score expression tests between groups of cells
	- Calculation of gene module expression variance within a group of cells

### Open-ST spatial transcriptomics analysis:
- to be added...

