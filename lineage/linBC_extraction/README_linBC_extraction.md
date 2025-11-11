This directory provides the scripts and example input and output data for scar extraction. Note that example input data is a small dummy set, which can be used to test the pipeline. Example output data was generated using a larger input dataset and does not correspond to the output one will obtain, when running the provided scripts on the dummy input data.

# Contents
- yaml-file with conda environment specifics used for this analysis. Crucial tools needed are bwa and samtools  - env2.yml
- Commands needed to run lineage barcode extraction from fastq files. Calls on other scripts and files in this directory  - commands_lineage_barcode_extraction.sh 
- Configuration file for lineage barcode extraction from endogenous target fastq files - config_endo.sh
- Configuration file for lineage barcode extraction from transgenic target fastq files - config_transgene.sh
- Script for lineage barcode extraction from endogenous target fastq files  - extract_endo_linbarcodes.sh
- Script for lineage barcode extraction from transgenic target fastq files - extract_transgene_linbarcodes.sh
- Reference fasta containing cDNA sequences of  endogenous target genes - oct_targ.fa
- Reference fasta containing the cDNA sequence of the transgenic target site - dsRedExpOpt_linRecorder.fa
- List of endogenous target genes used - gene_list.txt
- Data: Example input data
	- Small example set of R1 and R2 fastq files containing data from sequencing of all targeted endogenous target libraries. Files from individual target genes were concatenated - bAllos4_tums_S1_endoTargets_Rx.fastq
	- Small example set of R1 and R2 fastq files containing data from sequencing of the targeted transgene library - bAllos4_tums_S1_dsRedRecCas_Rx.fastq
	- An example barcode whitelist as written out from Seurat or Scanpy. This provides all barcodes that belong to a valid cell that passed scRNA-seq data filtering  - barcode_whitelist_tums_s1.txt
- Output: Example output data that is used as input for downstream analysis
	- Small example tables output from endogenous lineage barcode extration to be used downstream - bAllos4_tums_S1_endoTargets_scars_genes_counts_UMIw.txt and bAllos4_tums_S1_endoTargets_scars_UMIw_read_counts.txt
	- Small example table output from transgene lineage barcode extraction to be used downstream - bAllos4_tums_S1_dsRedRecCas_filtered_scars.csv and bAllos4_tums_S1_dsRedRecCas_reads_over1.txt
	- Example stattistics file generated during transgene lineage barcode extraction - bAllos4_tums_S1_dsRedRecCas_stats.txt

# How to run scar extraction
- Prepare your sequencing data. Fastq-files must be unzipped. Fastq from multiple endogenous target genes should be concatenated into one fastq-file (per read).
- Edit the configuration files to match your input file names, R2 read lengths to consider for lineage barcode extraction and your desired prefix for output file names.
- Execute the commands in commands_lineage_barcode_extraction.sh. This includes activation of a conda environment with bwa and samtools, but this can be ommited, if you have another means of making bwa and samtools available.

# Output
- Intermediate files as well as those used for further analysis are kept. Only files that are empty after running the lineage barcode extraction are deleted.

