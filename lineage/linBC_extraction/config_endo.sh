#!/bin/bash

# Define the common prefix for all output files
PREFIX="bAllos4_tums_S1_endoTargets"

# Define paths and filenames for external files
REF_FA="./oct_targ.fa"   # Path to the reference fasta file
BC_WL_FILE="./data/barcode_whitelist_tums_s1.txt" # Path to the barcode whitelist file from seurat or scanpy
BAM_FILE="/local/users/nfresma/remap_DR_NB/multiseq_bAllos4_PTs_s1/outs/possorted_genome_bam.bam" # Path to cellranger output BAM
GENE_FILE="./gene_list.txt"   # List of genes in an external file for easy updates
SUBSTR_LENGTH=120   # Set desired length of R2 input for mapping
R1_FASTQ="./data/bAllos4_tums_S1_endoTargets_R1.fastq"
R2_FASTQ="./data/bAllos4_tums_S1_endoTargets_R2.fastq"

