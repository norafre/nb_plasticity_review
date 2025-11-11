#!/bin/bash

# Input FASTQ files
R1_FASTQ="./data/bAllos4_tums_S1_dsRedRecCas_R1.fastq"
R2_FASTQ="./data/bAllos4_tums_S1_dsRedRecCas_R2.fastq"

# Reference FASTA and name
REFERENCE_FA="dsRedExpOpt_linRecorder.fa"
REFERENCE_NAME="dsRedExpOpt_linRecorder"

# Output prefix
OUTPUT_PREFIX="bAllos4_tums_S1_dsRedRecCas"

# Barcode threshold
BC_READ_THRESHOLD=100

# Sequence extraction parameters
KMER_LEN=140
SEQ_LEN=140

# Primer sequence
PRIMSEQ="CCACCACCTGTTCCTG"

# LINNAEUS script directory
SCRIPT_DIR="../LINNAEUS_resources/scar_extraction/"

