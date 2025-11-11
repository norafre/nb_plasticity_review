#!/bin/bash

# The code here was adapted from code written by Bastiaan Spanjaard.

# config file (first arg or CONFIG env)
if [[ ${1-} == "--config" ]]; then CONFIG="$2"; shift 2; fi
CONFIG=${CONFIG:-"./config_endo.sh"}
if [[ -f "$CONFIG" ]]; then source "$CONFIG"; fi

# Define other variables for file names using the prefix
INTERMEDIATE_SAM="${PREFIX}_genes.sam"
WHITELIST_FILE="${PREFIX}_whitelist.txt"
RESTR_WL_FILE="${PREFIX}_restr_whitelist.txt"
RESTR_WL_FILE_P="${PREFIX}_restr_whitelist_paste.txt"
UMIBC_FASTQ="${PREFIX}_scar_UMIBC.fastq"
MT_INT_FASTQ="${PREFIX}_scar_mt_int.fastq"
MT_FASTQ="${PREFIX}_scar_mt.fastq"
SAM_OUTPUT="${PREFIX}_scar_mt.sam"
SCARS_ALL_TXT="${PREFIX}_scars_all.txt"
SCARS_WTBC_TXT="${PREFIX}_scars_wtbc.txt"
SCARS_GENES_WTBC_TXT="${PREFIX}_scars_genes_wtbc.txt"
SCARS_GENES_WTBC_TXT_P="${PREFIX}_scars_genes_wtbc_paste.txt"
SCARS_UMIW_TXT="${PREFIX}_scars_UMIw.txt"
SCARS_COUNTS_TXT="${PREFIX}_scars_genes_counts_UMIw.txt"
READ_COUNTS_TXT="${PREFIX}_scars_UMIw_read_counts.txt"

# Index combined reference
bwa index "$REF_FA"

# Check if GENE_FILE exists
if [[ ! -f "$GENE_FILE" ]]; then
    echo "Gene list file '$GENE_FILE' not found!"
    exit 1
fi

# Check if the reference fasta file exists
if [[ ! -f "$REF_FA" ]]; then
    echo "Reference genome file '$REF_FA' not found!"
    exit 1
fi

# Check if the barcode whitelist file exists
if [[ ! -f "$BC_WL_FILE" ]]; then
    echo "Barcode whitelist file '$BC_WL_FILE' not found!"
    exit 1
fi

# Extract all reads mapped to genes in the list
echo "Extracting reads for the following genes:"
cat "$GENE_FILE" | while read gene; do echo "$gene"; done

IFS=$'\n'
if [[ ! -f "$INTERMEDIATE_SAM" ]]; then
    samtools view "$BAM_FILE" | grep -F -f "$GENE_FILE" > "$INTERMEDIATE_SAM"
else
    echo "Intermediate SAM file exists. Skipping extraction step."
fi

# For those reads, extract gene name, barcode, UMI, and write to file
if [[ ! -f "$WHITELIST_FILE" ]]; then
    awk '{
        for(i=1; i<=NF; i++) {
            if($i ~ /GN:Z:/){ gene=$i }
            if($i ~ /CR:Z:/){ bc=$i }
            if($i ~ /CB:Z:/){ cbc=$i }
            if($i ~ /UR:Z:/){ umi=$i }
            if($i ~ /UB:Z:/){ cumi=$i }
        }
        if(gene && bc && cbc && umi && cumi){
            printf("%s\t%s\t%s\t%s\t%s\n", substr(gene, 6), substr(bc, 6), substr(cbc, 6), substr(umi, 6), substr(cumi, 6))
        }
    }' "$INTERMEDIATE_SAM" | sort | uniq >> "$WHITELIST_FILE"
else
    echo "Whitelist file exists. Skipping AWK step."
fi

# Keep only valid barcodes in the whitelist
if [[ ! -f "$RESTR_WL_FILE" ]]; then
    join -t $'\t' -1 1 -2 2 <(sort -k 1 "$BC_WL_FILE") <(sort -k 2 "$WHITELIST_FILE") -o 2.1,2.2,2.3,2.4,2.5 > "$RESTR_WL_FILE"
else
    echo "Restricted whitelist file exists. Skipping join step."
fi

# Create mappable fastq with UMIs and barcodes in sequence headers
if [[ ! -f "$UMIBC_FASTQ" ]]; then
    awk '(NR % 4 == 1) {print $1, $2} (NR % 4 == 3) {print $1} (NR % 2 == 0) {printf("%s%s\n", substr($1, 17, 12), substr($1, 1, 16))}' "$R1_FASTQ" > "$UMIBC_FASTQ"
else
    echo "UMI and BC fastq already exists. Skipping this step."
fi

if [[ ! -f "$MT_INT_FASTQ" ]]; then
    paste <(awk "(NR%2==0){printf(\"%s\n\", substr(\$0, 0, ${SUBSTR_LENGTH}))} (NR%2==1){print \$0}" "$R2_FASTQ") <(awk 'NR%4==2{printf("BC:Z:%s\n\n\n\n", $0)}' "$UMIBC_FASTQ") > "$MT_INT_FASTQ"
else
    echo "Intermediate MT fastq exists. Skipping this step."
fi


if [[ ! -f "$MT_FASTQ" ]]; then
    awk '{if (NR%4==1){print $1, $3} else{print $0}}' "$MT_INT_FASTQ" > "$MT_FASTQ"
else
    echo "Final MT fastq exists. Skipping this step."
fi


# Map!
if [[ ! -f "$SAM_OUTPUT" ]]; then
    bwa mem -t 20 -v 0 -C "$REF_FA" "$MT_FASTQ" > "$SAM_OUTPUT"
else
    echo "SAM output file already exists. Skipping mapping step."
fi

# Extract all mapped scars
if [[ ! -f "$SCARS_ALL_TXT" ]]; then
    grep -v '^@' "$SAM_OUTPUT" | awk '$3!="*"{
        for(i=1; i<=NF; i++) {
            if($i ~ /BC:Z:/){
                umibc=$i
            }
        }
        if(umibc){printf("%s\t%s\t%s\t%s\t%s\n", substr(umibc, 18, 16), substr(umibc, 6, 12), $3, $6, $10)}}' > "$SCARS_ALL_TXT"
else
    echo "Mapped scars already extracted. Skipping this step."
fi



# Cross-reference with allowed barcodes
if [[ ! -f "$SCARS_WTBC_TXT" ]]; then
    join -t $'\t' -1 1 -2 1 <(sort -k 1 "$BC_WL_FILE") <(sort -k 1 "$SCARS_ALL_TXT") -o 2.1,2.2,2.3,2.4,2.5 > "$SCARS_GENES_WTBC_TXT"
else
    echo "Cross-referenced scars with allowed barcodes already exist. Skipping this step."
fi



# Count reads for all transcribed scars
if [[ ! -f "$SCARS_UMIW_TXT" ]]; then
      # If you don't want to filter to only keep UMIs that are also in the whole transcriptome library
#     awk '{printf("%s_%s_%s\n", $1, $2, $4)}' "$RESTR_WL_FILE" > "$RESTR_WL_FILE_P"
#     awk '{printf("%s\t%s_%s_%s\n", $0, $3, $1, $2)}' "$SCARS_GENES_WTBC_TXT" > "$SCARS_GENES_WTBC_TXT_P"
#     join -t $'\t' -1 1 -2 6 <(sort "$RESTR_WL_FILE_P" | uniq) <(sort -k 6 "$SCARS_GENES_WTBC_TXT_P") -o 2.1,2.2,2.3,2.4,2.5 > "$SCARS_UMIW_TXT"
     
     # If you don't want to filter to only keep UMIs that are also in the whole transcriptome library
       awk '{printf("%s\n", $2)}' "$RESTR_WL_FILE" > "$RESTR_WL_FILE_P"
       join -t $'\t' -1 1 -2 1 <(sort "$RESTR_WL_FILE_P" | uniq) <(sort -k 1 "$SCARS_GENES_WTBC_TXT") -o 2.1,2.2,2.3,2.4,2.5 > "$SCARS_UMIW_TXT"
else
    echo "UMI-wise scar file exists. Skipping this step."
fi

# Generate read count summary
if [[ ! -f "$SCARS_COUNTS_TXT" ]]; then
    cut -f1,3,5 "$SCARS_UMIW_TXT" | sort | uniq | cut -f1,2 | uniq -c > "$SCARS_COUNTS_TXT"
else
    echo "Scar gene counts already exist. Skipping this step."
fi

if [[ ! -f "$READ_COUNTS_TXT" ]]; then
    sort "$SCARS_UMIW_TXT" | uniq -c > "$READ_COUNTS_TXT"
else
    echo "Read counts already exist. Skipping this step."
fi

# Remove all empty files at the end of the script
find . -type f -empty -delete
