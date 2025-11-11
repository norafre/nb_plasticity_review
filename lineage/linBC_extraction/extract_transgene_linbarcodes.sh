#!/bin/bash
# This code is adapted from Bastiaan Spanjaard's code in https://github.com/Bastiaanspanjaard/LINNAEUS.

# config file (first arg or CONFIG env)
if [[ ${1-} == "--config" ]]; then CONFIG="$2"; shift 2; fi
CONFIG=${CONFIG:-"./config_transgene.sh"}
if [[ -f "$CONFIG" ]]; then source "$CONFIG"; fi

# index reference fasta
bwa index "$REFERENCE_FA"

# === STEP 1: Count reads per barcode ===
BC_FREQS="${OUTPUT_PREFIX}_bc_freqs.txt"
if [[ ! -f "$BC_FREQS" ]]; then
  echo "Generating barcode frequency file..."
  awk 'NR % 4 == 2 {printf("%s\n", substr($1, 1, 16))}' "$R1_FASTQ" | sort | uniq -c > "$BC_FREQS"
else
  echo "Skipping barcode frequency (already exists: $BC_FREQS)"
fi

# === STEP 2: Select barcodes over threshold ===
BARCODE_CSV="${OUTPUT_PREFIX}_barcodes.csv"
if [[ ! -f "$BARCODE_CSV" ]]; then
  echo "Filtering barcodes above threshold..."
  awk -v thresh="$BC_READ_THRESHOLD" 'BEGIN{OFS = "\t"; bc = 1} $1 > thresh {print bc, $2; bc++}' "$BC_FREQS" > "$BARCODE_CSV"
else
  echo "Skipping barcode filtering (already exists: $BARCODE_CSV)"
fi

# === STEP 3: Create UMIBC file ===
UMIBC_FASTQ="${OUTPUT_PREFIX}_UMIBC.fastq"
if [[ ! -f "$UMIBC_FASTQ" ]]; then
  echo "Creating UMIBC fastq..."
  awk '(NR % 4 == 1) {print $1, $2} (NR % 4 == 3) {print $1} (NR % 2 == 0) {printf("%s%s\n", substr($1, 17, 10), substr($1, 1, 16))}' "$R1_FASTQ" > "$UMIBC_FASTQ"
else
  echo "Skipping UMIBC creation (already exists: $UMIBC_FASTQ)"
fi

# === STEP 4: Map and extract scars ===
SCARS_TXT="${OUTPUT_PREFIX}_scars.txt"
if [[ ! -f "$SCARS_TXT" ]]; then
  echo "Running scar extraction..."
  "$SCRIPT_DIR/scar_CIGAR_sc_10X_v2.pl" \
    -R1="$UMIBC_FASTQ" \
    -R2="$R2_FASTQ" \
    -op="$OUTPUT_PREFIX" \
    -t=1 \
    -r="$REFERENCE_FA" \
    -bc="$BARCODE_CSV" \
    -g="$REFERENCE_NAME" \
    -k="$KMER_LEN" \
    -mU=1 \
    -l="$SEQ_LEN" \
    -ps="$PRIMSEQ"
else
  echo "Skipping scar extraction (already exists: $SCARS_TXT)"
fi

# === STEP 5: Count and filter scars ===
SCAR_READS="${OUTPUT_PREFIX}_reads.txt"
SCAR_OVER1="${OUTPUT_PREFIX}_reads_over1.txt"

if [[ ! -f "$SCAR_READS" ]]; then
  echo "Generating scar read counts..."
  tail -n +2 "$SCARS_TXT" | sort -k2,2 -k3,3 -k1,1 -k7,7 | uniq -c | cut -f1-4,7 > "$SCAR_READS"
else
  echo "Skipping scar read counts (already exists: $SCAR_READS)"
fi

if [[ ! -f "$SCAR_OVER1" ]]; then
  echo "Filtering scars with multiple reads..."
  awk '$1>1' "$SCAR_READS" > "$SCAR_OVER1"
else
  echo "Skipping filtering of scars >1 (already exists: $SCAR_OVER1)"
fi

# === STEP 6: Final filtering using Python script ===
FILTERED_SCARS="${OUTPUT_PREFIX}_filtered_scars.csv"
if [[ ! -f "$FILTERED_SCARS" ]]; then
  echo "Running final scar filtering..."
  "$SCRIPT_DIR/scar_filter.py" -i "$SCAR_OVER1" -o "$FILTERED_SCARS"
else
  echo "Skipping final scar filtering (already exists: $FILTERED_SCARS)"
fi

# Remove all empty files at the end of the script
find . -type f -empty -delete

