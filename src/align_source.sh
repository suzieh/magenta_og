#!/usr/bin/env bash

# Usage description
usage() {
	echo "Usage: $0 -i INPUT_DIR -m MAPPING_FILE -o OUTPUT_DIR -t THREADS -s ALIGN_SCORE"
	echo "Valid options for ALIGN_SCORE: 1.00, 0.99, 0.98"
	exit 1
}

# Default threads and alignment score in case not specified
ALIGN_SCORE=1.00
THREADS=4

# Parse arguments
while getopts "i:m:o:t:s:h" opt; do
    case ${opt} in
        i) INPUT_DIR="$OPTARG" ;;
	m) MAPPING_FILE="$OPTARG" ;;
        o) OUTPUT_DIR="$OPTARG" ;;
        t) THREADS="$OPTARG" ;;
	s) ALIGN_SCORE="$OPTARG" ;;
        h) usage ;;
        *) usage ;;
    esac
done

# --- CHECK FILEPATHS --- #
# Check required directories are provided, else send usage message
if [[ -z "$INPUT_DIR" || -z "$OUTPUT_DIR" ]]; then
    usage
fi
# Remove trailing slashes to avoid "//" in paths
INPUT_DIR="${INPUT_DIR%/}"
OUTPUT_DIR="${OUTPUT_DIR%/}"

# Ensure the mapping file exists and is a .txt file
if [[ ! -f "$MAPPING_FILE" || "${MAPPING_FILE##*.}" != "txt" ]]; then
    echo "Error: Mapping file must be a valid .txt file."
    exit 1
fi
if ! awk -F'\t' '{if (NF < 8) exit 1}' "$MAPPING_FILE"; then
    echo "Error: Mapping file '$MAPPING_FILE' does not contain at least 8 tab-separated columns."
    exit 1
fi

# --- CHECK ALIGNMENT IDENTITY SCORE --- #
# Set ALIGN_SCORE_MIN based on ALIGN_SCORE
case "$ALIGN_SCORE" in
    1.00) ALIGN_SCORE_MIN="C,0" ;;
    0.99) ALIGN_SCORE_MIN="L,0,-.01" ;;
    0.98) ALIGN_SCORE_MIN="L,0,-.02" ;;
    *)
        echo "Error: Invalid ALIGN_SCORE value '$ALIGN_SCORE'."
        echo "Valid options for ALIGN_SCORE: 1.00, 0.99, 0.98"
        exit 1
        ;;
esac

# --- CHECK SOFTWARE --- #
# Check if bowtie2-align-s is in the current path and usable
if ! command -v bowtie2-align-s &> /dev/null; then
    echo "Error: bowtie2-align-s not found in the current path, see installation instructions."
    exit 1
fi

# --- PREPARE OUTPUT DIR --- #
# Create output directories if not already present
mkdir -p "$OUTPUT_DIR/aligned"
mkdir -p "$OUTPUT_DIR/counts"

# --- BOWTIE_ALIGN --- #
# Use bowtie2-align-s to determine the alignment overlap between databases and samples
echo "MAGEnTa: Running Bowtie2-Align to align post samples with MAGs DBs"
THREADS=$(( THREADS < 16 ? THREADS : 16 ))
i=1
while ICF=$'\t' read -r -a myArray
do
	# skip header for now
	test $i -eq 1 && ((i=i+1)) && continue
    # store unique donor and pre filenames
    donor_id="${myArray[2]}"
    donor_file="${myArray[3]}"
    pre_id="${myArray[5]}"
    pre_file="${myArray[6]}"
    unique_donor_alignments["$donor_id|$donor_file|$pre_id"]=1
    unique_pre_alignments["$pre_id|$pre_file|$donor_id"]=1
	# donor alignment
	bowtie2-align-s -x "$OUTPUT_DIR/${myArray[2]}_db/${myArray[2]}_db" \
    -S "$OUTPUT_DIR/aligned/${myArray[1]%.fastq}_${myArray[2]}_aligned.txt" \
    --very-sensitive --np 1 --mp "1,1" --rdg "0,1" --rfg "0,1" --score-min "$ALIGN_SCORE_MIN" \
    -k 1  -p "$THREADS" --no-hd --no-unal -q "$INPUT_DIR/${myArray[1]}"
	# pre (baseline) alignment
	bowtie2-align-s -x "$OUTPUT_DIR/${myArray[5]}_db/${myArray[5]}_db" \
    -S "$OUTPUT_DIR/aligned/${myArray[1]%.fastq}_${myArray[5]}_aligned.txt" \
    --very-sensitive --np 1 --mp "1,1" --rdg "0,1" --rfg "0,1" --score-min "$ALIGN_SCORE_MIN" \
    -k 1  -p "$THREADS" --no-hd --no-unal -q "$INPUT_DIR/${myArray[1]}"
	# parse alignment file
	python3 parse_bt2_output.py ${myArray[1]} "$ALIGN_SCORE" "$OUTPUT_DIR/counts/count_table_$ALIGN_SCORE.txt"
done <  "$MAPPING_FILE"

# --- BOWTIE_ALIGN DONOR & PRE SAMPLES TO THEMSELVES / EACH OTHER--- #
# Align donor files
for set in "${!unique_donor_alignments[@]}"; do
    IFS='|' read -r donor_id donor_file pre_id <<< "$set"
    # donor self alignment
    bowtie2-align-s -x "$OUTPUT_DIR/${donor_id}_db/${donor_id}_db" \
        -S "$OUTPUT_DIR/aligned/${donor_file%.fastq}_${donor_id}_aligned.txt" \
        --very-sensitive --np 1 --mp "1,1" --rdg "0,1" --rfg "0,1" \
        --score-min "$ALIGN_SCORE_MIN" -k 1 -p "$THREADS" --no-hd --no-unal \
        -q "$INPUT_DIR/${donor_file}"
    # pre alignment
    bowtie2-align-s -x "$OUTPUT_DIR/${pre_id}_db/${pre_id}_db" \
        -S "$OUTPUT_DIR/aligned/${donor_file%.fastq}_${pre_id}_aligned.txt" \
        --very-sensitive --np 1 --mp "1,1" --rdg "0,1" --rfg "0,1" \
        --score-min "$ALIGN_SCORE_MIN" -k 1 -p "$THREADS" --no-hd --no-unal \
        -q "$INPUT_DIR/${donor_file}"
done
# Align pre-treatment files
for set in "${!unique_pre_alignments[@]}"; do
    IFS='|' read -r donor_id donor_file pre_id <<< "$set"
    # donor alignment
    bowtie2-align-s -x "$OUTPUT_DIR/${donor_id}_db/${donor_id}_db" \
        -S "$OUTPUT_DIR/aligned/${pre_file%.fastq}_${donor_id}_aligned.txt" \
        --very-sensitive --np 1 --mp "1,1" --rdg "0,1" --rfg "0,1" \
        --score-min "$ALIGN_SCORE_MIN" -k 1 -p "$THREADS" --no-hd --no-unal \
        -q "$INPUT_DIR/${pre_file}"
    # pre self alignment
    bowtie2-align-s -x "$OUTPUT_DIR/${pre_id}_db/${pre_id}_db" \
        -S "$OUTPUT_DIR/aligned/${pre_file%.fastq}_${pre_id}_aligned.txt" \
        --very-sensitive --np 1 --mp "1,1" --rdg "0,1" --rfg "0,1" \
        --score-min "$ALIGN_SCORE_MIN" -k 1 -p "$THREADS" --no-hd --no-unal \
        -q "$INPUT_DIR/${pre_file}"
done


# --- DONE --- #
echo "MAGEnTa: Completed Alignments."

