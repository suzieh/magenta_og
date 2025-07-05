#!/usr/bin/env bash

# Usage description
usage() {
	echo "Usage: $0 -i INPUT_DIR -m MAPPING_FILE -o OUTPUT_DIR -t THREADS"
	exit 1
}

# Default threads in case not specified
THREADS=4

# Parse arguments
while getopts "i:m:o:t:h" opt; do
    case ${opt} in
        i) INPUT_DIR="$OPTARG" ;;
	m) MAPPING_FILE="$OPTARG" ;;
        o) OUTPUT_DIR="$OPTARG" ;;
        t) THREADS="$OPTARG" ;;
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

# --- CHECK SOFTWARE --- #
## Check if metaspades.py is in the current path and usable
if ! command -v metaspades.py &> /dev/null; then
    echo "Error: metaspades.py not found in the current path, see installation instructions."
    exit 1
fi
# Check if bowtie2-build is in the current path and usable
if ! command -v bowtie2-build &> /dev/null; then
    echo "Error: bowtie2-build not found in the current path, see installation instructions."
    exit 1
fi

# --- PREPARE OUTPUT DIR --- #
# Create output directory if not already present
mkdir -p "$OUTPUT_DIR"

# ---  METASPADES --- #
# Use metaSPAdes to obtain MAGs per donor and pre sample
echo "MAGEnTa: Running metaSPAdes to get MAGs for each donor and pre sample..."
awk -F'\t' 'NR > 1 && !seen[$3, $4, $5]++ {print $3, $4, $5}' "$MAPPING_FILE" | while read -r donor_id donor_1 donor_2; do
	echo "--- Extracting ${donor_id} MAGs..."
	metaspades.py -t "$THREADS" -1 "$INPUT_DIR/${donor_1}" -2 "$INPUT_DIR/${donor_2}" -o "$OUTPUT_DIR/${donor_id}_mags"
done
awk -F'\t' 'NR > 1 && !seen[$6, $7, $8]++ {print $6, $7, $8}' "$MAPPING_FILE" | while read -r pre_id pre_1 pre_2; do
	echo "--- Extracting ${pre_id} MAGs..."
	metaspades.py -t "$THREADS" -1 "$INPUT_DIR/${pre_1}" -2 "$INPUT_DIR/${pre_2}" -o "$OUTPUT_DIR/${pre_id}_mags"
done

# --- BOWTIE_BUILD --- #
# Use Bowtie2-build to generate database format of scaffolds output from metaspades
echo "MAGEnTa: Running Bowtie2-Build to build databases from MAG output..."
THREADS=$(( THREADS < 16 ? THREADS : 16 ))
awk -F'\t' 'NR > 1 && !seen[$3]++ {print $3}' "$MAPPING_FILE" | while read -r donor_id; do
        echo "--- Building ${donor_id} MAGs DB..."
        bowtie2-build --seed 125 --threads "$THREADS" "$OUTPUT_DIR/${donor_id}_mags/scaffolds.fasta" "$OUTPUT_DIR/${donor_id}_db/${donor_id}_db"
done
awk -F'\t' 'NR > 1 && !seen[$6]++ {print $6}' "$MAPPING_FILE" | while read -r pre_id; do
        echo "--- Building ${pre_id} MAGs DB..."
        bowtie2-build --seed 125 --threads "$THREADS" "$OUTPUT_DIR/${pre_id}_mags/scaffolds.fasta" "$OUTPUT_DIR/${pre_id}_db/${pre_id}_db"
done

# --- DONE --- #
echo "MAGEnTa: Completed MAGs DB build."
