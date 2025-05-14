#!/bin/bash -l
#
# Dart-seq Data Processing Workflow
# This script orchestrates the entire Dart-seq data processing pipeline
#

# Display usage
function usage {
    echo "Usage: $0 [options]"
    echo "Options:"
    echo "  -r1, --read1 FILE         Read1 FASTQ file (required)"
    echo "  -r2, --read2 FILE         Read2 FASTQ file (required)"
    echo "  -g, --genome FILE         Reference genome file (required)"
    echo "  -o, --output DIR          Output directory (required)"
    echo "  -c, --min-coverage INT    Minimum coverage for Bullseye (default: 10)"
    echo "  -q, --min-quality INT     Minimum mapping quality (default: 20)"
    echo "  -h, --help                Display this help message"
    exit 1
}

# Parse command line arguments
READ1=""
READ2=""
GENOME=""
OUTPUT_DIR=""
MIN_COVERAGE=10
MIN_QUALITY=20

while [[ $# -gt 0 ]]; do
    key="$1"
    case $key in
        -r1|--read1)
            READ1="$2"
            shift 2
            ;;
        -r2|--read2)
            READ2="$2"
            shift 2
            ;;
        -g|--genome)
            GENOME="$2"
            shift 2
            ;;
        -o|--output)
            OUTPUT_DIR="$2"
            shift 2
            ;;
        -c|--min-coverage)
            MIN_COVERAGE="$2"
            shift 2
            ;;
        -q|--min-quality)
            MIN_QUALITY="$2"
            shift 2
            ;;
        -h|--help)
            usage
            ;;
        *)
            echo "Unknown option: $1"
            usage
            ;;
    esac
done

# Validate required parameters
if [[ -z "$READ1" || -z "$READ2" || -z "$GENOME" || -z "$OUTPUT_DIR" ]]; then
    echo "Error: Missing required parameters"
    usage
fi

# Create output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"
mkdir -p "$OUTPUT_DIR/trimmed"
mkdir -p "$OUTPUT_DIR/aligned"
mkdir -p "$OUTPUT_DIR/rm_dup"
mkdir -p "$OUTPUT_DIR/bullseye"

# Set script paths
SCRIPT_DIR="$(dirname "$0")"
TRIMMOMATIC_SCRIPT="$SCRIPT_DIR/preprocessing/Trimmomatic.sh"
STAR_SCRIPT="$SCRIPT_DIR/preprocessing/STAR_alignment.sh"
REMOVE_DUP_SCRIPT="$SCRIPT_DIR/preprocessing/remove_duplicates.sh"
BULLSEYE_SCRIPT="$SCRIPT_DIR/bullseye_calling/bullseye_caller.py"

# Get base name for output files
BASE_NAME=$(basename "$READ1" | cut -d_ -f1-4)

echo "=== Dart-seq Data Processing Workflow ==="
echo "Starting processing for: $BASE_NAME"
echo "Read1: $READ1"
echo "Read2: $READ2"
echo "Genome: $GENOME"
echo "Output directory: $OUTPUT_DIR"
echo

# Step 1: Trim reads
echo "=== Step 1: Trimming reads ==="
TRIM_OUT_PAIRED_R1="$OUTPUT_DIR/trimmed/${BASE_NAME}_R1_paired.fastq.gz"
TRIM_OUT_PAIRED_R2="$OUTPUT_DIR/trimmed/${BASE_NAME}_R2_paired.fastq.gz"
TRIM_OUT_UNPAIRED_R1="$OUTPUT_DIR/trimmed/${BASE_NAME}_R1_unpaired.fastq.gz"
TRIM_OUT_UNPAIRED_R2="$OUTPUT_DIR/trimmed/${BASE_NAME}_R2_unpaired.fastq.gz"

bash "$TRIMMOMATIC_SCRIPT" \
    "$READ1" "$READ2" \
    "$TRIM_OUT_PAIRED_R1" "$TRIM_OUT_PAIRED_R2" \
    "$TRIM_OUT_UNPAIRED_R1" "$TRIM_OUT_UNPAIRED_R2"

echo "Trimming completed"
echo

# Step 2: Align reads
echo "=== Step 2: Aligning reads ==="
bash "$STAR_SCRIPT" \
    "$TRIM_OUT_PAIRED_R1" "$TRIM_OUT_PAIRED_R2" \
    "$GENOME" "$OUTPUT_DIR/aligned"

# Find the BAM file created by STAR
BAM_FILE="$OUTPUT_DIR/aligned/${BASE_NAME}_Aligned.sortedByCoord.out.bam"

echo "Alignment completed"
echo

# Step 3: Remove duplicates
echo "=== Step 3: Removing duplicates ==="
# Temporarily modify the remove_duplicates.sh script with the correct paths
TMP_SCRIPT="$OUTPUT_DIR/tmp_remove_duplicates.sh"
cp "$REMOVE_DUP_SCRIPT" "$TMP_SCRIPT"

# Update the script with the correct input and output directories
sed -i.bak "s|input_dir=\"path/to/input/bams\"|input_dir=\"$OUTPUT_DIR/aligned\"|g" "$TMP_SCRIPT"
sed -i.bak "s|output_dir=\"path/to/output/bams/rm_dup\"|output_dir=\"$OUTPUT_DIR/rm_dup\"|g" "$TMP_SCRIPT"

# Execute the modified script
bash "$TMP_SCRIPT"

# Clean up
rm "$TMP_SCRIPT" "$TMP_SCRIPT.bak"

# Updated BAM file path
BAM_FILE="$OUTPUT_DIR/rm_dup/${BASE_NAME}_Aligned.sortedByCoord.rm_dup.bam"

echo "Duplicate removal completed"
echo

# Step 4: Build BAM index
echo "=== Step 4: Building BAM index ==="
echo "Indexing BAM file: $BAM_FILE"

# Check if samtools is available
if command -v samtools &> /dev/null; then
    samtools index "$BAM_FILE"
    echo "BAM indexing completed"
else
    echo "Warning: samtools not found in PATH. Skipping BAM indexing."
    echo "Please install samtools or add it to your PATH to enable BAM indexing."
fi
echo

# Step 5: Call Bullseye sites
echo "=== Step 5: Calling Bullseye sites ==="
python "$BULLSEYE_SCRIPT" \
    --input "$BAM_FILE" \
    --output "$OUTPUT_DIR/bullseye" \
    --genome "$GENOME" \
    --min-coverage "$MIN_COVERAGE" \
    --min-quality "$MIN_QUALITY"

echo "Bullseye calling completed"
echo

echo "=== Workflow completed successfully ==="
echo "Results are in: $OUTPUT_DIR/bullseye"
echo "Indexed BAM file: $BAM_FILE"
echo "BAM index file: ${BAM_FILE}.bai" 