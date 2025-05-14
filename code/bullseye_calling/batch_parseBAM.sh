#!/bin/bash -l
#
# Batch Process BAM Files with parseBAM.pl
# This script processes multiple BAM files with the parseBAM.pl script
#

# Display usage
function usage {
    echo "Usage: $0 [options]"
    echo "Options:"
    echo "  -i, --input-dir DIR       Directory containing BAM files (required)"
    echo "  -o, --output-dir DIR      Output directory for matrix files (required)"
    echo "  -p, --pattern PATTERN     Pattern to match BAM files (default: \"*.bam\")"
    echo "  -c, --cpu INT             Number of CPU cores to use (default: 4)"
    echo "  -m, --min-coverage INT    Minimum coverage threshold (default: 10)"
    echo "  --no-remove-duplicates    Do not remove duplicates"
    echo "  -h, --help                Display this help message"
    exit 1
}

# Parse command line arguments
INPUT_DIR=""
OUTPUT_DIR=""
PATTERN="*.bam"
CPU=4
MIN_COVERAGE=10
REMOVE_DUPLICATES="--removeDuplicates"

while [[ $# -gt 0 ]]; do
    key="$1"
    case $key in
        -i|--input-dir)
            INPUT_DIR="$2"
            shift 2
            ;;
        -o|--output-dir)
            OUTPUT_DIR="$2"
            shift 2
            ;;
        -p|--pattern)
            PATTERN="$2"
            shift 2
            ;;
        -c|--cpu)
            CPU="$2"
            shift 2
            ;;
        -m|--min-coverage)
            MIN_COVERAGE="$2"
            shift 2
            ;;
        --no-remove-duplicates)
            REMOVE_DUPLICATES=""
            shift
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
if [[ -z "$INPUT_DIR" || -z "$OUTPUT_DIR" ]]; then
    echo "Error: Missing required parameters"
    usage
fi

# Check if input directory exists
if [ ! -d "$INPUT_DIR" ]; then
    echo "Error: Input directory does not exist: $INPUT_DIR"
    exit 1
fi

# Create output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

# Find all BAM files matching the pattern
BAM_FILES=( $(find "$INPUT_DIR" -name "$PATTERN") )

# Check if any BAM files were found
if [ ${#BAM_FILES[@]} -eq 0 ]; then
    echo "Error: No BAM files found in $INPUT_DIR matching pattern: $PATTERN"
    exit 1
fi

# Load conda environment
module load miniconda
conda activate Bullseye

echo "=== Batch Processing BAM Files with parseBAM.pl ==="
echo "Input directory: $INPUT_DIR"
echo "Output directory: $OUTPUT_DIR"
echo "Number of files to process: ${#BAM_FILES[@]}"
echo "CPU cores per job: $CPU"
echo "Minimum coverage: $MIN_COVERAGE"
if [ -n "$REMOVE_DUPLICATES" ]; then
    echo "Remove duplicates: Yes"
else
    echo "Remove duplicates: No"
fi
echo

# Process each BAM file
for bam_file in "${BAM_FILES[@]}"; do
    # Extract sample name
    sample=$(basename "$bam_file" .bam)
    
    echo "Processing file: $sample.bam"
    
    # Run parseBAM.pl
    cmd="perl parseBAM.pl"
    cmd+=" --input \"$bam_file\""
    cmd+=" --output \"$OUTPUT_DIR/${sample}.matrix\""
    cmd+=" --cpu $CPU"
    cmd+=" --minCoverage $MIN_COVERAGE"
    
    if [ -n "$REMOVE_DUPLICATES" ]; then
        cmd+=" $REMOVE_DUPLICATES"
    fi
    
    cmd+=" --verbose"
    
    echo "Running command: $cmd"
    eval $cmd
    
    # Check if the command was successful
    if [ $? -eq 0 ]; then
        echo "ParseBAM completed successfully for $sample"
    else
        echo "Error: ParseBAM failed for $sample"
    fi
    
    echo
done

echo "=== Batch processing completed ==="
echo "All matrix files are in: $OUTPUT_DIR" 