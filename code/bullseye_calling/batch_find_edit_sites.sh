#!/bin/bash -l
#
# Batch Processing Script for Finding Edit Sites
# This script processes multiple matrix files to identify m6A sites
#

# Display usage
function usage {
    echo "Usage: $0 [options]"
    echo "Options:"
    echo "  -c, --control FILE        Control matrix file (required)"
    echo "  -i, --input-dir DIR       Directory containing edited matrix files (required)"
    echo "  -o, --output-dir DIR      Output directory for results (required)"
    echo "  -a, --annotation FILE     Annotation file in refFlat format (required)"
    echo "  -p, --pattern PATTERN     Pattern to match matrix files (default: \"*.matrix\")"
    echo "  -t, --cpu INT             Number of CPU cores to use (default: 4)"
    echo "  --min-edit INT            Minimum edit percent (default: 10)"
    echo "  --max-edit INT            Maximum edit percent (default: 95)"
    echo "  --fold FLOAT              Edit fold threshold (default: 1.5)"
    echo "  --min-sites INT           Minimum edit sites (default: 2)"
    echo "  --control-cov INT         Control minimum coverage (default: 10)"
    echo "  --edited-cov INT          Edited minimum coverage (default: 10)"
    echo "  --intron                  Include introns in analysis"
    echo "  --ext-utr INT             Extended UTR size (for protein coding genes)"
    echo "  -h, --help                Display this help message"
    exit 1
}

# Parse command line arguments
CONTROL_FILE=""
INPUT_DIR=""
OUTPUT_DIR=""
ANNOTATION_FILE=""
PATTERN="*.matrix"
CPU=4
MIN_EDIT=10
MAX_EDIT=95
FOLD_THRESHOLD=1.5
MIN_SITES=2
CONTROL_COV=10
EDITED_COV=10
INTRON=""
EXT_UTR=""

while [[ $# -gt 0 ]]; do
    key="$1"
    case $key in
        -c|--control)
            CONTROL_FILE="$2"
            shift 2
            ;;
        -i|--input-dir)
            INPUT_DIR="$2"
            shift 2
            ;;
        -o|--output-dir)
            OUTPUT_DIR="$2"
            shift 2
            ;;
        -a|--annotation)
            ANNOTATION_FILE="$2"
            shift 2
            ;;
        -p|--pattern)
            PATTERN="$2"
            shift 2
            ;;
        -t|--cpu)
            CPU="$2"
            shift 2
            ;;
        --min-edit)
            MIN_EDIT="$2"
            shift 2
            ;;
        --max-edit)
            MAX_EDIT="$2"
            shift 2
            ;;
        --fold)
            FOLD_THRESHOLD="$2"
            shift 2
            ;;
        --min-sites)
            MIN_SITES="$2"
            shift 2
            ;;
        --control-cov)
            CONTROL_COV="$2"
            shift 2
            ;;
        --edited-cov)
            EDITED_COV="$2"
            shift 2
            ;;
        --intron)
            INTRON="--intron"
            shift
            ;;
        --ext-utr)
            EXT_UTR="--extUTR $2"
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
if [[ -z "$CONTROL_FILE" || -z "$INPUT_DIR" || -z "$OUTPUT_DIR" || -z "$ANNOTATION_FILE" ]]; then
    echo "Error: Missing required parameters"
    usage
fi

# Check if control file exists
if [ ! -e "$CONTROL_FILE" ]; then
    echo "Error: Control matrix file not found: $CONTROL_FILE"
    exit 1
fi

# Check if annotation file exists
if [ ! -e "$ANNOTATION_FILE" ]; then
    echo "Error: Annotation file not found: $ANNOTATION_FILE"
    exit 1
fi

# Check if input directory exists
if [ ! -d "$INPUT_DIR" ]; then
    echo "Error: Input directory does not exist: $INPUT_DIR"
    exit 1
fi

# Create output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

# Find all matrix files matching the pattern
MATRIX_FILES=( $(find "$INPUT_DIR" -name "$PATTERN") )

# Check if any matrix files were found
if [ ${#MATRIX_FILES[@]} -eq 0 ]; then
    echo "Error: No matrix files found in $INPUT_DIR matching pattern: $PATTERN"
    exit 1
fi

# Activate Bullseye using conda
module load miniconda
conda activate Bullseye

echo "=== Batch Processing Edit Site Finding ==="
echo "Control matrix: $CONTROL_FILE"
echo "Input directory: $INPUT_DIR"
echo "Output directory: $OUTPUT_DIR"
echo "Annotation file: $ANNOTATION_FILE"
echo "Number of files to process: ${#MATRIX_FILES[@]}"
echo "Settings:"
echo "  CPU cores: $CPU"
echo "  Min edit percent: $MIN_EDIT"
echo "  Max edit percent: $MAX_EDIT"
echo "  Edit fold threshold: $FOLD_THRESHOLD"
echo "  Min edit sites: $MIN_SITES"
echo "  Control min coverage: $CONTROL_COV"
echo "  Edited min coverage: $EDITED_COV"
if [ -n "$INTRON" ]; then
    echo "  Include introns: Yes"
fi
if [ -n "$EXT_UTR" ]; then
    echo "  Extended UTR size: ${EXT_UTR#--extUTR }"
fi
echo

# Process each matrix file
for matrix_file in "${MATRIX_FILES[@]}"; do
    # Extract sample name
    sample=$(basename "$matrix_file" .matrix)
    
    echo "Processing file: $sample.matrix"
    
    # Output file name
    output_file="$OUTPUT_DIR/${sample}_edit_sites.tsv"
    
    # Build Find_edit_site.pl command
    cmd="perl Find_edit_site.pl"
    cmd+=" --annotationFile \"$ANNOTATION_FILE\""
    cmd+=" --controlMatrix \"$CONTROL_FILE\""
    cmd+=" --EditedMatrix \"$matrix_file\""
    cmd+=" --outfile \"$output_file\""
    cmd+=" --minEdit $MIN_EDIT"
    cmd+=" --maxEdit $MAX_EDIT"
    cmd+=" --editFoldThreshold $FOLD_THRESHOLD"
    cmd+=" --MinEditSites $MIN_SITES"
    cmd+=" --cpu $CPU"
    cmd+=" --ControlMinCoverage $CONTROL_COV"
    cmd+=" --EditedMinCoverage $EDITED_COV"
    
    # Add optional parameters
    if [ -n "$INTRON" ]; then
        cmd+=" $INTRON"
    fi
    
    if [ -n "$EXT_UTR" ]; then
        cmd+=" $EXT_UTR"
    fi
    
    cmd+=" --verbose"
    
    echo "Running command: $cmd"
    eval $cmd
    
    # Check if the command was successful
    if [ $? -eq 0 ]; then
        echo "Edit site finding completed successfully for $sample"
        echo "Results saved to: $output_file"
    else
        echo "Error: Edit site finding failed for $sample"
    fi
    
    echo
done

echo "=== Batch processing completed ==="
echo "All result files are in: $OUTPUT_DIR" 