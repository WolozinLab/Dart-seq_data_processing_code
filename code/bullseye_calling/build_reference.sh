#!/bin/bash -l
#
# Build Reference Files for Bullseye Pipeline
# This script constructs the necessary reference files for the Bullseye analysis
#

# Display usage
function usage {
    echo "Usage: $0 [options]"
    echo "Options:"
    echo "  -g, --gtf FILE           GTF annotation file (required)"
    echo "  -o, --output DIR         Output directory (required)"
    echo "  -m, --method [ucsc|custom] Method to convert GTF to refFlat (default: custom)"
    echo "  -h, --help               Display this help message"
    exit 1
}

# Parse command line arguments
GTF_FILE=""
OUTPUT_DIR=""
METHOD="custom"

while [[ $# -gt 0 ]]; do
    key="$1"
    case $key in
        -g|--gtf)
            GTF_FILE="$2"
            shift 2
            ;;
        -o|--output)
            OUTPUT_DIR="$2"
            shift 2
            ;;
        -m|--method)
            METHOD="$2"
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
if [[ -z "$GTF_FILE" || -z "$OUTPUT_DIR" ]]; then
    echo "Error: Missing required parameters"
    usage
fi

# Check if GTF file exists
if [ ! -e "$GTF_FILE" ]; then
    echo "Error: GTF file not found: $GTF_FILE"
    exit 1
fi

# Create output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

# Load conda environment
module load miniconda
conda activate Bullseye

echo "=== Building Reference Files for Bullseye Pipeline ==="
echo "GTF file: $GTF_FILE"
echo "Output directory: $OUTPUT_DIR"
echo "Method: $METHOD"
echo

# Set output filenames
GENEPRED_FILE="$OUTPUT_DIR/$(basename "$GTF_FILE" .gtf).genepredext"
REFFLAT_FILE="$OUTPUT_DIR/$(basename "$GTF_FILE" .gtf).refFlat"

# Method 1: Using UCSC gtfToGenePred utility
if [ "$METHOD" == "ucsc" ]; then
    echo "Using UCSC gtfToGenePred utility..."
    
    # Run gtfToGenePred
    echo "Converting GTF to GenePred..."
    gtfToGenePred -genePredExt "$GTF_FILE" "$GENEPRED_FILE"
    
    # Check if the command was successful
    if [ $? -eq 0 ]; then
        echo "GTF to GenePred conversion completed successfully"
        
        # Convert GenePred to refFlat
        echo "Generating refFlat file..."
        perl -lanE 'if ($_=~ /^\#/){say $_}else{say join("\t",$F[11] ,@F[0..9])}' "$GENEPRED_FILE" > "$REFFLAT_FILE"
        
        if [ $? -eq 0 ]; then
            echo "refFlat file generated successfully: $REFFLAT_FILE"
        else
            echo "Error generating refFlat file"
            exit 1
        fi
    else
        echo "Error: GTF to GenePred conversion failed"
        exit 1
    fi

# Method 2: Using custom gtf2genepred.pl script
elif [ "$METHOD" == "custom" ]; then
    echo "Using custom gtf2genepred.pl script..."
    
    # Run gtf2genepred.pl
    echo "Converting GTF to refFlat..."
    perl gtf2genepred.pl --gtf "$GTF_FILE" --out "$REFFLAT_FILE"
    
    # Check if the command was successful
    if [ $? -eq 0 ]; then
        echo "GTF to refFlat conversion completed successfully: $REFFLAT_FILE"
    else
        echo "Error: GTF to refFlat conversion failed"
        exit 1
    fi

# For handling hgTables.txt format
elif [ "$METHOD" == "hgtables" ]; then
    echo "Processing hgTables.txt format..."
    
    # Convert hgTables.txt to refFlat
    echo "Generating refFlat file..."
    perl -lanE 'if ($_=~ /^\#/){say $_;next;}else{$name=pop(@F)}; print join("\t",$name,@F);' "$GTF_FILE" > "$REFFLAT_FILE"
    
    # Check if the command was successful
    if [ $? -eq 0 ]; then
        echo "refFlat file generated successfully: $REFFLAT_FILE"
    else
        echo "Error generating refFlat file"
        exit 1
    fi

else
    echo "Error: Unknown method $METHOD. Use 'ucsc', 'custom', or 'hgtables'."
    exit 1
fi

echo
echo "=== Reference file construction completed ===" 