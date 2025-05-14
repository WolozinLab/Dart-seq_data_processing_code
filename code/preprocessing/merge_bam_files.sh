#!/bin/bash
#
# Merge BAM Files by Group
# This script merges multiple BAM files into a single file based on a list
#

# Display usage
function usage {
    echo "Usage: $0 [options]"
    echo "Options:"
    echo "  -l, --list FILE          File containing a list of BAM files to merge (required)"
    echo "  -o, --output FILE        Output merged BAM file path (required)"
    echo "  -i, --index              Create index for the merged BAM file"
    echo "  -h, --help               Display this help message"
    exit 1
}

# Parse command line arguments
FILE_LIST=""
OUTPUT_BAM=""
CREATE_INDEX=0

while [[ $# -gt 0 ]]; do
    key="$1"
    case $key in
        -l|--list)
            FILE_LIST="$2"
            shift 2
            ;;
        -o|--output)
            OUTPUT_BAM="$2"
            shift 2
            ;;
        -i|--index)
            CREATE_INDEX=1
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
if [[ -z "$FILE_LIST" || -z "$OUTPUT_BAM" ]]; then
    echo "Error: Missing required parameters"
    usage
fi

# Check if the file list exists
if [ ! -e "$FILE_LIST" ]; then
    echo "Error: File list not found: $FILE_LIST"
    exit 1
fi

# Create output directory if it doesn't exist
OUTPUT_DIR=$(dirname "$OUTPUT_BAM")
mkdir -p "$OUTPUT_DIR"

echo "=== Merging BAM Files ==="
echo "File list: $FILE_LIST"
echo "Output BAM: $OUTPUT_BAM"

# Count the number of files to merge
FILE_COUNT=$(wc -l < "$FILE_LIST")
echo "Number of files to merge: $FILE_COUNT"

# Use samtools to merge the BAM files
echo "Starting merge operation..."
samtools merge -b "$FILE_LIST" "$OUTPUT_BAM"

# Check if the merge was successful
if [ $? -eq 0 ]; then
    echo "BAM files merged successfully into: $OUTPUT_BAM"
    
    # Create index if requested
    if [ $CREATE_INDEX -eq 1 ]; then
        echo "Creating index for merged BAM file..."
        samtools index "$OUTPUT_BAM"
        
        if [ $? -eq 0 ]; then
            echo "Index created successfully: ${OUTPUT_BAM}.bai"
        else
            echo "Error creating index for merged BAM file."
            exit 1
        fi
    fi
else
    echo "Error during merging of BAM files."
    exit 1
fi

echo "=== Merge operation completed ===" 