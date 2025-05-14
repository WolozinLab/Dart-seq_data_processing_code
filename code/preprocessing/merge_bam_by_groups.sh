#!/bin/bash
#
# Merge BAM Files by Multiple Groups
# This script merges BAM files within each group into separate merged files
#

# Display usage
function usage {
    echo "Usage: $0 [options]"
    echo "Options:"
    echo "  -i, --input-dir DIR      Directory containing BAM files (required)"
    echo "  -o, --output-dir DIR     Output directory for merged BAM files (required)"
    echo "  -g, --groups FILE        Tab-delimited file defining groups and their samples (required)"
    echo "                           Format: group_name    pattern_or_prefix"
    echo "  -p, --pattern PATTERN    Pattern for BAM files (default: \"*.bam\")"
    echo "  -s, --suffix SUFFIX      Suffix to add to output files (default: \"combined\")"
    echo "  -x, --index              Create index for merged BAM files"
    echo "  -h, --help               Display this help message"
    exit 1
}

# Parse command line arguments
INPUT_DIR=""
OUTPUT_DIR=""
GROUPS_FILE=""
PATTERN="*.bam"
SUFFIX="combined"
CREATE_INDEX=0

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
        -g|--groups)
            GROUPS_FILE="$2"
            shift 2
            ;;
        -p|--pattern)
            PATTERN="$2"
            shift 2
            ;;
        -s|--suffix)
            SUFFIX="$2"
            shift 2
            ;;
        -x|--index)
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
if [[ -z "$INPUT_DIR" || -z "$OUTPUT_DIR" || -z "$GROUPS_FILE" ]]; then
    echo "Error: Missing required parameters"
    usage
fi

# Check if input directory exists
if [ ! -d "$INPUT_DIR" ]; then
    echo "Error: Input directory does not exist: $INPUT_DIR"
    exit 1
fi

# Check if groups file exists
if [ ! -e "$GROUPS_FILE" ]; then
    echo "Error: Groups file not found: $GROUPS_FILE"
    exit 1
fi

# Create output directory
mkdir -p "$OUTPUT_DIR"
mkdir -p "$OUTPUT_DIR/lists"

# Path to the merge script
SCRIPT_DIR="$(dirname "$0")"
MERGE_SCRIPT="$SCRIPT_DIR/merge_bam_files.sh"

# Make sure the merge script is executable
chmod +x "$MERGE_SCRIPT"

# Process each group
echo "=== Merging BAM Files by Groups ==="
echo "Input directory: $INPUT_DIR"
echo "Output directory: $OUTPUT_DIR"
echo "Groups file: $GROUPS_FILE"
echo

while IFS=$'\t' read -r group_name pattern_or_prefix; do
    # Skip comments and empty lines
    [[ "$group_name" =~ ^#.*$ || -z "$group_name" ]] && continue
    
    echo "Processing group: $group_name"
    echo "Pattern/Prefix: $pattern_or_prefix"
    
    # Create a list file for this group
    list_file="$OUTPUT_DIR/lists/${group_name}_file_list.txt"
    
    # Find all matching BAM files for this group
    if [[ "$pattern_or_prefix" == *"*"* ]]; then
        # It's a pattern
        find "$INPUT_DIR" -name "$pattern_or_prefix" -type f | sort > "$list_file"
    else
        # It's a prefix
        find "$INPUT_DIR" -name "${pattern_or_prefix}*${PATTERN#*\*}" -type f | sort > "$list_file"
    fi
    
    # Check if any files were found
    file_count=$(wc -l < "$list_file")
    if [ $file_count -eq 0 ]; then
        echo "Warning: No files found for group $group_name with pattern/prefix $pattern_or_prefix"
        continue
    fi
    
    echo "Found $file_count files for group $group_name"
    
    # Output merged BAM file name
    output_bam="$OUTPUT_DIR/${group_name}-${SUFFIX}.bam"
    
    # Run the merge script
    merge_cmd="$MERGE_SCRIPT --list \"$list_file\" --output \"$output_bam\""
    if [ $CREATE_INDEX -eq 1 ]; then
        merge_cmd+=" --index"
    fi
    
    echo "Running merge command: $merge_cmd"
    eval $merge_cmd
    
    # Check result
    if [ $? -eq 0 ]; then
        echo "Group $group_name processed successfully"
    else
        echo "Error processing group $group_name"
    fi
    
    echo
done < "$GROUPS_FILE"

echo "=== Group merging completed ==="
echo "All merged BAM files are in: $OUTPUT_DIR" 