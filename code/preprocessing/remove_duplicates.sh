#!/usr/bin/bash
  
# Script to remove duplicates from BAM files using Picard

# Path to input and output directories
input_dir="path/to/input/bams"
output_dir="path/to/output/bams/rm_dup"

# Path to Picard
picard="path/to/picard.jar"

# Ensure output directory exists; create a new directory if it doesn't
mkdir -p "$output_dir"

# Counter for checkpoint
count=0

# Loop over input directory and search for bam files only
for file in "$input_dir"/*out.bam; do
    if [ -e "$file" ]; then

        # Extract output name
        name=$(basename $file | cut -d_ -f1-4)
        ((count++))
        echo "Mark duplicate on file: $name"

        # Run picard - mark duplicates
        java -jar $picard MarkDuplicates \
            --INPUT $file \
            --OUTPUT ${output_dir}/${name}_Aligned.sortedByCoord.rm_dup.bam \
            --METRICS_FILE ${output_dir}/${name}_Aligned.sortedByCoord.rm_dup.csv \
            --REMOVE_DUPLICATES true \
            --REMOVE_SEQUENCING_DUPLICATES true
    fi
done

echo "Duplicate removal completed for all files." 