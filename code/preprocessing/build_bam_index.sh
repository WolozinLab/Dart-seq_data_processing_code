#!/usr/bin/bash
  
# Script to build index for BAM files

# Path of file directory
input_dir="path/to/bam/files"

# Load samtools if needed (comment out if samtools is already in PATH)
# module load samtools

# Loop over files in directory
for file in "$input_dir"/*.bam; do
    # Check if file exists
    if [ -e "$file" ]; then
        name=$(basename "$file")
        echo "Currently indexing file: $name"
        
        # Index BAM files
        samtools index $file
        
        echo "Indexing completed for: $name"
    fi
done

echo "All BAM files have been indexed successfully." 