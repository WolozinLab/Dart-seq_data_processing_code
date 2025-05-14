#!/bin/bash -l
#
# Run parseBAM.pl script
# This script will parse the aligned and sorted BAM files to output a tab delimited file
# with the count of each nucleotide at each position in the genome
#

# Check if the input arguments are provided
if [ $# -lt 2 ]; then
    echo "Usage: $0 <input.bam> <output_prefix>"
    echo "Example: $0 sample.bam sample_output"
    exit 1
fi

# Load conda environment
module load miniconda
conda activate Bullseye

# Run parseBAM.pl
perl parseBAM.pl --input ${1} --output ${2}.matrix --cpu 26 \
        --minCoverage 10 \
        --removeDuplicates \
        --verbose

# Check if the command was successful
if [ $? -eq 0 ]; then
    echo "ParseBAM completed successfully. Output file: ${2}.matrix"
else
    echo "Error: ParseBAM failed to complete."
    exit 1
fi 