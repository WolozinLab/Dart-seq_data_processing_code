#!/bin/bash -l
#
# Run Find_edit_site.pl script
# This script compares editing in a matrix file from parseBAM.pl to a control matrix
# (Mettl3 KO or YTHmut-APOBEC1) or to the genomic sequence to identify m6A sites
#

# Check if the input arguments are provided
if [ $# -lt 4 ]; then
    echo "Syntax: $0 {control-file} {edit-file} {output} {annotation-file}"
    echo "Example: $0 control.matrix edited.matrix results.out annotation.refFlat"
    exit 1
fi

# Activate Bullseye using conda
module load miniconda
conda activate Bullseye

# Run Find_edit_site.pl
perl Find_edit_site.pl --annotationFile ${4} \
        --controlMatrix ${1} \
        --EditedMatrix ${2} \
        --outfile ${3} \
        --minEdit 10 \
        --maxEdit 95 \
        --editFoldThreshold 1.5 \
        --MinEditSites 2 \
        --cpu 26 \
        --ControlMinCoverage 10 \
        --EditedMinCoverage 10 \
        --verbose

# Check if the command was successful
if [ $? -eq 0 ]; then
    echo "Find_edit_site.pl completed successfully. Output file: ${3}"
else
    echo "Error: Find_edit_site.pl failed to complete."
    exit 1
fi 