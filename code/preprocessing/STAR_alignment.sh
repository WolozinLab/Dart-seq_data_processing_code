#!/bin/bash -l
  
## Syntax: STAR_alignment.sh {read1} {read2} {reference-genome} {output-dir}

## Load STAR
STAR="path/to/STAR/STAR"

## Extract file name
name=$(basename ${1} | cut -d_ -f1-4)

$STAR --genomeDir ${3} \
    --readFilesCommand gunzip -c \
    --readFilesIn ${1} ${2} \
    --outFileNamePrefix ${4}/${name}_ \
    --outSAMtype BAM SortedByCoordinate

echo "Finished running STAR on $name" 