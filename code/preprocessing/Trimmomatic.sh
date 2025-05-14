#!/bin/bash -l
  
## NOTE: UNTESTED
## Syntax: ./Trimmomatic.sh {input-R1} {input-R2} {output-paired-R1} {output-paired-R2} {output-unpaired-R1} {output-unpaired-R2}

# Path to Trimmomatic
Trimmomatic="path/to/Trimmomatic/trimmomatic.jar"
adapter="path/to/Trimmomatic/adapters/TruSeq3-PE.fa"

# Run Trimmomatic
java -jar $Trimmomatic PE \
        ${1} ${2} ${3} ${4} ${5} ${6} \
        ILLUMINACLIP:${adapter}:2:30:10:8:true #\
        #LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:36 