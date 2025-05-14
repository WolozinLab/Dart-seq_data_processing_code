#!/bin/bash -l
#
# DRACH Filter Script
# This script filters Bullseye sites by the DRACH motif (where D=A/G/U, R=A/G, A=A, C=C, H=A/C/U)
#

# Check if input arguments are provided
if [ $# -lt 2 ]; then
    echo "Usage: $0 <input.bed> <output_directory>"
    echo "Example: $0 unfilter_raw_bullseye.bed /path/to/output"
    exit 1
fi

INPUT_BED=$1
OUTPUT_DIR=$2

# Make sure the output directory exists
mkdir -p "$OUTPUT_DIR"

# Determine the genome path (use default or from environment if available)
if [ -z "$GENOME" ]; then
    # Check if a genome path is set in environment
    if [ -z "$GENOME_PATH" ]; then
        # Use a default path - update this with your default genome path
        GENOME="/restricted/projectnb/benwol/jmh/script/ref/human/hg38_111/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa"
        echo "Using default genome: $GENOME"
    else
        GENOME="$GENOME_PATH"
        echo "Using genome from environment: $GENOME"
    fi
else
    echo "Using specified genome: $GENOME"
fi

# Check if the genome file exists
if [ ! -f "$GENOME" ]; then
    echo "Error: Genome file not found at $GENOME"
    echo "Please set a valid genome path with GENOME_PATH environment variable."
    exit 1
fi

# Set output filenames
RAC_OUT="$OUTPUT_DIR/RAC.bed"
DRACH_OUT="$OUTPUT_DIR/DRACH.bed"
TEMP_FASTA="$OUTPUT_DIR/temp.fa"
TEMP_EXTEND="$OUTPUT_DIR/temp_extend.bed"

echo "Input BED file: $INPUT_BED"
echo "Output directory: $OUTPUT_DIR"
echo "Genome file: $GENOME"
echo "RAC output: $RAC_OUT"
echo "DRACH output: $DRACH_OUT"

# Function to extract sequence context
extract_sequence_context() {
    local input_bed=$1
    local genome=$2
    local output_fasta=$3
    local extend_bed=$4
    
    echo "Extending coordinates to get sequence context..."
    # Extend BED coordinates by 2 bases on each side (for DRACH/RRACH motifs)
    awk 'BEGIN{OFS="\t"} {
        # Handle strand-specific extension
        if($6 == "+") {
            start = $2 - 2;
            end = $3 + 2;
            if(start < 0) start = 0;
            print $1, start, end, $4, $5, $6;
        } else {
            start = $2 - 2;
            end = $3 + 2;
            if(start < 0) start = 0;
            print $1, start, end, $4, $5, $6;
        }
    }' "$input_bed" > "$extend_bed"
    
    echo "Extracting sequences..."
    # Extract DNA sequences for the extended regions
    bedtools getfasta -fi "$genome" -bed "$extend_bed" -fo "$output_fasta" -s -name
}

# Function to filter for RAC and DRACH motifs
filter_for_motifs() {
    local input_bed=$1
    local fasta=$2
    local rac_out=$3
    local drach_out=$4
    
    echo "Filtering for RAC and DRACH motifs..."
    
    # Get the total number of sites
    local total_sites=$(wc -l < "$input_bed")
    echo "Total sites before filtering: $total_sites"
    
    # Extract sequences and match motifs
    paste "$input_bed" "$fasta" | \
    awk '
    BEGIN {OFS="\t"; count_rac=0; count_drach=0;}
    {
        # Extract the sequence from the FASTA part of the pasted file
        seq = toupper($NF);
        
        # For + strand, we look for motifs directly
        if($6 == "+") {
            # Get the center base (should be the A in RAC/DRACH)
            center_pos = length(seq)/2;
            center_base = substr(seq, center_pos, 1);
            
            # Check for RAC motif (R=A/G, A=A, C=C)
            if (center_base == "A" && 
                (substr(seq, center_pos+1, 1) == "A" || substr(seq, center_pos+1, 1) == "G") && 
                substr(seq, center_pos+2, 1) == "C") {
                count_rac++;
                print $0, seq >> "'$rac_out'";
                
                # Further check for DRACH motif (D=A/G/T not C, R=A/G, A=A, C=C, H=A/C/T not G)
                if ((substr(seq, center_pos-2, 1) == "A" || 
                     substr(seq, center_pos-2, 1) == "G" || 
                     substr(seq, center_pos-2, 1) == "T") && 
                    (substr(seq, center_pos+3, 1) == "A" || 
                     substr(seq, center_pos+3, 1) == "C" || 
                     substr(seq, center_pos+3, 1) == "T")) {
                    count_drach++;
                    print $0, seq >> "'$drach_out'";
                }
            }
        }
        # For - strand, we need to reverse complement for proper motif checking
        else {
            # Reverse complement logic can be complex - for simplicity we just check 
            # if the center base is T (complement of A), followed by proper RAC/DRACH pattern
            center_pos = length(seq)/2;
            center_base = substr(seq, center_pos, 1);
            
            # Check for RAC motif on - strand (complement is: GYT where Y=C/T)
            if (center_base == "T" && 
                (substr(seq, center_pos-1, 1) == "C" || substr(seq, center_pos-1, 1) == "T") && 
                substr(seq, center_pos-2, 1) == "G") {
                count_rac++;
                print $0, seq >> "'$rac_out'";
                
                # Further check for DRACH motif (complement is: HYTGD where D=A/G/T not C, H=A/C/T not G)
                if ((substr(seq, center_pos+2, 1) == "A" || 
                     substr(seq, center_pos+2, 1) == "C" || 
                     substr(seq, center_pos+2, 1) == "T") && 
                    (substr(seq, center_pos-3, 1) == "A" || 
                     substr(seq, center_pos-3, 1) == "G" || 
                     substr(seq, center_pos-3, 1) == "T")) {
                    count_drach++;
                    print $0, seq >> "'$drach_out'";
                }
            }
        }
    }
    END {
        print "RAC motifs found: " count_rac > "/dev/stderr";
        print "DRACH motifs found: " count_drach > "/dev/stderr";
    }'
    
    # Count the filtered sites
    local rac_count=$(wc -l < "$rac_out")
    local drach_count=$(wc -l < "$drach_out")
    
    echo "Sites with RAC motif: $rac_count ($(echo "scale=1; 100*$rac_count/$total_sites" | bc)%)"
    echo "Sites with DRACH motif: $drach_count ($(echo "scale=1; 100*$drach_count/$total_sites" | bc)%)"
}

# Main execution
extract_sequence_context "$INPUT_BED" "$GENOME" "$TEMP_FASTA" "$TEMP_EXTEND"
filter_for_motifs "$INPUT_BED" "$TEMP_FASTA" "$RAC_OUT" "$DRACH_OUT"

# Clean up temporary files
echo "Cleaning up temporary files..."
rm -f "$TEMP_FASTA" "$TEMP_EXTEND"

echo "DRACH filtering completed successfully!"
echo "Results are in: $RAC_OUT (RAC motifs) and $DRACH_OUT (DRACH motifs)" 