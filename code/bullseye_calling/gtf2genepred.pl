#!/usr/bin/env perl
#
# GTF to GenePred/refFlat Converter
# This script converts GTF format files to GenePred/refFlat format for use with Bullseye
#

use strict;
use warnings;
use Getopt::Long;
use IO::Zlib;
use File::Basename;

# Process command line arguments
my $gtf_file = '';
my $output_file = '';
my $help = 0;

GetOptions(
    "gtf=s" => \$gtf_file,
    "out=s" => \$output_file,
    "help" => \$help
) or die("Error in command line arguments\n");

# Display usage
if ($help || !$gtf_file || !$output_file) {
    print "Usage: $0 --gtf <gtf_file> --out <output_file>\n";
    print "Options:\n";
    print "  --gtf FILE    Input GTF file (can be gzipped)\n";
    print "  --out FILE    Output refFlat file\n";
    print "  --help        Display this help message\n";
    exit(1);
}

# Check if the GTF file exists
if (!-e $gtf_file) {
    die "Error: GTF file not found: $gtf_file\n";
}

# Data structures to store the GTF information
my %transcripts;
my %genes;
my %exons;
my %gene_names;
my %gene_ids;

# Open the GTF file
my $fh;
if ($gtf_file =~ /\.gz$/) {
    $fh = IO::Zlib->new($gtf_file, "rb") or die "Cannot open gzipped file $gtf_file: $!\n";
} else {
    open($fh, "<", $gtf_file) or die "Cannot open file $gtf_file: $!\n";
}

print "Parsing GTF file: $gtf_file\n";

# Parse the GTF file
while (my $line = <$fh>) {
    chomp $line;
    
    # Skip comment lines
    next if $line =~ /^#/;
    
    # Parse GTF line
    my ($chr, $source, $feature, $start, $end, $score, $strand, $frame, $attributes) = split(/\t/, $line);
    
    # Skip if not a valid feature
    next unless ($feature eq 'exon' || $feature eq 'transcript' || $feature eq 'gene');
    
    # Parse attributes
    my %attr;
    while ($attributes =~ /(\w+)\s+"([^"]+)";/g) {
        $attr{$1} = $2;
    }
    
    # Check for required attributes
    my $gene_id = $attr{'gene_id'} || "";
    my $transcript_id = $attr{'transcript_id'} || "";
    my $gene_name = $attr{'gene_name'} || $attr{'gene_id'} || "";
    
    if (!$gene_id || !$transcript_id) {
        warn "Warning: Missing gene_id or transcript_id in line: $line\n";
        next;
    }
    
    # Store gene name
    $gene_names{$gene_id} = $gene_name;
    $gene_ids{$transcript_id} = $gene_id;
    
    # Store feature information
    if ($feature eq 'exon') {
        push @{$exons{$transcript_id}}, [$start, $end];
    } elsif ($feature eq 'transcript') {
        $transcripts{$transcript_id} = {
            'chr' => $chr,
            'start' => $start,
            'end' => $end,
            'strand' => $strand,
            'gene_id' => $gene_id
        };
    } elsif ($feature eq 'gene') {
        $genes{$gene_id} = {
            'chr' => $chr,
            'start' => $start,
            'end' => $end,
            'strand' => $strand,
            'name' => $gene_name
        };
    }
}

close($fh);

# Sort exons for each transcript
foreach my $tid (keys %exons) {
    # Sort exons by start position
    my @sorted_exons = sort { $a->[0] <=> $b->[0] } @{$exons{$tid}};
    $exons{$tid} = \@sorted_exons;
}

# Open output file
open(my $out_fh, ">", $output_file) or die "Cannot open output file $output_file: $!\n";

print "Generating refFlat file: $output_file\n";

# Generate refFlat output
foreach my $tid (sort keys %transcripts) {
    next unless exists $exons{$tid};
    
    my $t = $transcripts{$tid};
    my $gid = $t->{'gene_id'};
    my $gene_name = $gene_names{$gid} || $gid;
    
    # Collect exon start and end positions
    my @exon_starts;
    my @exon_ends;
    
    foreach my $exon (@{$exons{$tid}}) {
        push @exon_starts, $exon->[0] - 1; # Convert to 0-based
        push @exon_ends, $exon->[1];
    }
    
    # Print refFlat format line
    print $out_fh join("\t",
        $gene_name,                      # gene name
        $tid,                            # transcript name/ID
        $t->{'chr'},                     # chromosome
        $t->{'strand'},                  # strand
        $exon_starts[0],                 # transcript start (0-based)
        $exon_ends[-1],                  # transcript end
        $exon_starts[0],                 # coding start (placeholder)
        $exon_ends[-1],                  # coding end (placeholder)
        scalar(@exon_starts),            # exon count
        join(',', @exon_starts) . ',',   # exon starts (comma-separated)
        join(',', @exon_ends) . ','      # exon ends (comma-separated)
    ) . "\n";
}

close($out_fh);

print "Conversion completed successfully.\n";
exit(0); 