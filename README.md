# Dart-seq Data Processing Code

This repository contains code for processing Dart-seq data, a method for detecting RNA modifications.

## Project Structure

- `code/` - Contains all data processing code
  - `preprocessing/` - Code for data preprocessing
    - `STAR_alignment.sh` - Script for aligning reads using STAR
    - `Trimmomatic.sh` - Script for trimming reads using Trimmomatic
    - `remove_duplicates.sh` - Script for removing duplicates from BAM files
    - `build_bam_index.sh` - Script for building indices for BAM files
    - `merge_bam_files.sh` - Script for merging multiple BAM files
    - `merge_bam_by_groups.sh` - Script for merging BAM files by condition groups
    - `sample_groups.txt` - Example configuration for defining sample groups
  - `bullseye_calling/` - Code for Bullseye calling
    - `run_parseBAM.sh` - Script for parsing BAM files to matrix format
    - `batch_parseBAM.sh` - Batch processing script for parsing multiple BAM files
    - `run_find_edit_site.sh` - Script for identifying m6A sites
    - `batch_find_edit_sites.sh` - Batch processing script for identifying m6A sites in multiple files
    - `build_reference.sh` - Script for creating reference files
    - `gtf2genepred.pl` - Perl script for converting GTF to genepred format
  - `DRACH_filter.sh` - Script for filtering sites by DRACH motif
  - `merge_filter_sites.R` - R script for merging and filtering Bullseye sites

## Data Processing Pipeline

The Dart-seq data processing pipeline consists of the following steps:

1. **Quality Control and Trimming**: Raw reads are processed using Trimmomatic to remove adapters and low-quality sequences.
2. **Alignment**: Trimmed reads are aligned to the reference genome using STAR aligner.
3. **Duplicate Removal**: PCR and sequencing duplicates are removed from aligned BAM files using Picard.
4. **BAM Indexing**: BAM files are indexed using samtools for faster access.
5. **Sample Merging** (optional): BAM files can be merged by condition/group for combined analysis.
6. **Bullseye Calling**: Processed BAM files are analyzed to identify Bullseye sites, which indicate RNA modifications.
7. **DRACH Filtering**: Called sites are filtered to retain those matching the DRACH motif (where D=A/G/U, R=A/G, A=A, C=C, H=A/C/U).
8. **Results Merging**: Results from multiple samples can be merged and summarized for downstream analysis.

## Installation

For detailed installation instructions of the basic pipeline components, please refer to the [Installation Guide](INSTALL.md).

### Bullseye Pipeline Installation

The Bullseye calling pipeline requires additional dependencies. The easiest way to install them is through conda:

#### Prerequisites

- Perl > 5.26
- Samtools
- Bedtools
- Tabix
- Perl modules: MCE, Math::CDF, Bio::DB::Fasta, and others
- R packages: GenomicRanges, BSgenome, Biostrings, gUtils, dplyr, purrr, stringr

#### Installation Steps

1. Create and activate the conda environment with the provided yml file:

```bash
conda env create -f bullseye.yml
conda activate Bullseye
```

2. Install XML::Parser (required for Bio::DB::Fasta):

```bash
wget http://www.cpan.org/authors/id/T/TO/TODDR/XML-Parser-2.46.tar.gz 
tar -xf XML-Parser-2.46.tar.gz
cd XML-Parser-2.46 
perl Makefile.PL EXPATLIBPATH=$CONDA_PREFIX/lib EXPATINCPATH=$CONDA_PREFIX/include
make
make install
```

3. Install the remaining Perl packages:

```bash
cpanm Bio::DB::Fasta
cpanm Text::NSP
cpanm Array::IntSpan
cpanm MCE
```

4. Install required R packages:

```bash
R -e 'if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager"); BiocManager::install(c("GenomicRanges", "BSgenome", "Biostrings", "gUtils"))'
R -e 'install.packages(c("dplyr", "purrr", "stringr"))'
```

## Usage

### Data Preprocessing Steps

#### STAR Alignment

To align reads to a reference genome using STAR:

```bash
./code/preprocessing/STAR_alignment.sh {read1.fastq.gz} {read2.fastq.gz} {reference-genome} {output-dir}
```

#### Trimmomatic

To trim reads using Trimmomatic:

```bash
./code/preprocessing/Trimmomatic.sh {input-R1.fastq.gz} {input-R2.fastq.gz} {output-paired-R1.fastq.gz} {output-paired-R2.fastq.gz} {output-unpaired-R1.fastq.gz} {output-unpaired-R2.fastq.gz}
```

#### Remove Duplicates

To remove PCR duplicates from BAM files:

```bash
./code/preprocessing/remove_duplicates.sh
```
Note: Before running this script, modify the `input_dir` and `output_dir` variables to point to your input and output directories, and the `picard` variable to point to your Picard JAR file.

#### Build BAM Index

To create indices for BAM files:

```bash
./code/preprocessing/build_bam_index.sh
```
Note: Before running this script, modify the `input_dir` variable to point to your BAM files directory and ensure samtools is in your PATH.

### Merging BAM Files

#### Basic BAM Merging

To merge multiple BAM files into a single file:

```bash
./code/preprocessing/merge_bam_files.sh --list <file_list.txt> --output <output.bam> [--index]
```

Required arguments:
- `--list`, `-l`: File containing a list of BAM files to merge (one per line)
- `--output`, `-o`: Output merged BAM file path

Additional options:
- `--index`, `-i`: Create index for the merged BAM file

#### Merging BAM Files by Groups

To merge multiple BAM files by condition groups (e.g., disease vs. control):

```bash
./code/preprocessing/merge_bam_by_groups.sh --input-dir <bam_dir> --output-dir <output_dir> --groups <groups_file.txt> [options]
```

Required arguments:
- `--input-dir`, `-i`: Directory containing BAM files
- `--output-dir`, `-o`: Output directory for merged BAM files
- `--groups`, `-g`: Tab-delimited file defining groups and their samples

Additional options:
- `--pattern`, `-p`: Pattern for BAM files (default: "*.bam")
- `--suffix`, `-s`: Suffix to add to output files (default: "combined")
- `--index`, `-x`: Create index for merged BAM files
- `--help`, `-h`: Display help message

Example groups file format:
```
# group_name <tab> pattern_or_prefix
DART-mut_AD     *DART-mut*AD*.bam
DART-mut_healthy     *DART-mut*healthy*.bam
```

### Bullseye Analysis

#### Parsing BAM Files

To parse a BAM file into a matrix format required for Bullseye:

```bash
./code/bullseye_calling/run_parseBAM.sh <input.bam> <output_prefix>
```

For batch processing multiple BAM files:

```bash
./code/bullseye_calling/batch_parseBAM.sh --input-dir <bam_directory> --output-dir <output_directory> [options]
```

Required arguments:
- `--input-dir`, `-i`: Directory containing BAM files
- `--output-dir`, `-o`: Output directory for matrix files

Additional options:
- `--pattern`, `-p`: Pattern to match BAM files (default: "*.bam")
- `--cpu`, `-c`: Number of CPU cores to use (default: 4)
- `--min-coverage`, `-m`: Minimum coverage threshold (default: 10)
- `--no-remove-duplicates`: Do not remove duplicates
- `--help`, `-h`: Display help message

#### Finding Edit Sites

To identify m6A sites from a matrix file:

```bash
./code/bullseye_calling/run_find_edit_site.sh <control-file> <edit-file> <output> <annotation-file>
```

For batch processing multiple matrix files:

```bash
./code/bullseye_calling/batch_find_edit_sites.sh --control <control.matrix> --input-dir <matrix_directory> --output-dir <output_directory> --annotation <annotation.refFlat> [options]
```

Required arguments:
- `--control`, `-c`: Control matrix file
- `--input-dir`, `-i`: Directory containing edited matrix files
- `--output-dir`, `-o`: Output directory for results
- `--annotation`, `-a`: Annotation file in refFlat format

Additional options:
- `--pattern`, `-p`: Pattern to match matrix files (default: "*.matrix")
- `--cpu`, `-t`: Number of CPU cores to use (default: 4)
- `--min-edit`: Minimum edit percent (default: 10)
- `--max-edit`: Maximum edit percent (default: 95)
- `--fold`: Edit fold threshold (default: 1.5)
- `--min-sites`: Minimum edit sites (default: 2)
- `--control-cov`: Control minimum coverage (default: 10)
- `--edited-cov`: Edited minimum coverage (default: 10)
- `--intron`: Include introns in analysis
- `--ext-utr`: Extended UTR size (for protein coding genes)
- `--help`, `-h`: Display help message

### Building Reference Files

To build reference files needed for Bullseye analysis:

```bash
./code/bullseye_calling/build_reference.sh --gtf <annotation.gtf> --output <output_directory> [options]
```

Required arguments:
- `--gtf`, `-g`: GTF annotation file
- `--output`, `-o`: Output directory

Additional options:
- `--method`, `-m`: Method to convert GTF to refFlat [ucsc|custom] (default: custom)
- `--help`, `-h`: Display help message

### DRACH Filtering

To filter Bullseye sites by the DRACH motif:

```bash
./code/DRACH_filter.sh <input.bed> <output_directory>
```

Required arguments:
- `input.bed`: Input BED file containing Bullseye sites
- `output_directory`: Directory to store the filtered results

The script will produce two output files:
- `RAC.bed`: Sites matching the RAC motif (R=A/G, A=A, C=C)
- `DRACH.bed`: Sites matching the stricter DRACH motif (D=A/G/U, R=A/G, A=A, C=C, H=A/C/U)

### Merging and Filtering Results

To merge and filter results from multiple Bullseye runs:

```bash
Rscript code/merge_filter_sites.R <bed_directory> <output_directory> <genome_path>
```

Required arguments:
- `bed_directory`: Directory containing BED files from Bullseye
- `output_directory`: Directory to store the merged and filtered results
- `genome_path`: Path to the reference genome FASTA file

This script:
1. Merges multiple BED files from Bullseye analysis
2. Filters sites to retain those matching the DRACH motif
3. Generates various output files for downstream analysis:
   - `Sum_of_all_no_DRACH_filter.csv`: All merged sites without filtering
   - `Sum_of_all_filter_by_DRACH.csv`: Sites filtered by DRACH motif
   - `matrix_Col1.csv`: Gene information for filtered sites
   - `matrix_Col6.csv`: Editing ratio information for filtered sites

## Suggested Workflow

While there is no single workflow script, here's a recommended sequence of steps to process Dart-seq data:

1. **Trim raw reads** using Trimmomatic
2. **Align trimmed reads** to the reference genome using STAR
3. **Remove duplicates** from aligned BAM files
4. **Index BAM files** for faster access
5. **Parse BAM files** to matrix format using run_parseBAM.sh or batch_parseBAM.sh
6. **Build reference files** from GTF annotation if needed
7. **Find edit sites** using run_find_edit_site.sh or batch_find_edit_sites.sh
8. **Filter sites by DRACH motif** using DRACH_filter.sh
9. **Merge and analyze results** using merge_filter_sites.R

Users can create their own workflow scripts by combining these individual steps according to their specific requirements.

## Requirements

The pipeline requires the following tools to be installed:

- Trimmomatic (v0.39 or later)
- STAR (v2.7 or later)
- Picard Tools
- Samtools
- Bedtools
- Python 3.6 or later
- Perl 5.26 or later
- Perl modules: MCE, Math::CDF, Bio::DB::Fasta
- R 4.0 or later with packages:
  - GenomicRanges
  - BSgenome
  - Biostrings
  - gUtils
  - dplyr
  - purrr
  - stringr

See the [Installation Guide](INSTALL.md) for detailed installation instructions.

## Changelog

- Initial project structure created
- Added STAR alignment script
- Added Trimmomatic script
- Added duplicate removal script
- Added BAM indexing script 
- Added BAM file merging functionality
- Added group-based BAM merging
- Added parseBAM scripts for matrix generation
- Added edit site finding scripts
- Added reference building utilities
- Added DRACH filter script
- Added R script for merging and filtering results
- Removed dart_seq_workflow.sh in favor of modular approach
