# Dart-seq Data Processing Code

This repository contains code for processing Dart-seq data, a method for detecting RNA modifications.

## Project Structure

- `code/` - Contains all data processing code
  - `preprocessing/` - Code for data preprocessing
    - `STAR_alignment.sh` - Script for aligning reads using STAR
    - `Trimmomatic.sh` - Script for trimming reads using Trimmomatic
    - `remove_duplicates.sh` - Script for removing duplicates from BAM files
    - `build_bam_index.sh` - Script for building indices for BAM files
  - `bullseye_calling/` - Code for Bullseye calling
    - `bullseye_caller.py` - Python script for identifying Bullseye sites from BAM files
  - `dart_seq_workflow.sh` - Master workflow script that runs the entire pipeline

## Data Processing Pipeline

The Dart-seq data processing pipeline consists of the following steps:

1. **Quality Control and Trimming**: Raw reads are processed using Trimmomatic to remove adapters and low-quality sequences.
2. **Alignment**: Trimmed reads are aligned to the reference genome using STAR aligner.
3. **Duplicate Removal**: PCR and sequencing duplicates are removed from aligned BAM files using Picard.
4. **BAM Indexing**: BAM files are indexed using samtools for faster access.
5. **Bullseye Calling**: Processed BAM files are analyzed to identify Bullseye sites, which indicate RNA modifications.

## Usage

### Complete Workflow

To run the entire Dart-seq data processing pipeline:

```
./code/dart_seq_workflow.sh --read1 <read1.fastq.gz> --read2 <read2.fastq.gz> --genome <reference.fa> --output <output_dir> [options]
```

Required arguments:
- `--read1`, `-r1`: Read1 FASTQ file
- `--read2`, `-r2`: Read2 FASTQ file
- `--genome`, `-g`: Reference genome file
- `--output`, `-o`: Output directory

Additional options:
- `--min-coverage`, `-c`: Minimum coverage for Bullseye (default: 10)
- `--min-quality`, `-q`: Minimum mapping quality (default: 20)
- `--help`, `-h`: Display help message

The workflow automatically handles all processing steps from raw reads to Bullseye calling.

### Individual Steps

If you prefer to run individual steps of the pipeline separately, you can use the following scripts:

#### STAR Alignment
```
./code/preprocessing/STAR_alignment.sh {read1.fastq.gz} {read2.fastq.gz} {reference-genome} {output-dir}
```

#### Trimmomatic
```
./code/preprocessing/Trimmomatic.sh {input-R1.fastq.gz} {input-R2.fastq.gz} {output-paired-R1.fastq.gz} {output-paired-R2.fastq.gz} {output-unpaired-R1.fastq.gz} {output-unpaired-R2.fastq.gz}
```

#### Remove Duplicates
```
./code/preprocessing/remove_duplicates.sh
```
Note: Before running this script, modify the `input_dir` and `output_dir` variables to point to your input and output directories, and the `picard` variable to point to your Picard JAR file.

#### Build BAM Index
```
./code/preprocessing/build_bam_index.sh
```
Note: Before running this script, modify the `input_dir` variable to point to your BAM files directory and ensure samtools is in your PATH.

#### Bullseye Caller
```
python ./code/bullseye_calling/bullseye_caller.py --input <input.bam> --output <output_dir> --genome <reference.fa> [options]
```

Optional arguments:
- `--min-coverage`, `-c`: Minimum coverage threshold (default: 10)
- `--min-quality`, `-q`: Minimum mapping quality (default: 20)

## Requirements

The pipeline requires the following tools to be installed:

- Trimmomatic (v0.39 or later)
- STAR (v2.7 or later)
- Picard Tools
- Samtools
- Python 3.6 or later

## Changelog

- Initial project structure created
- Added STAR alignment script
- Added Trimmomatic script
- Added Bullseye caller script
- Added complete workflow script
- Added duplicate removal script
- Added BAM indexing script
