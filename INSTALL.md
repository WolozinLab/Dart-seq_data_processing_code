# Installation Guide for Dart-seq Data Processing Pipeline

This guide provides instructions for installing and configuring the tools required for the Dart-seq data processing pipeline.

## Prerequisites

- Unix-like operating system (Linux or macOS)
- At least 16GB RAM (32GB or more recommended for large genomes)
- Sufficient disk space for sequence data and analysis results
- Administrative privileges or ability to install software

## Required Software

### 1. Trimmomatic

Trimmomatic is used for trimming and quality filtering of sequencing reads.

```bash
# Download Trimmomatic v0.39
wget http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.39.zip
unzip Trimmomatic-0.39.zip
```

### 2. STAR Aligner

STAR is used for aligning sequencing reads to the reference genome.

```bash
# Clone STAR repository
git clone https://github.com/alexdobin/STAR.git
cd STAR/source

# Compile
make STAR

# Add to PATH (optional)
export PATH=$PATH:$(pwd)
```

### 3. Picard Tools

Picard is used for processing BAM files and removing duplicates.

```bash
# Download Picard
wget https://github.com/broadinstitute/picard/releases/download/2.27.5/picard.jar
```

### 4. Samtools

Samtools is used for manipulating SAM/BAM files and building indices.

```bash
# Install dependencies
sudo apt-get install libncurses5-dev libbz2-dev liblzma-dev # Debian/Ubuntu
# OR
brew install xz ncurses # macOS with Homebrew

# Download and install samtools
wget https://github.com/samtools/samtools/releases/download/1.16.1/samtools-1.16.1.tar.bz2
tar -xjf samtools-1.16.1.tar.bz2
cd samtools-1.16.1

./configure --prefix=/usr/local
make
sudo make install
```

### 5. Python 3

Python 3.6 or later is required for running analysis scripts.

```bash
# Install Python 3 on Debian/Ubuntu
sudo apt-get install python3 python3-pip

# OR on macOS with Homebrew
brew install python

# Install required Python packages
pip3 install numpy pandas pysam matplotlib
```

### 6. Bedtools

Bedtools is required for DRACH filtering and genomic feature analysis.

```bash
# Install on Debian/Ubuntu
sudo apt-get install bedtools

# OR on macOS with Homebrew
brew install bedtools
```

### 7. R and Bioconductor packages

R and several Bioconductor packages are required for merging and filtering sites.

```bash
# Install R
# Debian/Ubuntu
sudo apt-get install r-base r-base-dev

# OR macOS with Homebrew
brew install r

# Install required R packages
R -e 'if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager"); BiocManager::install(c("GenomicRanges", "BSgenome", "Biostrings", "gUtils"))'
R -e 'install.packages(c("dplyr", "purrr", "stringr"))'
```

## Configuring the Pipeline

After installing the required software, you need to update the path variables in the scripts:

1. Edit `code/preprocessing/Trimmomatic.sh` to set the correct path for Trimmomatic:
   ```bash
   Trimmomatic="path/to/Trimmomatic-0.39/trimmomatic-0.39.jar"
   adapter="path/to/Trimmomatic-0.39/adapters/TruSeq3-PE.fa"
   ```

2. Edit `code/preprocessing/STAR_alignment.sh` to set the correct path for STAR:
   ```bash
   STAR="path/to/STAR/source/STAR"
   ```

3. Edit `code/preprocessing/remove_duplicates.sh` to set the correct paths:
   ```bash
   picard="path/to/picard.jar"
   ```

4. Edit `code/preprocessing/build_bam_index.sh` to ensure samtools is in your PATH or uncomment and modify the module loading line if needed.

5. Ensure `code/DRACH_filter.sh` has the correct path to your reference genome or set the `GENOME_PATH` environment variable.

## Troubleshooting

- **Java errors**: Ensure you have Java 8 or later installed.
- **Memory issues**: Adjust memory allocation for Java tools by modifying the `-Xmx` parameter.
- **Path issues**: Verify that all paths in the scripts point to the correct locations of the installed tools.
- **File permissions**: Ensure the script files are executable (`chmod +x code/*.sh code/preprocessing/*.sh code/bullseye_calling/*.sh`).

For additional help, please open an issue in the repository. 