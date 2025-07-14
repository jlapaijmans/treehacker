<!-- badges: start -->

[![Treehacker CI](https://github.com/jlapaijmans/treehacker/actions/workflows/treehacker_test.yml/badge.svg)](https://github.com/jlapaijmans/treehacker/actions/workflows/treehacker_test.yml)
<!-- badges: end -->


# TreeHacker v1.0

This tool was developed by Johanna Paijmans and Axel Barlow to do topology tests in sliding windows along the genome. The method was originally developed in [Barlow et al 2018](https://doi.org/10.1038/s41559-018-0654-8) (Nature Ecology & Evolution).

**NEW in v1.0:** Enhanced with parallel processing, comprehensive file checking, command-line options, and somewhat improved user experience!

TreeHacker v1.0 is designed to be user-friendly with comprehensive validation, clear error messages, and flexible command-line options. Feel free to contact us if you have any questions.

You are free to use treehacker and adapt it to your needs, if you do please cite:

*Lucena-Perez M, et al (2024)* 

And/or the original development of the method:

*Barlow A, et al (2018) Partial genomic survival of cave bears in living brown bears. Nature Ecology and Evolution 2: 1563–1570. https://doi.org/10.1038/s41559-018-0654-8*

## Quick Start

```bash
# Basic usage with default settings
./TREEhacker_1.0.sh fastafiles.txt results 1000 1000 0.3 BOTH

# With custom parallel processing
./TREEhacker_1.0.sh --parallel 8 --threads 4 fastafiles.txt results 1000 1000 0.3 DNA

# Get help
./TREEhacker_1.0.sh --help
```

## Key Features

- **Parallel Processing**: Configurable parallel RAxML jobs for faster analysis
- **File Checking**: Validates some of the input files and dependencies before running
- **Flexible Analysis**: Support for DNA-only, Binary-only, or both analyses
- **Missing Data Filtering**: Configurable missing data thresholds per window

## Installation & Dependencies

### Required Software

All of these tools must be installed and available in your PATH:

- **samtools** - For sequence indexing and extraction
- **bedtools** - For generating sliding windows  
- **bioawk** - For FASTA processing
- **seqtk** - For sequence manipulation
- **raxml-ng** - For phylogenetic inference (replaces older RAxML)
- **bc** - For mathematical calculations
- **nw_utils** (newick-utils) - For tree topology analysis

### Installation Examples

```bash
# Using Conda (untested, use with caution)
conda install -c bioconda samtools bedtools seqtk bioawk raxml-ng bc newick-utils

# Using Ubuntu/Debian
sudo apt install samtools bedtools seqtk bioawk bc
# Note: You will need to install newick-utils and raxml separately
```

## Usage

```
./TREEhacker_1.0.sh [OPTIONS] <fastafiles> <output> <winsize> <stepsize> <prop_missing> [analysis_type]
```

### Parameters

- **fastafiles**: File containing full filenames of fasta files (one per line)
  - Example: `ls *.fasta > fastafiles.txt`
  - **Important**: Outgroup must be the last entry in this file
- **output**: Output base name for all generated files
- **winsize**: Size of sliding windows (in base pairs)
- **stepsize**: Step size for sliding windows (in base pairs)
  - For non-overlapping windows: stepsize = winsize
- **prop_missing**: Maximum proportion of missing data allowed per sample (0.0 to 1.0)
- **analysis_type** [optional]: DNA, BIN, or BOTH (default: BOTH)

### Options

- `--help, -h`: Show detailed help message
- `--threads N`: Number of threads per RAxML job (default: 2)
- `--parallel N`: Number of parallel RAxML jobs (default: 4)

### Examples

```bash
# Basic analysis with 1000bp non-overlapping windows, 30% missing data allowed
./TREEhacker_1.0.sh fastafiles.txt myproject 1000 1000 0.3

# DNA analysis only with high-performance settings
./TREEhacker_1.0.sh --parallel 8 --threads 4 fastafiles.txt dna_analysis 500 500 0.2 DNA

# Binary analysis with overlapping windows
./TREEhacker_1.0.sh fastafiles.txt binary_trees 1000 500 0.4 BIN
```

## Notes

- The script includes some file validation and should provide clear error messages if files are missing or dependencies are not installed
- Performance can be optimized by adjusting `--parallel` (number of concurrent RAxML jobs) and `--threads` (CPU threads per job)
- For large datasets, monitor memory usage when increasing parallel jobs

## Description

TreeHacker takes a number of haploidised genomes in fasta format, and then moves along the genome in non-overlapping sliding windows to generate maximum likelihood trees using RaXML (ref). The TreeHacker tool currently does a tree from the DNA data directly using the GTRGAMMA model, and also converts the alignment to a binary format for running the phylogeny in a BINGAMMA model.

## Input
$1 	   fastanames.txt: basenames fasta files. Outgroup needs to be last individual
E.g. `ls fasta.fa | cut -d'.' -f1 > fastanames.txt`
$2 	  output base name
$3	window size. TreeHacker is currently hardcoded for non-overlapping windows

## Tool functionality in bullet points:
1. Generates windows of X bp (user defined)
2. Remove tail ends of the scaffolds/chromosomes that are less than X bp
3. Generates report of the number of missing data (N’s) in each window: "$outname"_missing_data_report.txt
4. Creates alignment files of each window, that contains one specific window from each individual in a single file
5. Runs RAXML with GTRGAMMA model
6. Converts the alignment file to binary format
7. Runs RAXML with BINGAMMA model
8. Creates a summary from all windows that gives all topologies and the number of times that topology occurred, for both the DNA and BIN alignments: "$outname"_RAxML_bestTrees.txt and "$outname"_RAxML_bestTrees_BIN.txt

## Current vs. Legacy Versions

**TreeHacker v1.0** includes major improvements over the original version while maintaining the same core functionality:

### New Features in v1.0
- **Parallel Processing**: Configurable parallel RAxML jobs (--parallel option)
- **Threading Control**: Adjustable threads per job (--threads option)
- **Comprehensive Validation**: Pre-run checks for files and dependencies
- **Flexible Analysis Types**: Choose DNA, BIN, or BOTH analyses
- **More User-Friendly Interface**: Clear help messages and error reporting
- **Modern Tools**: Uses RAxML-NG instead of legacy RAxML

### Preserved from Original
- Core sliding window algorithm
- GTR+G and BIN+G phylogenetic models  
- Binary sequence conversion (A/G=0, C/T=1)
- Topology frequency analysis
- Missing data filtering
- Multi-individual alignment processing

The legacy input format and hardcoded variables have been replaced with a flexible command-line interface while maintaining full backward compatibility of the analysis methodology.

