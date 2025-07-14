<!-- badges: start -->

[![conda-single-platform](https://github.com/treehacker/actions/workflows/treehacker_test.yml/badge.svg)](https://github.com/treehacker/actions/workflows/treehacker_test.yml)
<!-- badges: end -->


# TreeHacker

This tool was developed by Johanna Paijmans and Axel Barlow to do topology tests in sliding windows along the genome. The method was originally developed in [Barlow et al 2018](https://doi.org/10.1038/s41559-018-0654-8) (Nature Ecology & Evolution).

In its current form, it is probably not very straightforward to use out-of-the-box for other projects, but we are (slowly) working to make it more user friendly. Feel free to contact us if you have any questions.

You are free to use treehacker and adapt it to your needs, if you do please cite:

*Lucena-Perez M, et al (2024)* 

And/or the original development of the method:

*Barlow A, et al (2018) Partial genomic survival of cave bears in living brown bears. Nature Ecology and Evolution 2: 1563–1570. https://doi.org/10.1038/s41559-018-0654-8*

## Notes

### required software:
* snp-sites (Page et al 2016) 	https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5320690/
* samtools (Li et al 2009)    	  https://www.ncbi.nlm.nih.gov/pubmed/19505943
* raxml (Stamatakis 2008)  	 https://cme.h-its.org/exelixis/web/software/raxml/
* nw_utils             		 http://cegg.unige.ch/newick_utils
* bioawk               		 https://github.com/lh3/bioawk
* bedtools   	 	https://bedtools.readthedocs.io/en/latest/
* seqtk   		 https://github.com/lh3/seqtk

Change the executables for the tools TreeHacker uses (`samtools`, `snp-sites`, `raxml`, `bedtools`) in the "tools" section on line 25. Also please note that a lot of variables are hardcoded at the moment, e.g. window step size and substitution model.

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

