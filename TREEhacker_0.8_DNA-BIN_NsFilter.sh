#!/bin/bash

### This script extracts informative positions from aligned genomes (multi-fasta) and computes a tree for each window
# Hacked together by Johanna Paijmans

#### required software:
# snp-sites (Page et al 2016) 	https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5320690/
# samtools (Li et al 2009)  		https://www.ncbi.nlm.nih.gov/pubmed/19505943
# raxml (Stamatakis 2008)       https://cme.h-its.org/exelixis/web/software/raxml/
# nw_utils                      http://cegg.unige.ch/newick_utils
# bioawk                        https://github.com/lh3/bioawk
# bedtools
# seqtk

### USAGE: ./TREE_hacker_v0.1.sh $1 $2 $3
# $1		fastanames.txt 	  basenames fasta files)			e.g.--- ls fasta.fa | cut -d'.' -f1 > fastanames.txt. Outgroup needs to be last individual
# $2		ouput 		        output base name
# $3    winsize

### variables
step_size=$winsize           		# step size (equal window size for non-overlapping windows)
threads=1                 # number of threads RAXML uses
max_ns=50000

### tools
samtools=/raid6/paijmans/software/samtools-1.11/samtools     		# path to samtools executable of choice
snpsites=snp-sites    		# path to snp-sites executable
raxml=raxmlHPC-PTHREADS
bedtools=bedtools-v2.23
seqtk=seqtk


### script variables for readability, do not change
fastanames="$1"
outname="$2"
win_size="$3"            		# window size
#win_size_HR=$(numfmt --to=si --suffix=b --format=%.0f $win_size)

### preamble:
# make necessary directory
[[ -d TREEhackerFiles_"$outname"/ ]] && { echo "Output folder TREEhackerFiles_"$outname"/ exists, remove or select different folder"; exit 1; }
mkdir TREEhackerFiles_"$outname"/

### the script
##1. generate windows from each fasta using bedtools, and then removing the small scaffolds (tail ends)
for ind in $(cat $fastanames);
do
  if [[ ! -f "$ind"_"$win_size"_wins.fa ]]; then
    bioawk -c fastx '{print $name"\t0\t"length($seq)}' $ind.fa | $bedtools makewindows -b - -w $win_size -s $step_size | $bedtools getfasta -fi $ind.fa -bed - -fo stdout | $seqtk seq -L $win_size - | sed 's/:/_/g'> "$ind"_"$win_size"_wins.fa
    $samtools faidx "$ind"_"$win_size"_wins.fa
  fi;
done;

### generate a report of the number of N's in each scaffold
if [[ ! -f "$outname"_missing_data_report.txt ]]; then
(for win in $(cat *_"$win_size"_wins.fa.fai | awk '{print $1}' | sort -n | uniq);
do
  #3. concatenate individuals together (1 file per scaffold)
   (for ind in $(cat $fastanames);
  do
	samtools-v0.1.19 faidx "$ind"_"$win_size"_wins.fa $win |  $seqtk comp | awk '{print $9}';
  done;) | paste <(echo "$win") - - - - -;
done;) > "$outname"_missing_data_report.txt
fi;

### filter for windows with more than "$max_ns" N's
[[ ! -f "$outname"_filtered_wins.txt ]] && awk -v max_ns=$max_ns '{ok=1;for (x=2;x<=NF;x++) if ($x>max_ns) ok=0; if (ok) print $0}' "$outname"_missing_data_report.txt > "$outname"_filtered_wins.txt


#2. generate alignments per window
for filtered_win in $(cat "$outname"_filtered_wins.txt | awk '{print $1}');
do
  #3. concatenate individuals together (1 file per scaffold)
  for ind in $(cat $fastanames);
  do
    $samtools faidx "$ind"_"$win_size"_wins.fa $filtered_win >> TREEhackerFiles_"$outname"/$filtered_win.concat.fa;
  done;

   #4. rename sequences
   bioawk -c fastx '{print $seq}' TREEhackerFiles_"$outname"/"$filtered_win".concat.fa | (while read l; do echo $l; done) | paste $fastanames - | awk '{print ">"$1"\n"$2}' > TREEhackerFiles_"$outname"/"$filtered_win".concat.RN.fa;


  ####################################################################################
  ############ NON BINARY ############
   #6. RAXML
   [[ -f TREEhackerFiles_"$outname"/"$filtered_win".concat.RN.fa ]] && $raxml -T $threads -m GTRGAMMA -o $(tail -n1 $fastanames) -p 12345 -s TREEhackerFiles_"$outname"/"$filtered_win".concat.RN.fa -w "$(pwd)"/TREEhackerFiles_"$outname"/ -n "$filtered_win".concat.RN.fa.tree >> "$(pwd)"/"$outname"_RaXML_run.log;

  ####################################################################################
  ############ BINARY ############c
   #4. convert to BINARY
   [[ -f TREEhackerFiles_"$outname"/"$filtered_win".concat.RN.fa ]] && bioawk -c fastx '(gsub(/[AGag]/,"0",$seq)) && (gsub(/[CTct]/,"1",$seq)) { print ">"$name; print $seq }' TREEhackerFiles_"$outname"/"$filtered_win".concat.RN.fa | bioawk -c fastx '$seq ~ /N/ (gsub(/[Nn]/,"?",$seq)) { print ">"$name; print $seq }' > TREEhackerFiles_"$outname"/"$filtered_win".concat.RN.BIN.fa

   #6. RAXML BIN
   [[ -f TREEhackerFiles_"$outname"/"$filtered_win".concat.RN.BIN.fa ]] && $raxml -T $threads -m BINGAMMA -o $(tail -n1 $fastanames) -p 12345 -s TREEhackerFiles_"$outname"/"$filtered_win".concat.RN.BIN.fa -w "$(pwd)"/TREEhackerFiles_"$outname"/ -n "$filtered_win".concat.RN.BIN.fa.tree >> "$(pwd)"/"$outname"_RaXML_BIN-run.log;

done;


####################################################################################
############ SUMMARISING ############

# NON-BINARY
find ./TREEhackerFiles_"$outname"/ -type f -name 'RAxML_bestTree.*.concat.RN.fa.tree' -exec cat {} + > "$outname"_RAxML_bestTrees.txt
nw_topology "$outname"_RAxML_bestTrees.txt | nw_order - | sort | uniq -c | sort -nr > "$outname"_RAxML_bestTrees_TOPO.txt
# BINARY
find ./TREEhackerFiles_"$outname"/ -type f -name 'RAxML_bestTree.*.concat.RN.BIN.fa.tree' -exec cat {} + > "$outname"_RAxML_bestTrees_BIN.txt
nw_topology "$outname"_RAxML_bestTrees_BIN.txt | nw_order - | sort | uniq -c | sort -nr > "$outname"_RAxML_bestTrees_TOPO_BIN.txt

