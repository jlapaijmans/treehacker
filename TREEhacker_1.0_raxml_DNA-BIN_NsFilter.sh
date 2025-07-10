#!/bin/bash

### This script extracts inf    ./TREEhacker_1.0_raxml_DNA-BIN_NsFilter.sh [OPTIONS] <fastafiles> <o> <winsize> <stepsize> <prop_missing> [analysis_type]rmative positions from aligned genomes (multi-fasta) and computes a tree for each window
# Hacked together by Johanna Paijmans
# Modified to run RAxML calculations in parallel for faster processing


### USAGE: ./TREEhacker_1.0_raxml_DNA-BIN_NsFilter.sh $1 $2 $3 $4 $5 $6
# $1	fastafiles.txt 	  full fasta filenames)			e.g.--- ls *.fasta > fastafiles.txt. Outgroup needs to be last individual
# $2	output 		        output base name
# $3    winsize
# $4    step size (for non overlapping sliding windows, this is the same as winsize)
# $5	max prop missing data
# $6    analysis type: DNA, BIN, or BOTH (default: BOTH if not specified)

# Note: Modify the 'max_parallel_jobs' variable below to control the number of parallel RAxML jobs.
# The default is 4 parallel jobs. Increase this based on your system's CPU cores and memory.
# 
# Logging: The script creates one log file per parallel thread (up to max_parallel_jobs files)
# which are combined at the end. This approach is efficient for large-scale analyses.

### Help function
show_help() {
    cat << 'EOF'
╔════════════════════════════════════════════════════════════════════════════════╗
║                              TREEhacker v1.0                                   ║
║                    Phylogenetic Analysis in Sliding Windows                    ║
╚════════════════════════════════════════════════════════════════════════════════╝

DESCRIPTION:
    This script extracts informative positions from aligned genomes (multi-fasta) 
    and computes phylogenetic trees for each sliding window using RAxML-NG.
    
    Features:
    • Parallel processing for faster analysis
    • Sliding window approach for phylogenetic analysis
    • Support for both DNA (GTR+G) and Binary (BIN+G) models
    • Missing data filtering
    • Optimized for large-scale genomic datasets

USAGE:
    ./TREEhacker_1.0_raxml_DNA-BIN_NsFilter.sh [OPTIONS] <fastanames> <output> <winsize> <stepsize> <prop_missing> [analysis_type]

PARAMETERS:
    fastafiles      File containing full filenames of fasta files (one per line)
                    Example: ls *.fasta > fastafiles.txt
                    Note: Outgroup must be the last entry in this file
    
    output          Output base name for all generated files
    
    winsize         Size of sliding windows (in base pairs)
    
    stepsize        Step size for sliding windows (in base pairs)
                    For non-overlapping windows: stepsize = winsize
    
    prop_missing    Maximum proportion of missing data allowed per sample
                    Range: 0.0 to 1.0 (e.g., 0.5 = 50% missing data allowed)
    
    analysis_type   Type of phylogenetic analysis [optional]
                    DNA  = DNA analysis only (GTR+G model)
                    BIN  = Binary analysis only (BIN+G model)  
                    BOTH = Both analyses (default if not specified)

OPTIONS:
    --help, -h      Show this help message and exit

EXAMPLES:
    # Basic usage with both DNA and binary analysis in 1000bp non-overlapping windows, allowing 30% missing data:
    ./TREEhacker_1.0_raxml_DNA-BIN_NsFilter.sh fastafiles.txt results 1000 1000 0.3 BOTH

REQUIRED SOFTWARE:
    • samtools       - For sequence indexing and extraction
    • bedtools       - For generating sliding windows
    • bioawk         - For FASTA processing
    • seqtk          - For sequence manipulation
    • raxml-ng       - For phylogenetic inference
    • bc             - For mathematical calculations
    • nw_utils       - For tree topology analysis

PERFORMANCE TUNING:
    Edit the script to adjust these variables:
    • max_parallel_jobs = 4    (number of parallel RAxML jobs)
    • threads = 2              (threads per RAxML job)
    
OUTPUT FILES:
    • [output]_filtered_wins.txt           - Windows passing missing data filter
    • [output]_missing_data_report.txt     - Missing data statistics
    • [output]_RaXML_run.log              - DNA analysis log (if DNA/BOTH)
    • [output]_RaXML_BIN-run.log          - Binary analysis log (if BIN/BOTH)
    • [output]_raxml_trees.txt            - DNA tree collection (if DNA/BOTH)
    • [output]_raxml_trees_BIN.txt        - Binary tree collection (if BIN/BOTH)
    • [output]_raxml_trees_TOPO.txt       - DNA topology summary
    • [output]_raxml_trees_TOPO_BIN.txt   - Binary topology summary 
    • output_TREEhackerFiles_[output]/    - Directory with individual tree files

AUTHORS:
    Original: Johanna Paijmans
    Parallel optimization and enhancements: GitHub Copilot

VERSION: 1.0
EOF
}

# Check for help flags or insufficient parameters
if [[ "$1" == "--help" || "$1" == "-h" ]]; then
    show_help
    exit 0
elif [[ $# -lt 5 ]]; then
    echo "   ERROR: Insufficient parameters provided."
    echo "   Required: fastafiles output winsize stepsize prop_missing [analysis_type]"
    echo "   Provided: $# parameters"
    echo ""
    echo "Run './TREEhacker_1.0_raxml_DNA-BIN_NsFilter.sh --help' for detailed usage information."
    exit 1
fi

### variables
threads=2                 # number of threads RAXML uses
max_parallel_jobs=4       # maximum number of parallel RAxML jobs

### tools
samtools=samtools     		# path to samtools executable of choice
snpsites=snp-sites    		# path to snp-sites executable
iqtree=iqtree
bedtools=bedtools
seqtk=seqtk
bioawk=bioawk
raxml=raxml-ng

### script variables for readability, do not change
fastafiles="$1"
outname="$2"
win_size="$3"            		# window size
step_size="$4"           		# step size (equal window size for non-overlapping windows)
prop_missing="$5"				# max proportion missing data per sample in a window
analysis_type="${6:-BOTH}"      # analysis type: DNA, BIN, or BOTH (default: BOTH)

# Validate analysis type parameter
analysis_type=$(echo "$analysis_type" | tr '[:lower:]' '[:upper:]')  # Convert to uppercase
if [[ "$analysis_type" != "DNA" && "$analysis_type" != "BIN" && "$analysis_type" != "BOTH" ]]; then
    echo "Error: Invalid analysis type '$6'. Must be DNA, BIN, or BOTH."
    echo "Usage: $0 fastafiles.txt output winsize stepsize prop_missing [DNA|BIN|BOTH]"
    exit 1
fi

echo "Analysis type: $analysis_type"

# use the one below if you prefer humanised numbering "10mb" over "10000000"
#win_size=$(numfmt --to=si --suffix=b --format=%.0f $3)
#step_size=$(numfmt --to=si --suffix=b --format=%.0f $4)

max_ns=$(echo "$3*$5" | bc | xargs printf "%.0f\n")

# Cache outgroup for efficiency (avoids repeated tail calls)
outgroup=$(tail -n1 "$fastafiles" | xargs basename -s .fasta)

### Function to wait for available job slot
wait_for_jobs() {
    while [ $(jobs -r | wc -l) -ge $max_parallel_jobs ]; do
        sleep 1
    done
}

### Function to run RAxML for a single window
run_raxml_for_window() {
    local filtered_win="$1"
    local outname="$2"
    local fastafiles="$3"
    local raxml="$4"
    local threads="$5"
    local bioawk="$6"
    local win_size="$7"
    local step_size="$8"
    local samtools="$9"
    local thread_id="${10}"
    local outgroup="${11}"
    local analysis_type="${12}"
    
    echo "Processing window: $filtered_win"
    
    #3. concatenate individuals together (1 file per scaffold)
    for ind in $(cat $fastafiles);
    do
        ind_name=$(basename "$ind" .fasta)
        $samtools faidx "$ind_name"_"$win_size"_"$step_size"s_wins.fasta $filtered_win >> output_TREEhackerFiles_"$outname"/$filtered_win.concat.fasta;
    done;

    #4. rename sequences - optimized to avoid subshell
    # Create temporary file with basenames for sequence naming
    temp_basenames=$(mktemp)
    temp_sequences=$(mktemp)
    
    # Extract basenames
    while read ind; do
        basename "$ind" .fasta
    done < $fastafiles > "$temp_basenames"
    
    # Extract sequences
    $bioawk -c fastx '{print $seq}' output_TREEhackerFiles_"$outname"/"$filtered_win".concat.fasta > "$temp_sequences"
    
    # Combine with proper error checking
    if [ $(wc -l < "$temp_basenames") -eq $(wc -l < "$temp_sequences") ]; then
        paste "$temp_basenames" "$temp_sequences" | awk '{print ">"$1"\n"$2}' > output_TREEhackerFiles_"$outname"/"$filtered_win".concat.RN.fasta
    else
        echo "Error: Mismatch between number of samples and sequences in window $filtered_win" >&2
        # Fallback: use original sequence headers
        cp output_TREEhackerFiles_"$outname"/"$filtered_win".concat.fasta output_TREEhackerFiles_"$outname"/"$filtered_win".concat.RN.fasta
    fi
    
    rm -f "$temp_basenames" "$temp_sequences"
 
    ####################################################################################
    ############ NON BINARY ############
    if [[ "$analysis_type" == "DNA" || "$analysis_type" == "BOTH" ]]; then
        #6. RAXML DNA
        [[ -f output_TREEhackerFiles_"$outname"/"$filtered_win".concat.RN.fasta ]] && $raxml --threads $threads --model GTR+G --outgroup "$outgroup" --seed 12345 --msa output_TREEhackerFiles_"$outname"/"$filtered_win".concat.RN.fasta --prefix "$(pwd)"/output_TREEhackerFiles_"$outname"/"$filtered_win".concat.RN.fasta >> "$(pwd)"/"$outname"_RaXML_run_thread"$thread_id".log 2>&1;
    fi

    ####################################################################################
    ############ BINARY ############
    if [[ "$analysis_type" == "BIN" || "$analysis_type" == "BOTH" ]]; then
        #4. convert to BINARY
        [[ -f output_TREEhackerFiles_"$outname"/"$filtered_win".concat.RN.fasta ]] && $bioawk -c fastx '(gsub(/[AGag]/,"0",$seq)) && (gsub(/[CTct]/,"1",$seq)) { print ">"$name; print $seq }' output_TREEhackerFiles_"$outname"/"$filtered_win".concat.RN.fasta | $bioawk -c fastx '$seq ~ /N/ (gsub(/[Nn]/,"?",$seq)) { print ">"$name; print $seq }' > output_TREEhackerFiles_"$outname"/"$filtered_win".concat.RN.BIN.fasta

        #6. RAXML BIN
        [[ -f output_TREEhackerFiles_"$outname"/"$filtered_win".concat.RN.BIN.fasta ]] && $raxml --threads $threads --model BIN+G --outgroup "$outgroup" --seed 12345 --msa output_TREEhackerFiles_"$outname"/"$filtered_win".concat.RN.BIN.fasta --prefix "$(pwd)"/output_TREEhackerFiles_"$outname"/"$filtered_win".concat.RN.BIN.fasta >> "$(pwd)"/"$outname"_RaXML_BIN-run_thread"$thread_id".log 2>&1;
    fi

    # mv some files to stay within quota
    #rm -f output_TREEhackerFiles_"$outname"/"$filtered_win".*.fasta
    rm -f output_TREEhackerFiles_"$outname"/"$filtered_win".*.log
    rm -f output_TREEhackerFiles_"$outname"/"$filtered_win".*.mlTrees
    rm -f output_TREEhackerFiles_"$outname"/"$filtered_win".*.bestModel
    rm -f output_TREEhackerFiles_"$outname"/"$filtered_win".*.rba
    rm -f output_TREEhackerFiles_"$outname"/"$filtered_win".*.phy
    rm -f output_TREEhackerFiles_"$outname"/"$filtered_win".*.startTree
    rm -f output_TREEhackerFiles_"$outname"/"$filtered_win".*.bestTreeCollapsed

    
    echo "Completed window: $filtered_win"
}

### preamble:
# make necessary directory
[[ -d output_TREEhackerFiles_"$outname"/ ]] && { echo "Output folder output_TREEhackerFiles_"$outname"/ exists, remove or select different folder"; exit 1; }
mkdir output_TREEhackerFiles_"$outname"/

### the script
##1. generate windows from each fasta using bedtools, and then removing the small scaffolds (tail ends)
echo "Generating sliding windows for each individual..."
for ind in $(cat $fastafiles);
do
    # Extract basename without extension (e.g., fasta1.fasta -> fasta1)
    ind_name=$(basename "$ind" .fasta)
    if [[ ! -f "$ind_name"_"$win_size"_"$step_size"s_wins.fasta ]]; then
        echo "Processing $ind..."
        $bioawk -c fastx '{print $name"\t0\t"length($seq)}' "$ind" | $bedtools makewindows -b - -w $win_size -s $step_size | $bedtools getfasta -fi "$ind" -bed - -fo stdout | $seqtk seq -L $win_size - | sed 's/:/_/g'> "$ind_name"_"$win_size"_"$step_size"s_wins.fasta
        $samtools faidx "$ind_name"_"$win_size"_"$step_size"s_wins.fasta
    else
        echo "$ind_name windows already exist, skipping..."
    fi;
done;

### generate a report of the number of N's in each scaffold (optimized)
if [[ ! -f "$outname"_missing_data_report.txt ]]; then
    echo "Generating missing data report..."
    # Get all unique windows first
    cat *_"$win_size"_"$step_size"s_wins.fasta.fai | awk '{print $1}' | sort -V | uniq > temp_windows_list.txt
    
    # Process all windows at once to reduce I/O
    while read win; do
        printf "%s" "$win"
        for ind in $(cat $fastafiles); do
            ind_name=$(basename "$ind" .fasta)
            n_count=$($samtools faidx "$ind_name"_"$win_size"_"$step_size"s_wins.fasta "$win" | $seqtk comp | awk '{print $9}')
            printf "\t%s" "$n_count"
        done
        printf "\n"
    done < temp_windows_list.txt > "$outname"_missing_data_report.txt
    
    rm -f temp_windows_list.txt
fi;

### filter for windows with more than "$max_ns" N's
[[ ! -f "$outname"_filtered_wins.txt ]] && awk -v max_ns=$max_ns '{ok=1;for (x=2;x<=NF;x++) if ($x>max_ns) ok=0; if (ok) print $0}' "$outname"_missing_data_report.txt > "$outname"_filtered_wins.txt


#2. generate alignments per window and run RAxML in parallel
total_windows=$(wc -l < "$outname"_filtered_wins.txt)
echo "Starting parallel processing of $total_windows windows with max $max_parallel_jobs parallel jobs..."

# Initialize thread counter and progress tracking
thread_counter=0
completed_windows=0

for filtered_win in $(cat "$outname"_filtered_wins.txt | awk '{print $1}');
do
    # Wait for available job slot
    wait_for_jobs
    
    # Assign thread ID (cycling through available threads)
    current_thread_id=$((thread_counter % max_parallel_jobs))
    thread_counter=$((thread_counter + 1))
    
    # Run RAxML for this window in background
    run_raxml_for_window "$filtered_win" "$outname" "$fastafiles" "$raxml" "$threads" "$bioawk" "$win_size" "$step_size" "$samtools" "$current_thread_id" "$outgroup" "$analysis_type" &
    
    # Progress reporting every 100 windows
    if (( thread_counter % 100 == 0 )); then
        echo "Submitted $thread_counter/$total_windows windows for processing..."
    fi
done

# Wait for all background jobs to complete
echo "Waiting for all RAxML jobs to complete..."
wait
echo "All RAxML jobs completed!"

# Combine thread-based log files into single files
echo "Combining log files..."
if [[ "$analysis_type" == "DNA" || "$analysis_type" == "BOTH" ]]; then
    cat "$outname"_RaXML_run_thread*.log > "$outname"_RaXML_run.log 2>/dev/null
fi
if [[ "$analysis_type" == "BIN" || "$analysis_type" == "BOTH" ]]; then
    cat "$outname"_RaXML_BIN-run_thread*.log > "$outname"_RaXML_BIN-run.log 2>/dev/null
fi

# Optionally remove individual thread log files (uncomment if you want to keep only the combined logs)
rm -f "$outname"_RaXML_run_thread*.log
rm -f "$outname"_RaXML_BIN-run_thread*.log


####################################################################################
############ SUMMARISING ############

# NON-BINARY
if [[ "$analysis_type" == "DNA" || "$analysis_type" == "BOTH" ]]; then
    echo "Summarizing DNA trees..."
    find ./output_TREEhackerFiles_"$outname"/ -type f -name '*.concat.RN.fasta.raxml.bestTree' -exec cat {} + > "$outname"_raxml_trees.txt 2>/dev/null
    if command -v nw_topology &> /dev/null && command -v nw_order &> /dev/null; then
        nw_topology "$outname"_raxml_trees.txt | nw_order - | sort | uniq -c > "$outname"_raxml_trees_TOPO.txt
    else
        echo "Warning: nw_utils not available, skipping topology analysis for DNA trees"
    fi
fi

# BINARY
if [[ "$analysis_type" == "BIN" || "$analysis_type" == "BOTH" ]]; then
    echo "Summarizing BIN trees..."
    find ./output_TREEhackerFiles_"$outname"/ -type f -name '*.concat.RN.BIN.fasta.raxml.bestTree' -exec cat {} + > "$outname"_raxml_trees_BIN.txt 2>/dev/null
    if command -v nw_topology &> /dev/null && command -v nw_order &> /dev/null; then
        nw_topology "$outname"_raxml_trees_BIN.txt | nw_order - | sort | uniq -c > "$outname"_raxml_trees_TOPO_BIN.txt
    else
        echo "Warning: nw_utils not available, skipping topology analysis for BIN trees"
    fi
fi


