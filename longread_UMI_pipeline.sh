#!/bin/bash
# DESCRIPTION
#    longread-UMI-pipeline script. 
#    
# IMPLEMENTATION
#    author	Søren Karst (sorenkarst@gmail.com)
#               Ryans Ziels (ziels@mail.ubc.ca)
#    version	0.1.0
#    license	GNU General Public License
#
# To-do:
# - Fix logging

### Source commands and subscripts -------------------------------------
. scripts/dependencies.sh # Path to dependencies script
THREADS=60 # Number of threads to use. Ps. Medaka bug means medaka uses all available..
UMI_SUBSET_N=1000000 # Number of UMI bins to use for post processing i.e. for testing subset to 100
UMI_COVERAGE_MIN=30 # UMI bin coverage cut off for variant phasing.

### Pipeline -----------------------------------------------------------
# Logging
LOG_NAME="ncec_log_$(date +"%Y-%m-%d-%T").txt"
echo "ncec log" >> $LOG_NAME
echo "Script start: $(date)" >> $LOG_NAME
ncec_version_dump $LOG_NAME
exec &> >(tee -a "$LOG_NAME")
exec 2>&1
echo ""
echo "### Settings:"
echo "Threads: $THREADS"
echo "UMI subsampling: $UMI_SUBSET_N"
echo "References: $REF"
echo "References unique: $REFU"
echo "Bin size cutoff: $UMI_COVERAGE_MIN"
echo ""

# Fetch data
# cat ./path/to/fastq/files/*fastq > ./reads.fq
cat ../data/ncec-7-guppy3.0.3/*.fastq > reads.fq

# Read filtering and UMI binning
$UMI_BINNING \
  # Raw nanopore data in fastq format
  reads.fq \
  # Output folder
  umi_binning \
  # Number of threads
  $THREADS \
  # Min read length
  3500 \
  # Max read length
  6000

# Sample UMI bins for testing
find umi_binning/read_binning/bins \
   -name 'umi*bins.fastq' | sed -e 's|^.*/||' -e 's|\..*||' |\
   head -n $UMI_SUBSET_N> sample$UMI_SUBSET_N.txt

# Consensus
$CONSENSUS_SRACON \
  # Path to UMI bins
  umi_binning/read_binning/bins \
  # Output folder
  sracon \
  # Number of threads
  $THREADS \
  # List of bins to process
  sample$UMI_SUBSET_N.txt

# Polishing
$POLISH_MEDAKA \
  # Path to consensus data
  sracon/consensus_*.fa \
  # Path to UMI bins
  umi_binning/read_binning/bins \
  # Output folder
  sracon_medaka \
  # Number of threads
  $THREADS \
  # List of bins to process
  sample$UMI_SUBSET_N.txt

$POLISH_MEDAKA \
  # Path to consensus data
  sracon_medaka/consensus_*.fa \
  # Path to UMI bins
  umi_binning/read_binning/bins \
  # Output folder
  sracon_medaka_medaka \
  # Number of threads
  $THREADS \
  # List of bins to process
  sample$UMI_SUBSET_N.txt

# Trim UMI consensus data
$TRIM_AMPLICON \
  # Path to consensus data
  sracon_medaka_medaka\
  # Consensus file pattern
  consensus*fa \
  # Output folder
  . \
  # Primers used
  rrna_8f2490r \
  $THREADS

# Generate phasing data
mkdir phasing

## Subset to UMI consensus sequences with min read coverage
awk -v bsco="$UMI_COVERAGE_MIN" '
  /^>/{
    s=$0
    gsub(".*size=", "", s)
    if (s+0 >= bsco+0){
      print
      getline
      print
    }
  }
' consensus_sracon_medaka_medaka.fa \
> consensus_sracon_medaka_medaka_${UMI_COVERAGE_MIN}.fa

## Variant calling of from UMI consensus sequences
$VARIANTS \
  # Path to consensus data
  consensus_sracon_medaka_medaka_${UMI_COVERAGE_MIN}.fa \
  # Output folder
  variants \
  # Number of threads
  $THREADS

## Copy phase
cp variants/variants*fa .


