#!/bin/bash
# DESCRIPTION
#    Script for finding consensus from UMI read bins. 
#    Raw read centroid found with usearch and used as seed for
#    <ROUNDS> x times racon polishing.
#
#    This script is a part of the longread-UMI-pipeline.
#
# IMPLEMENTATION
#    author	Søren Karst (sorenkarst@gmail.com)
#               Ryan Ziels (ziels@mail.ubc.ca)
#    license	GNU General Public License

### Description ----------------------------------------------------------------

USAGE="
-- longread_umi consensus_racon: Generate UMI consensus sequence with racon
   Raw read centroid found with usearch and used as seed for
   (r) x times racon polishing.

usage: $(basename "$0" .sh) [-h] (-d dir -o dir -p string -r value -t value -n file -a string) 

where:
    -h  Show this help text.
    -d  Directory containing UMI read bins in the format
        'umi*bins.fastq'. Recursive search.
    -o  Output directory.
    -p  Minimap2 preset. 'map-ont' for Nanopore and 'asm20' for PacBio CCS.
    -a  Additional arguments for racon.
    -r  Number of racon polishing rounds.
    -t  Number of threads to use.
    -n  Process n number of bins. If not defined all bins
        are processed.
"

### Terminal Arguments ---------------------------------------------------------

# Import user arguments
while getopts ':hzd:o:p:a:r:t:n:' OPTION; do
  case $OPTION in
    h) echo "$USAGE"; exit 1;;
    d) IN=$OPTARG;;
    o) OUT=$OPTARG;;
    p) PRESET=$OPTARG;;
    a) RACON_ARG=$OPTARG;;
    r) ROUNDS=$OPTARG;;
    t) THREADS=$OPTARG;;
    n) SAMPLE=$OPTARG;;
    :) printf "missing argument for -$OPTARG\n" >&2; exit 1;;
    \?) printf "invalid option for -$OPTARG\n" >&2; exit 1;;
  esac
done

# Check missing arguments
MISSING="is missing but required. Exiting."
if [ -z ${IN+x} ]; then echo "-d $MISSING"; echo "$USAGE"; exit 1; fi; 
if [ -z ${OUT+x} ]; then echo "-o $MISSING"; echo "$USAGE"; exit 1; fi; 
if [ -z ${PRESET+x} ]; then echo "-p $MISSING"; echo "$USAGE"; exit 1; fi;
if [ -z ${RACON_ARG+x} ]; then RACON_ARG=""; fi; 
if [ -z ${ROUNDS+x} ]; then echo "-r $MISSING"; echo "$USAGE"; exit 1; fi; 
if [ -z ${THREADS+x} ]; then echo "-t is missing. Defaulting to 1 thread."; THREADS=1; fi;

### Source commands and subscripts -------------------------------------
. $LONGREAD_UMI_PATH/scripts/dependencies.sh # Path to dependencies script

### Usearch centroid and racon polishing --------------------------------------

# Prepare output folders
if [ -d "$OUT" ]; then
  echo "Output folder exists. Exiting..."
  exit 0
fi
mkdir $OUT

# Wrapper
seed_racon () {
  # Input
  local RB=$1
  local OUT=$2
  local ROUNDS=$3
  local PRESET=$4
  local RACON_ARG=$5

  # Name format
  local UMINO=${RB##*/}
  local UMINO=${UMINO%.*}
  local OUT=$OUT/$UMINO

  # Create dirs    
  mkdir $OUT

  # Count binsize
  local BINSIZE=$($GAWK 'NR%4==1 {N++} END {print N}' $RB)

  # Find seed read
  $USEARCH \
    -cluster_fast \
    $RB \
    -id 0.75 \
    -strand both\
    -sizeout \
    -centroids $OUT/${UMINO}_centroids.fa

  $USEARCH \
    -sortbysize $OUT/${UMINO}_centroids.fa \
    -topn 1 \
    -fastaout $OUT/${UMINO}_sr.fa \
    -relabel seed

  # Racon polishing
  for i in `seq 1 $ROUNDS`; do
    $MINIMAP2 \
      -t 1 \
      -x $PRESET \
      $OUT/${UMINO}_sr.fa \
      $RB > $OUT/ovlp.paf

    $RACON \
      -t 1 \
      -m 8 \
      -x -6 \
      -g -8 \
      -w 25000 \
      $RACON_ARG \
      $RB \
      $OUT/ovlp.paf \
      $OUT/${UMINO}_sr.fa > $OUT/${UMINO}_tmp.fa
    mv $OUT/${UMINO}_tmp.fa $OUT/${UMINO}_sr.fa
  done
  
  # Rename and add binsize
  sed -i "s/^>.*/>$UMINO;ubs=$BINSIZE/" $OUT/${UMINO}_sr.fa                    
} 

export -f seed_racon

# Perform assembly in parallel
find $IN -name 'umi*bins.fastq'  |\
  ( [[ -f "${SAMPLE}" ]] && grep -Ff $SAMPLE || cat ) |\
  $GNUPARALLEL \
    --env seed_racon \
    --progress \
    -j $THREADS \
    "seed_racon {} $OUT $ROUNDS $PRESET $RACON_ARG"

#Collect seed-racon consensus sequences
OUT_NAME=${OUT##*/}
OUT_NAME=${OUT_NAME%.*}

find $OUT/ \
  -mindepth 2 \
  -maxdepth 2 \
  -name "*_sr*.fa" \
  -exec cat {} \; \
  > $OUT/consensus_${OUT_NAME}.fa
