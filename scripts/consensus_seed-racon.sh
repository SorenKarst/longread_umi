#!/bin/bash
# DESCRIPTION
#    Script for finding consensus from UMI binned Nanopore reads. 
#    Raw read centroid found with usearch and used as seed for
#    4 x times racon polishing.
#
#    This script is a part of the longread-UMI-pipeline.
#
# IMPLEMENTATION
#    author	Søren Karst (sorenkarst@gmail.com)
#               Ryans Ziels (ziels@mail.ubc.ca)
#    version	0.1.0
#    license	GNU General Public License

### Terminal input ------------------------------------------------------------
IN=$1 #input folder name with bins (IN=read_binning/bins/)
OUT=$2 #output folder name
THREADS=$3 #number of threads
SAMPLE=$4 # List of bins to process

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

  # Name format
  local UMINO=${RB##*/}
  local UMINO=${UMINO%.*}
  local OUT=$OUT/$UMINO

  # Create dirs    
  mkdir $OUT

  # Find seed read
  $USEARCH -cluster_fast $RB -id 0.75 -strand both\
    -sizeout -centroids $OUT/${UMINO}_centroids.fa
  $USEARCH -sortbysize $OUT/${UMINO}_centroids.fa -topn 1 \
    -fastaout $OUT/${UMINO}_sr.fa -relabel seed

  # Racon polishing
  for i in `seq 1 5`; do
    $MINIMAP2 \
      -t 1 \
      -x ava-ont \
      $OUT/${UMINO}_sr.fa \
      $RB > $OUT/ovlp.paf

    $RACON \
      -t 1 \
      -m 8 \
      -x -6 \
      -g -8 \
      -w 500 \
      $RB \
      $OUT/ovlp.paf \
      $OUT/${UMINO}_sr.fa > $OUT/${UMINO}_tmp.fa
    mv $OUT/${UMINO}_tmp.fa $OUT/${UMINO}_sr.fa
  done
  
  # Rename
  sed -i "s/^>.*/>$UMINO/" $OUT/${UMINO}_sr.fa                    
} 

export -f seed_racon

# Perform assembly in parallel
find $IN -name 'umi*bins.fastq'  |\
  ( [[ "${SAMPLE}" ]] && grep -Ff $SAMPLE || cat ) |\
  $GNUPARALLEL --progress -j $THREADS "seed_racon {} $OUT"

#Collect seed-racon consensus sequences
find $OUT/ \
  -mindepth 2 \
  -maxdepth 2 \
  -name "*_sr*.fa" \
  -exec cat {} \; \
  > $OUT/consensus_${OUT}.fa
