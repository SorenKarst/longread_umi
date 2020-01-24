#!/bin/bash
# DESCRIPTION
#    Consensus sequence polishing using medaka.
#    This script is a part of the longread-UMI-pipeline.
#    
# IMPLEMENTATION
#    author   Søren Karst (sorenkarst@gmail.com)
#             Ryan Ziels (ziels@mail.ubc.ca)
#    license  GNU General Public License
# TODO
#

### Description ----------------------------------------------------------------

USAGE="
-- longread_umi polish_medaka: Nanopore UMI consensus polishing with Medaka
   
usage: $(basename "$0" .sh) [-h] [-l value T value] 
(-c file -m string -d dir -o dir -t value -n file -T value)

where:
    -h  Show this help text.
    -c  File containing consensus sequences.
    -m  Medaka model.
    -l  Expected minimum chunk size. [Default = 6000]
    -d  Directory containing UMI read bins in the format
        'umi*bins.fastq'. Recursive search.
    -o  Output directory.
    -t  Number of threads to use.
    -n  Process n number of bins. If not defined all bins
        are processed.
    -t  Number of Medaka jobs to run. [Default = 1].
"

### Terminal Arguments ---------------------------------------------------------

# Import user arguments
while getopts ':hzc:m:l:d:o:t:n:T:' OPTION; do
  case $OPTION in
    h) echo "$USAGE"; exit 1;;
    c) CONSENSUS_FILE=$OPTARG;;
    m) MEDAKA_MODEL=$OPTARG;;
    l) CHUNK_SIZE=$OPTARG;;
    d) BINNING_DIR=$OPTARG;;
    o) OUT_DIR=$OPTARG;;
    t) THREADS=$OPTARG;;
    n) SAMPLE=$OPTARG;;
    T) MEDAKA_JOBS=$OPTARG;;
    :) printf "missing argument for -$OPTARG\n" >&2; exit 1;;
    \?) printf "invalid option for -$OPTARG\n" >&2; exit 1;;
  esac
done

# Check missing arguments
MISSING="is missing but required. Exiting."
if [ -z ${CONSENSUS_FILE+x} ]; then echo "-c $MISSING"; echo "$USAGE"; exit 1; fi; 
if [ -z ${MEDAKA_MODEL+x} ]; then echo "-m $MISSING"; echo "$USAGE"; exit 1; fi; 
if [ -z ${CHUNK_SIZE+x} ]; then echo "-l missing. Defaulting to 6000."; CHUNK_SIZE=6000; fi;
if [ -z ${BINNING_DIR+x} ]; then echo "-d $MISSING"; echo "$USAGE"; exit 1; fi; 
if [ -z ${OUT_DIR+x} ]; then echo "-o $MISSING"; echo "$USAGE"; exit 1; fi; 
if [ -z ${THREADS+x} ]; then echo "-t $MISSING"; echo "$USAGE"; exit 1; fi; 
if [ -z ${MEDAKA_JOBS+x} ]; then echo "-T is missing. Defaulting to 1 Medaka job."; MEDAKA_JOBS=1; fi;

### Source commands and subscripts -------------------------------------
. $LONGREAD_UMI_PATH/scripts/dependencies.sh # Path to dependencies script

### Medaka polishing assembly -------------------------------------------------

# Format names
OUT_NAME=${OUT_DIR##*/}

# Medaka jobs
MEDAKA_THREADS=$(( THREADS/MEDAKA_JOBS ))

# Start medaka environment if relevant
eval "$MEDAKA_ENV_START"

# Prepare output folders
if [ -d "$OUT_DIR" ]; then
  echo "Output folder exists. Exiting..."
  exit 0
fi
mkdir $OUT_DIR

# Individual mapping of UMI bins to consensus

mkdir $OUT_DIR/mapping 

medaka_align() {
  # Input
  local IN=$(cat)
  local BINNING_DIR=$1
  local OUT_DIR=$2

  # Name format
  local UMI_NAME=$(echo "$IN" | grep -o "umi.*bins")
  local UMI_BIN=$(find $BINNING_DIR -name ${UMI_NAME}.fastq)

  # Setup working directory
  mkdir $OUT_DIR/$UMI_NAME
  echo "$IN" > $OUT_DIR/$UMI_NAME/$UMI_NAME.fa

  # Map UMI reads to consensus
  mini_align \
    -i $UMI_BIN \
    -r $OUT_DIR/$UMI_NAME/$UMI_NAME.fa \
    -m \
    -p $OUT_DIR/$UMI_NAME/${UMI_NAME} \
    -t 1
}

export -f medaka_align

cat $CONSENSUS_FILE |\
  $SEQTK seq -l0 - |\
  ( [[ -f "${SAMPLE}" ]] && grep -A1 -Ff $SAMPLE | sed '/^--$/d' || cat ) |\
  $GNUPARALLEL \
    --env medaka_align \
    --progress  \
    -j $THREADS \
    --recstart ">" \
    -N 1 \
    --pipe \
    "medaka_align \
       $BINNING_DIR \
       $OUT_DIR/mapping
    "

# Calculate consensus probabilities
mkdir $OUT_DIR/consensus

consensus_wrapper() {
  # Input
  local JOB_NR=$1
  local OUT_DIR=$2
  local MEDAKA_MODEL=$3
  local MEDAKA_THREADS=$4
  local CHUNK_SIZE=$5

  # Merge bam files
  
  ### Custom merge function required for > 1024 files
  bam_merge() {
    # Input
    local OUT_DIR=$1
	local JOB_NR=$2
  
    awk -v OUT_BAM="$OUT_DIR/${JOB_NR}" '
      # Print HD line
      (NR==1 && $1 ~ /^@HD$/){print $0 > OUT_BAM ".header"}
      # Print reference lines
      ($1 ~ /^@SQ$/){print $0 > OUT_BAM ".header"}
      # Print reads lines
      ($1 ~ !/^@/){print $0 > OUT_BAM ".body"}
    '
  }  
  export -f bam_merge
  
  ### View bam files in parallel and pipe into merge function
  cat |\
    $GNUPARALLEL \
      --env bam_merge \
      -j $MEDAKA_THREADS \
      $SAMTOOLS view -h {} |\
      bam_merge $OUT_DIR $JOB_NR

  ### Build bam file
  cat $OUT_DIR/${JOB_NR}.header $OUT_DIR/${JOB_NR}.body  |\
    samtools view -b - > $OUT_DIR/${JOB_NR}.bam
  rm $OUT_DIR/${JOB_NR}.header $OUT_DIR/${JOB_NR}.body

  # Index bam
  $SAMTOOLS index \
    $OUT_DIR/${JOB_NR}.bam

  # Medaka consensus
  medaka consensus \
    $OUT_DIR/${JOB_NR}.bam \
    $OUT_DIR/${JOB_NR}_consensus.hdf \
    --threads $MEDAKA_THREADS \
    --model $MEDAKA_MODEL \
	  --chunk_len $CHUNK_SIZE
}

export -f consensus_wrapper

find $OUT_DIR/mapping/ \
  -mindepth 2 \
  -maxdepth 2 \
  -type f \
  -name "umi*bins.bam" |\
$GNUPARALLEL \
  --env consensus_wrapper \
  --progress \
  -j $MEDAKA_JOBS \
  -N1 \
  --roundrobin \
  --pipe \
  "consensus_wrapper \
     {#} \
     $OUT_DIR/consensus \
     $MEDAKA_MODEL \
     $MEDAKA_THREADS \
     $CHUNK_SIZE"

# Stitch consensus sequences
medaka stitch \
  $OUT_DIR/consensus/*_consensus.hdf \
  $OUT_DIR/consensus_${OUT_NAME}.fa

# Clean consensus header
sed -i -e "s/:.*//" -e "s/_segment.*//" $OUT_DIR/consensus_${OUT_NAME}.fa

# Deactivate medaka environment if relevant
eval "$MEDAKA_ENV_STOP"

# Testing
exit 0
/space/users/smk/Software/longread-UMI-pipeline/scripts/polish_medaka.sh \
  sracon/consensus_sracon.fa \
  umi_binning/read_binning/bins \
  test_new \
  24 \
  sample500.txt

CONSENSUS_FILE=racon/consensus_racon.fa
BINNING_DIR=umi_binning/read_binning/bins
OUT_DIR=racon_medaka
THREADS=1
SAMPLE=10
