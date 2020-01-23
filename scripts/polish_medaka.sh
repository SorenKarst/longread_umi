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

### Terminal input ------------------------------------------------------------
CONSENSUS_FILE=$1 #Consensus file path name
MEDAKA_MODEL=$2 #Medaka model to use.
CHUNK_SIZE=$3 # Sensible chunk size for amplicon type
BINNING_DIR=$4 #Raw read bins file path name
OUT_DIR=$5 #output folder name
THREADS=$6 #Number of threads
SAMPLE=$7 # List of bins to process
MEDAKA_JOBS=${8:-1} # Number medaka jobs to start

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
