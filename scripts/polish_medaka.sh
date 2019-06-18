#!/bin/bash
# DESCRIPTION
#    Consensus sequence polishing using medaka.
#    This script is a part of the longread-UMI-pipeline.
#    
# IMPLEMENTATION
#    author	Søren Karst (sorenkarst@gmail.com)
#               Ryans Ziels (ziels@mail.ubc.ca)
#    license	GNU General Public License
# TODO
#

### Terminal input ------------------------------------------------------------
CONSENSUS_FILE=$1 #Consensus file path name
BINNING_DIR=$2 #Raw read bins file path name
OUT_DIR=$3 #output folder name
THREADS=$4 #Number of threads
SAMPLE=$5 # List of bins to process

### Medaka polishing assembly -------------------------------------------------

# Format names
CONSENSUS_NAME=${CONSENSUS_FILE##*/}
CONSENSUS_NAME=${CONSENSUS_NAME%.*}

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
  local UMI_BIN=$BINNING_DIR/*/${UMI_NAME}.fa*
  local RC=$(( $(wc -l < $UMI_BIN)/4 ))  

  # Setup working directory
  mkdir $OUT_DIR/$UMI_NAME
  echo "$IN" |\
    awk -v rc="$RC" '
      /^>/{print $0 ":size=" rc}
      !/^>/
    ' > $OUT_DIR/$UMI_NAME/$UMI_NAME.fa

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
  ( [[ "${SAMPLE}" ]] && grep -A1 -Ff $SAMPLE | sed '/^--$/d' || cat ) |\
  $GNUPARALLEL \
    --progress  \
    -j $(( THREADS * 10 )) \
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
  local MODEL=$3
  local CON_THREADS=$4

  # Merge bams
  $SAMTOOLS merge \
    -b <(cat) \
    $OUT_DIR/${JOB_NR}.bam

  # Index bam
  $SAMTOOLS index \
    $OUT_DIR/${JOB_NR}.bam

  # Medaka consensus
  medaka consensus \
    $OUT_DIR/${JOB_NR}.bam \
    $OUT_DIR/${JOB_NR}_consensus.hdf \
    --threads 1 \
    --model $MODEL
}

export -f consensus_wrapper

find $OUT_DIR/mapping/ \
  -mindepth 2 \
  -maxdepth 2 \
  -type f \
  -name "umi*bins.bam" |\
$GNUPARALLEL \
  --progress \
  -j $THREADS \
  -N1 \
  --roundrobin \
  --pipe \
  "consensus_wrapper \
     {#} \
     $OUT_DIR/consensus \
     $MEDAKA_MODEL
  "

# Stitch consensus sequences
medaka stitch \
  $OUT_DIR/consensus/*_consensus.hdf \
  $OUT_DIR/${CONSENSUS_NAME}_medaka.fa 

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
