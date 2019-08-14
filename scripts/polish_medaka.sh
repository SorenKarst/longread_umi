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
CHUNK_SIZE=$2 # Sensible chunk size for amplicon type
BINNING_DIR=$3 #Raw read bins file path name
OUT_DIR=$4 #output folder name
THREADS=$5 #Number of threads
SAMPLE=$6 # List of bins to process
TURBO=${7:-NO} # Use all cores 

### Source commands and subscripts -------------------------------------
. $LONGREAD_UMI_PATH/scripts/dependencies.sh # Path to dependencies script

### Medaka polishing assembly -------------------------------------------------

# Format names
CONSENSUS_NAME=${CONSENSUS_FILE##*/}
CONSENSUS_NAME=${CONSENSUS_NAME%.*}

# Turbo mode
if [ "$TURBO" == "YES" ]; then
  CON_THREADS=$THREADS
  CON_NICE=10
elif [ "$TURBO" == "NO" ]; then
  CON_THREADS=1
  CON_NICE=0
fi

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
    --progress  \
    -j $(( THREADS * 5 )) \
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
  local CHUNK_SIZE=$5

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
    --threads $CON_THREADS \
    --model $MODEL \
	--chunk_len $CHUNK_SIZE
}

export -f consensus_wrapper

find $OUT_DIR/mapping/ \
  -mindepth 2 \
  -maxdepth 2 \
  -type f \
  -name "umi*bins.bam" |\
$GNUPARALLEL \
  --progress \
  --nice $CON_NICE \
  -j $THREADS \
  -N1 \
  --roundrobin \
  --pipe \
  "consensus_wrapper \
     {#} \
     $OUT_DIR/consensus \
     $MEDAKA_MODEL \
     $CON_THREADS \
	 $CHUNK_SIZE
  "

# Stitch consensus sequences
medaka stitch \
  $OUT_DIR/consensus/*_consensus.hdf \
  $OUT_DIR/${CONSENSUS_NAME}_medaka.fa

# Clean consensus header
sed -i "s/:.*//" $OUT_DIR/${CONSENSUS_NAME}_medaka.fa   

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
