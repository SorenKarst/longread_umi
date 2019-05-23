#!/bin/bash
# DESCRIPTION
#    Consensus sequence polishing using medaka.
#    This script is a part of the longread-UMI-pipeline.
#    
# IMPLEMENTATION
#    author	Søren Karst (sorenkarst@gmail.com)
#               Ryans Ziels (ziels@mail.ubc.ca)
#    version	0.1.0
#    license	GNU General Public License
# TODO
#

### Terminal input ------------------------------------------------------------
CONSENSUS_FILE=$1 #Consensus file path name
BINNING_DIR=$2 #Raw read bins file path name
OUT_DIR=$3 #output folder name
THREADS=$4 #Number of threads
SAMPLE=$5 # List of bins to process
RV=${6:-no} # Reverse complement reference?

### Medaka polishing assembly -------------------------------------------------

# Format names
CONSENSUS_NAME=${CONSENSUS_FILE##*/}
CONSENSUS_NAME=${CONSENSUS_NAME%.*}
[[ $RV == "yes" ]] && RV_PREFIX="rc"

# Start medaka environment if relevant
eval "$MEDAKA_ENV_START"

# Prepare output folders
if [ -d "$OUT_DIR" ]; then
  echo "Output folder exists. Exiting..."
  exit 0
fi
mkdir $OUT_DIR

# Custom function
medaka_polish() {

  # Input
  local IN=$(cat)
  local BINNING_DIR=$1
  local OUT_DIR=$2
  local MEDAKA_MODEL=$3
  local RV=$4

  # Name format
  local UMI_NAME=$(echo "$IN" | grep -o "umi.*bins")
  local UMI_BIN=$(find $BINNING_DIR/ -name ${UMI_NAME}.fa*)
  local RC=$(( $(wc -l < $UMI_BIN)/4 )) 

  # Prepare output folder
  local UMI_OUT=$OUT_DIR/$UMI_NAME
  mkdir $UMI_OUT
  echo "$IN" |\
    awk -v umi_name="$UMI_NAME" '
      /^>/{print ">" umi_name}
      !/^>/
    ' > $UMI_OUT/$UMI_NAME.fa
  if [ $RV == "yes" ]; then
    RVSEQ=$($SEQTK seq -r $UMI_OUT/$UMI_NAME.fa)
    echo "$RVSEQ" > $UMI_OUT/$UMI_NAME.fa
  fi

  # Perform polishing
  medaka_consensus -i $UMI_BIN \
   -d $UMI_OUT/$UMI_NAME.fa \
   -o $UMI_OUT \
   -m $MEDAKA_MODEL \
   -t 1

  if [ -f $UMI_OUT/consensus.fasta ]; then
    awk -v rc="$RC" '
      /^>/{print $0 ";size=" rc}
      !/^>/
    ' $UMI_OUT/consensus.fasta > $UMI_OUT/${UMI_NAME}_md.fa
    rm $UMI_OUT/consensus.fasta
  fi
}

export -f medaka_polish

# Perform medaka polishing in parallel
cat $CONSENSUS_FILE | $SEQTK seq -l0 - |\
  ( [[ "${SAMPLE}" ]] && grep -A1 -Ff $SAMPLE | sed '/^--$/d' || cat ) |\
  $GNUPARALLEL --progress -j$THREADS --recstart ">" -N 1 --pipe \
  "cat | medaka_polish $BINNING_DIR $OUT_DIR $MEDAKA_MODEL $RV"

find $OUT_DIR -name '*_md.fa' -exec cat {} \; \
  > $OUT_DIR/${CONSENSUS_NAME}_medaka${RV_PREFIX}.fa

# Deactivate medaka environment if relevant
eval "$MEDAKA_ENV_STOP"

