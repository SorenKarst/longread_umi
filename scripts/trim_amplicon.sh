#!/bin/bash
# DESCRIPTION
#    Script for trimming sequences based on gene primers.
#    This script is a part of the longread-UMI-pipeline.
#    
# IMPLEMENTATION
#    author	Søren Karst (sorenkarst@gmail.com)
#               Ryan Ziels (ziels@mail.ubc.ca)
#    license	GNU General Public License

### Terminal input ------------------------------------------------------------
IN_DIR=$1
IN_REGEX=$2
OUT_DIR=$3
FW2=$4
RV2=$5
MIN_LENGTH=$6
MAX_LENGTH=$7
THREADS=$8
LOG_DIR=${9:-logs}

### Source commands and subscripts -------------------------------------
. $LONGREAD_UMI_PATH/scripts/dependencies.sh # Path to dependencies script

#Format input
IN_DIR_F=$(echo $IN_DIR | sed -e 's/[,;\t]/ /g')
IN_REGEX_F=$(echo $IN_REGEX |\
  sed -e 's/^/-name /g' \
      -e 's/[,;\t]/ -o -name /g')

### Primer formating
revcom() {
  echo $1 |\
  $GAWK '{print ">dummy\n" $0}' |\
  $SEQTK seq -r - |\
  $GAWK '!/^>/'  
}
RV2R=$(revcom "$RV2")


# Custom functions
cutadapt_wrapper(){
  # Input
  local IN=$1
  local OUT_DIR=$2
  local FW2=$3
  local RV2R=$4
  local MIN_LENGTH=$5
  local MAX_LENGTH=$6
  local LOG_DIR=$7
  local THREADS=$8

  # Format input name
  local IN_NAME=${IN##*/}
  local IN_NAME=${IN_NAME%.*}
  
  # Run cutadapt
  ### Trims all sequences with correct orientation
  $CUTADAPT \
    --untrimmed-output $OUT_DIR/${IN_NAME}.tmp \
    -m $MIN_LENGTH \
    -M $MAX_LENGTH \
    -g $FW2...$RV2R \
    $IN > $OUT_DIR/${IN_NAME}.fa 2> $LOG_DIR/${IN_NAME}_amplicon_trim_log.txt
  ### Reverse/complement untrimmed and trim again
  $SEQTK seq -r $OUT_DIR/$IN_NAME.tmp |\
    $CUTADAPT \
      --discard-untrimmed \
      -m $MIN_LENGTH \
      -M $MAX_LENGTH \
      -g $FW2...$RV2R - >> $OUT_DIR/${IN_NAME}.fa 2>> $LOG_DIR/${IN_NAME}_amplicon_trim_log.txt
  rm $OUT_DIR/$IN_NAME.tmp
}

export -f cutadapt_wrapper

# Trim sequences
mkdir -p $OUT_DIR
FIND_CMD="find $IN_DIR_F -mindepth 1 -maxdepth 1 -type f \( $IN_REGEX_F \)"
eval $FIND_CMD |\
  $GNUPARALLEL \
    --env cutadapt_wrapper \ 
    -j $THREADS \
    cutadapt_wrapper \
      {} \
      $OUT_DIR \
      $FW2 \
      $RV2R \
      $MIN_LENGTH \
      $MAX_LENGTH \
      $LOG_DIR \
      $THREADS

