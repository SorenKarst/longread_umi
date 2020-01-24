#!/bin/bash
# DESCRIPTION
#    Script for trimming sequences based on gene primers.
#    This script is a part of the longread-UMI-pipeline.
#    
# IMPLEMENTATION
#    author   Søren Karst (sorenkarst@gmail.com)
#             Ryan Ziels (ziels@mail.ubc.ca)
#    license  GNU General Public License

### Description ----------------------------------------------------------------

USAGE="
-- longread_umi trim_amplicon: Trimming sequences based on primers
   
usage: $(basename "$0" .sh) [-h] (-d dir(s) -p string(s) -o dir )
(-F string -R string -m value -M value -t value -l dir) 

where:
    -h  Show this help text.
    -d  Directory to look for sequence files using pattern (-p).
        Mutliple directories can be seperated by ';'.
    -p  File pattern(s) to look for. Multiple patterns
        can be separared by ';' and flanked by '\"...\"'.
    -o  Output directory.
    -F  Forward primer sequence.
    -R  Reverse primer sequence.
    -m  Minimum read length.
    -M  Maximum read length.
    -t  Number of threads to use.
    -l  Log directory

"

### Terminal Arguments ---------------------------------------------------------

# Import user arguments
while getopts ':hzd:p:o:F:R:m:M:t:l:' OPTION; do
  case $OPTION in
    h) echo "$USAGE"; exit 1;;
    d) IN_DIR=$OPTARG;;
    p) IN_REGEX=$OPTARG;;
    o) OUT_DIR=$OPTARG;;
    F) FW2=$OPTARG;;
    R) RV2=$OPTARG;;
    m) MIN_LENGTH=$OPTARG;;
    M) MAX_LENGTH=$OPTARG;;
    t) THREADS=$OPTARG;;
    l) LOG_DIR=$OPTARG;;
    :) printf "missing argument for -$OPTARG\n" >&2; exit 1;;
    \?) printf "invalid option for -$OPTARG\n" >&2; exit 1;;
  esac
done

# Check missing arguments
MISSING="is missing but required. Exiting."
if [ -z ${IN_DIR+x} ]; then echo "-d $MISSING"; echo "$USAGE"; exit 1; fi; 
if [ -z ${IN_REGEX+x} ]; then echo "-p $MISSING"; echo "$USAGE"; exit 1; fi; 
if [ -z ${OUT_DIR+x} ]; then echo "-o $MISSING"; echo "$USAGE"; exit 1; fi; 
if [ -z ${FW2+x} ]; then echo "-F $MISSING"; echo "$USAGE"; exit 1; fi; 
if [ -z ${RV2+x} ]; then echo "-R $MISSING"; echo "$USAGE"; exit 1; fi; 
if [ -z ${MIN_LENGTH+x} ]; then echo "-m $MISSING"; echo "$USAGE"; exit 1; fi; 
if [ -z ${MAX_LENGTH+x} ]; then echo "-M $MISSING"; echo "$USAGE"; exit 1; fi; 
if [ -z ${THREADS+x} ]; then echo "-t is missing. Defaulting to 1 thread."; THREADS=1; fi;
if [ -z ${LOG_DIR+x} ]; then echo "-l is missing. Defaulting to ./logs."; THREADS=1; fi;

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
mkdir -p $LOG_DIR
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

