#!/bin/bash

# DESCRIPTION
# Script for splitting reads based on custom sequences
#   
# IMPLEMENTATION
#    author   Søren Karst (sorenkarst@gmail.com)
#             Ryan Ziels (ziels@mail.ubc.ca)
#    license  GNU General Public License
#
### Description ----------------------------------------------------------------

USAGE="
-- longread_umi split_reads: Locates custom sequence(s) in fastq reads and splits
   reads at middle position
   
usage: $(basename "$0" .sh) [-h] [-e value -n value ] (-d value -o dir -t value)
(-a file -p value -f string -k)

where:
    -h  Show this help text.
    -d  Raw fastq reads.
    -o  Output directory
    -t  Number of threads to use.
    -a  Fasta file containing adapter sequences to split reads by. 
    -p  Shift split positions according to first base in adapter sequences.
        E.g. Use -12 to split 12 bp before adapter sequence. Use 10 
        to split in the middle of adapter that is 20 bp long.
    -k  Keep temporary files for debugging.
"

### Terminal Arguments ---------------------------------------------------------

# Import user arguments
while getopts ':hzd:o:t:a:p:f:k' OPTION; do
  case $OPTION in
    h) echo "$USAGE"; exit 1;;
    d) SEQ_IN=$OPTARG;;
    o) OUT_DIR=$OPTARG;;
    t) THREADS=$OPTARG;;
    a) ADP_FILE=$OPTARG;;
    p) POS_COR=$OPTARG;;
    f) SEQ_TYPE=$OPTARG;;
    k) KEEP_TMP="YES";;
    :) printf "missing argument for -$OPTARG\n" >&2; exit 1;;
    \?) printf "invalid option for -$OPTARG\n" >&2; exit 1;;
  esac
done

# Check missing arguments
MISSING="is missing but required. Exiting."
if [ -z ${SEQ_IN+x} ]; then echo "-d $MISSING"; echo "$USAGE"; exit 1; fi; 
if [ -z ${OUT_DIR+x} ]; then echo "-o $MISSING"; echo "$USAGE"; exit 1; fi; 
if [ -z ${THREADS+x} ]; then echo "-t is missing. Defaulting to 1 thread."; THREADS=1; fi;
if [ -z ${ADP_FILE+x} ]; then echo "-a $MISSING"; echo "$USAGE"; exit 1; fi; 
if [ -z ${POS_COR+x} ]; then echo "-p $MISSING"; echo "$USAGE"; exit 1; fi;
if [ -z ${SEQ_TYPE+x} ]; then echo "-f $MISSING"; echo "$USAGE"; exit 1; fi; 

### Source commands and subscripts -------------------------------------
. $LONGREAD_UMI_PATH/scripts/dependencies.sh # Path to dependencies script

# Create output dir
mkdir $OUT_DIR

# Logging
LOG_DIR=$OUT_DIR

LOG_NAME="$LOG_DIR/longread_umi_split_seq_log_$(date +"%Y-%m-%d-%T").txt"
echo "longread_umi split_seq log" >> $LOG_NAME
longread_umi_version_dump $LOG_NAME
exec &> >(tee -a "$LOG_NAME")
exec 2>&1

### Split reads

# Format input name
SEQ_NAME=${SEQ_IN##*/}
SEQ_NAME=${SEQ_NAME%.*}

# Define split function
seq_split(){
  #######################################
  # Format input
  local SEQ_IN=$(cat)
  local OUT_DIR=$2
  local ADP_FILE=$3
  local SEQ_TYPE=$4
  local POS_COR=$5
  
  # Create working dir
  mkdir $OUT_DIR
  echo \
    "$SEQ_IN" \
    > $OUT_DIR/input_sequence
  local SEQ_IN=$OUT_DIR/input_sequence

  # Map adapters to sequences
  $BWA index $SEQ_IN

  $BWA \
    aln \
    $SEQ_IN \
    $ADP_FILE \
    -i 5 \
    -d 16 \
    -n 4 \
    -t 1 \
    -N \
    > $OUT_DIR/split.sai


  $BWA \
    samse \
    -n 100000000 \
    $SEQ_IN \
    $OUT_DIR/split.sai \
    $ADP_FILE |\
    $SAMTOOLS view -F 4 - \
    > $OUT_DIR/split.sam
  
  # Determine split positions
  cat \
    $OUT_DIR/split.sam |\
  $GAWK \
    -v PC="$POS_COR" \
    '
    # Cigar to end pos function
    function maplen(CIGAR,  LEN, CIGV, cigv, CIGT, cigt, i) {
      CIGV=CIGAR
      gsub("[A-Z]", ";", CIGV)
      gsub("^;|;$", "", CIGV)
      split(CIGV, cigv, ";")
      CIGT=CIGAR
      gsub("[0-9]*", ";", CIGT)
      gsub("^;|;$", "", CIGT)
      split(CIGT, cigt, ";")
    
      for (i in cigt){
        if (cigt[i] != "S" && cigt[i] != "I"){
          LEN=LEN+cigv[i]
        }
      }
      return(LEN)
    }
  
    # Terminal message
    NR==1 {
      print "[" strftime("%T") "] ### Splitting reads by adapter matches  ###" > "/dev/stderr";
      print "[" strftime("%T") "] Reading adapter mapping file..." > "/dev/stderr";
    }

    # Adapter mapping file
    NR==FNR{

      # Extract data from optional fields
      for (i=12; i <= NF; i++){
        # Find NM field and remove prefix (primary hit err)
        if($i ~ /^NM:i:/){sub("NM:i:", "", $i); ERR = $i};
        # Find secondary hit field, remove prefix and split hits
        if($i ~ /^XA:Z:/){sub("XA:Z:", "", $i); split($i, shits, ";")};
      }

      # Add primary hit to err and strand arrays
      RNAME=$3
      SPOS=$4
      STRAND=$2
      CIG=$6
      # Check if match is empty
      if(RNAME != ""){
        # Create position or update position if match error is lower
        if (!(RNAME in err && SPOS in err[RNAME] && ERR > err[RNAME][SPOS])){
          if(STRAND == "0" ){
            err[RNAME][SPOS] = ERR
            str[RNAME][SPOS] = "+"
            cig[RNAME][SPOS] = CIG
          } else if (STRAND == "16"){
            err[RNAME][SPOS] = ERR
            str[RNAME][SPOS] = "-"
            cig[RNAME][SPOS] = CIG
          }
        }
      }

      # Cleanup vars
      RNAME=""
      SPOS=""
      STRAND=""
      CIG=""
      ERR=""

      # Add secondary hit to split matrix
      for (i in shits){
        # Split hit in subdata (read, pos, cigar, err) and format
        split(shits[i], tmp, ",")
        RNAME=tmp[1]
        POS_STRAND=tmp[2]""
        STRAND=match(POS_STRAND,"-")
        sub("+|-", "", POS_STRAND)
        SPOS=POS_STRAND+0
        CIG=tmp[3]
        ERR=tmp[4]+0

        # Check if match is empty
        if(RNAME != ""){
          # Create position or update position if match error is lower
          if (!(RNAME in err && SPOS in err[RNAME] && ERR > err[RNAME][SPOS])){
            if(STRAND == 0 ){
              err[RNAME][SPOS] = ERR
              str[RNAME][SPOS] = "+"
              cig[RNAME][SPOS] = CIG
            } else if (STRAND == 1){
              err[RNAME][SPOS] = ERR
              str[RNAME][SPOS] = "-"
              cig[RNAME][SPOS] = CIG
            }
          }
        }

        # Cleanup vars
        RNAME=""
        SPOS=""
        STRAND=""
        CIG=""
        ERR=""
      }
      next
    } END {
    
      # Remove redundant matches
      for (r in str){
        for (s in str[r]){
          # Define vars
          RNAME = r
          POS = s
          LEN = maplen(cig[r][s])
          ERR = err[r][s]
          STRAND = str[r][s]

          # Format position
          if (STRAND == "+"){
            POS = POS + PC
          } else if (STRAND == "-"){
            POS = POS + LEN - 1 - PC
          }

          # Search for overlapping mappings and choose best
          for (i=POS-LEN; i <= POS+LEN; i++){
            if (RNAME in err_filt && i in err_filt[RNAME] && ERR > err_filt[RNAME][i]){
              break
            } else {
              delete err_filt[RNAME][i]
              delete str_filt[RNAME][i]
              delete cig_filt[RNAME][i]
            }
            if (i == POS+LEN){
              err_filt[RNAME][POS] = err[r][s]
              str_filt[RNAME][POS] = str[r][s]
              cig_filt[RNAME][POS] = cig[r][s]
            }
          }
        }
      }

      # Print clean split list
      for (r in err_filt){
        for (s in err_filt[r]){
          print r, s , err_filt[r][s], str_filt[r][s], cig_filt[r][s]
        }
      }

    }' \
  > $OUT_DIR/split_pos.txt
  
  # Split reads
  #-- Assumes split_pos.txt is sorted!
  #-- After split there will be 1 bp overlap between subseqs
  
  $GAWK \
    -v SPLIT_POS="$OUT_DIR/split_pos.txt" \
    -v SEQ_IN="$SEQ_IN" \
    -v SEQ_TYPE="$SEQ_TYPE" \
    '
      FILENAME == SPLIT_POS{
        split_array[$1][$2]++
      }
      FILENAME == SEQ_IN && FNR%4==1 && SEQ_TYPE == "fastq" {
        SEQ_NAME=substr($1,2)
        if (SEQ_NAME in split_array){
          # Store sequence in memory
          getline
          SEQ = $0
          getline
          getline
          QUAL = $0
          # Define start pos
          START_POS=1
          for (i in split_array[SEQ_NAME]){
            # Name sub seq
            SUB_SEQ_NAME="@" SEQ_NAME "_" START_POS "_" i
            # Subset sequence and quality
            LENGTH=i - START_POS + 1
            SUB_SEQ=substr(SEQ, START_POS, LENGTH)
            SUB_QUAL=substr(QUAL, START_POS, LENGTH)
            # Print subs sequence
            print SUB_SEQ_NAME "\n" SUB_SEQ "\n" "+" "\n" SUB_QUAL
            START_POS=i
          }
        }
      }
      FILENAME == SEQ_IN && FNR%2==1 && SEQ_TYPE == "fasta" {
        SEQ_NAME=substr($1,2)
        if (SEQ_NAME in split_array){
          # Store sequence in memory
          getline
          SEQ = $0
          # Define start pos
          START_POS=1
          for (i in split_array[SEQ_NAME]){
            # Name sub seq
            SUB_SEQ_NAME=">" SEQ_NAME "_" START_POS "_" i
            # Subset sequence
            LENGTH=i - START_POS + 1
            SUB_SEQ=substr(SEQ, START_POS, LENGTH)
            # Print subs sequence
            print SUB_SEQ_NAME "\n" SUB_SEQ 
            START_POS=i
          }
        }
      }
    ' \
    $OUT_DIR/split_pos.txt \
    $SEQ_IN
} ## End split_seq

export -f seq_split

cat $SEQ_IN |\
$GNUPARALLEL \
  --env seq_split\
  -N 100 \
  -L 4 \
  -j $THREADS \
  --pipe \
  "
  seq_split \
    $SEQ_IN \
    $OUT_DIR/{#} \
    $ADP_FILE \
    $SEQ_TYPE \
    $POS_COR; \
  if [ -z ${KEEP_TMP+x} ]; then
    rm -rf $OUT_DIR/{#}
  fi
  " \
> $OUT_DIR/${SEQ_NAME}_split.${SEQ_TYPE}
