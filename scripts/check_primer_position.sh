#!/bin/bash

# DESCRIPTION
#    Script for checking position of gene specific primers
#    in read terminals when used with custom library preparation.
#    Custom library preparation includes change or modification to
#    any adaptors and/or use of Nanopore barcoding kits. 
#   
# IMPLEMENTATION
#    author	Søren Karst (sorenkarst@gmail.com)
#               Ryans Ziels (ziels@mail.ubc.ca)
#    license	GNU General Public License
#

### Source commands and subscripts -------------------------------------
export PIPELINE_PATH="$(dirname "$(readlink -f "$0")")"
. $PIPELINE_PATH/dependencies.sh # Path to dependencies script

### Terminal input ------------------------------------------------------------
READ_IN=${1:-reads.fq}
OUT_DIR=${2:-primer_position}
THREADS=${3:-60}
FW2=${4:-AGRGTTYGATYMTGGCTCAG} #RC: CTGAGCCAKRATCRAACYCT
RV2=${5:-CGACATCGAGGTGCCAAAC} #RC: GTTTGGCACCTCGATGTCG
TERMINAL_LENGTH=${6:-500}

### Determine primer mapping positions

# Create output dir
mkdir $OUT_DIR

# Extract adaptor region
$GAWK -v BD="$OUT_DIR" -v TL="$TERMINAL_LENGTH" '
    NR%4==1{
       print ">" substr($1,2) > BD"/reads_t1.fa";
       print ">" substr($1,2) > BD"/reads_t2.fa"; 
     }
     NR%4==2{
       print substr($0, 1, TL) > BD"/reads_t1.fa";
       print substr($0, length($0) - TL + 1, TL)  > BD"/reads_t2.fa"; 
     }
' $READ_IN

# Determine primer start position from either end
primer_start_pos(){
  $CUTADAPT \
    -j $THREADS \
    -O 10 \
    -a $FW2 \
    -a $RV2 \
    - \
    --discard-untrimmed \
    --quiet |\
  awk '
    NR%2==0{
      # Position
      p=length($0)+0
      c[p]++
      if(pm+0 < p){pm = p}
      # SD and MEAN
      sum+=p
      sums+=p*p 
      n++
    } END {
      for(j=0; j <= pm; j++){
        print j, c[j]+0 
      }
      print "\nMean\tsd"
      print sum/n"\t"sqrt(sums/n -(sum/n)^2)   
    }
  '
}

cat $OUT_DIR/reads_t1.fa |\
  primer_start_pos \
  > $OUT_DIR/reads_t1_pos.txt

$SEQTK seq -r $OUT_DIR/reads_t2.fa |\
  primer_start_pos \
  > $OUT_DIR/reads_t2_pos.txt

