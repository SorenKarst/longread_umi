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

### Terminal input ------------------------------------------------------------
READ_IN=${1:-reads.fq}
OUT_DIR=${2:-primer_position}
THREADS=${3:-60}
FW2=${4:-AGRGTTYGATYMTGGCTCAG} #RC: CTGAGCCAKRATCRAACYCT
RV2=${5:-CGACATCGAGGTGCCAAAC} #RC: GTTTGGCACCTCGATGTCG
TERMINAL_LENGTH=${6:-500}
SUBSET=${7:-10000000}
FW1=${8:-CAAGCAGAAGACGGCATACGAGAT}
RV1=${9:-AATGATACGGCGACCACCGAGATC}

### Source commands and subscripts -------------------------------------
. $LONGREAD_UMI_PATH/scripts/dependencies.sh # Path to dependencies script

### Determine primer mapping positions

# Create output dir
mkdir $OUT_DIR

# Extract adaptor region
head -n $SUBSET $READ_IN |\
$GAWK -v BD="$OUT_DIR" -v TL="$TERMINAL_LENGTH" '
    NR%4==1{
       print ">" substr($1,2) > BD"/reads_t1.fa";
       print ">" substr($1,2) > BD"/reads_t2.fa"; 
     }
     NR%4==2{
       print substr($0, 1, TL) > BD"/reads_t1.fa";
       print substr($0, length($0) - TL + 1, TL)  > BD"/reads_t2.fa"; 
     }
'
# Determine primer start position from either end
primer_start_pos(){
  $CUTADAPT \
    -j $THREADS \
    -O 10 \
	-e 0.2\
    -a $FW2 \
    -a $RV2 \
    - \
    --discard-untrimmed \
	2>> $OUT_DIR/umi_endpos_log.txt |\
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
  > $OUT_DIR/reads_t1_umi_endpos.txt

$SEQTK seq -r $OUT_DIR/reads_t2.fa |\
  primer_start_pos \
  > $OUT_DIR/reads_t2_umi_endpos.txt
  
# Determine fraction of reads with terminal adaptors
adp_presence(){
  $CUTADAPT \
    -j 1 \
    -O 10 \
	-e 0.2\
    -g $FW1 \
    -g $RV1 \
    -G $FW1 \
    -G $RV1 \
    $1 \
	$2 \
    --discard-untrimmed \
	-p /dev/null \
	-o /dev/null \
    > $OUT_DIR/adp_presence_log.txt
}

adp_presence \
  $OUT_DIR/reads_t1.fa \
  <($SEQTK seq -r $OUT_DIR/reads_t2.fa)

#adp_presence2(){
#  $CUTADAPT \
#    -j $THREADS \
#    -O 10 \
#	-e 0.2\
#    -a $FW1 \
#    -a $RV1 \
#    - \
#    --discard-untrimmed \
#	1> /dev/null \
#	2>> $OUT_DIR/adp_presence_log.txt
#}

#cat $OUT_DIR/reads_t1.fa |\
#  adp_presence2 
#
#$SEQTK seq -r $OUT_DIR/reads_t2.fa |\
#  adp_presence2 
