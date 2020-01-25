#!/bin/bash

# DESCRIPTION
#    Script for checking position of gene specific primers
#    in read terminals when used with custom library preparation.
#    Custom library preparation includes change or modification to
#    any adaptors and/or use of Nanopore barcoding kits. 
#   
# IMPLEMENTATION
#    author   Søren Karst (sorenkarst@gmail.com)
#             Ryan Ziels (ziels@mail.ubc.ca)
#    license  GNU General Public License
#
### Description ----------------------------------------------------------------

USAGE="
-- longread_umi primer_position: Locate adapter and primer positions in read data
    Script for checking position of adapters and gene specific primers flanking
    UMI sequences in read terminals. Relevant if using custom UMI adapters/primers,
    sample barcoding or if basecalling/processing truncates reads.
   
usage: $(basename "$0" .sh) [-h] [-e value -n value ] (-d value -o dir -t value)
(-f string -F string -r string -R string ) 

where:
    -h  Show this help text.
    -d  Raw fastq reads.
    -o  Output directory
    -t  Number of threads to use.
    -f  Forward adaptor sequence. 
    -F  Forward primer sequence.
    -r  Reverse adaptor sequence.
    -R  Reverse primer sequence.
    -e  Length of terminal end to search for primers. [Default = 500]
    -n  Subset reads before search. [Default = 100000]
"

### Terminal Arguments ---------------------------------------------------------

# Import user arguments
while getopts ':hzd:o:t:f:F:r:R:e:n:' OPTION; do
  case $OPTION in
    h) echo "$USAGE"; exit 1;;
    d) READ_IN=$OPTARG;;
    o) OUT_DIR=$OPTARG;;
    t) THREADS=$OPTARG;;
    f) FW1=$OPTARG;;
    F) FW2=$OPTARG;;
    r) RV1=$OPTARG;;
    R) RV2=$OPTARG;;
    e) TERMINAL_LENGTH=$OPTARG;;
    n) SUBSET=$OPTARG;;
    :) printf "missing argument for -$OPTARG\n" >&2; exit 1;;
    \?) printf "invalid option for -$OPTARG\n" >&2; exit 1;;
  esac
done

# Check missing arguments
MISSING="is missing but required. Exiting."
if [ -z ${READ_IN+x} ]; then echo "-d $MISSING"; echo "$USAGE"; exit 1; fi; 
if [ -z ${OUT_DIR+x} ]; then echo "-o $MISSING"; echo "$USAGE"; exit 1; fi; 
if [ -z ${THREADS+x} ]; then echo "-t is missing. Defaulting to 1 thread."; THREADS=1; fi;
if [ -z ${FW1+x} ]; then echo "-f $MISSING"; echo "$USAGE"; exit 1; fi; 
if [ -z ${FW2+x} ]; then echo "-F $MISSING"; echo "$USAGE"; exit 1; fi; 
if [ -z ${RV1+x} ]; then echo "-r $MISSING"; echo "$USAGE"; exit 1; fi; 
if [ -z ${RV2+x} ]; then echo "-R $MISSING"; echo "$USAGE"; exit 1; fi; 
if [ -z ${TERMINAL_LENGTH+x} ]; then echo "-e missing. Defaulting to 500 bp."; TERMINAL_LENGTH=500; fi;
if [ -z ${SUBSET+x} ]; then echo "-n missing. Defaulting to 100,000 sequences"; SUBSET=100000; fi;

### Source commands and subscripts -------------------------------------
. $LONGREAD_UMI_PATH/scripts/dependencies.sh # Path to dependencies script

### Determine primer mapping positions

# Create output dir
mkdir $OUT_DIR

# Extract adaptor region
head -n $(( SUBSET*4 )) $READ_IN |\
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
$SEQTK seq -r $OUT_DIR/reads_t2.fa > $OUT_DIR/reads_t2_rc.fa

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

cat $OUT_DIR/reads_t2_rc.fa |\
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
  $OUT_DIR/reads_t2_rc.fa

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
