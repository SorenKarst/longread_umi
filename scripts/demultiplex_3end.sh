#!/bin/bash

# DESCRIPTION
#    Script for demultiplexing of UMI consensus sequences based on 
#    custom barcodes. The script demultiplexes raw read data
#    and creates an assignment based on a consensus of the raw read 
#    assignments. This approach is the most effective strategy for 
#    optimizing yield. The script expects dual barcodes in a barcode file.
#    If the same barcode is used in both ends simply repeat barcode. This
#    version of the script only looks for barcodes in the 3' end. This is
#    if the raw data is truncated in the 5' end. Dual barcoding is still used.
#    
# IMPLEMENTATION
#    author  Søren Karst (sorenkarst@gmail.com)
#            Ryan Ziels (ziels@mail.ubc.ca)
#    license  GNU General Public License
#
# TO DO

### Source commands and subscripts -------------------------------------
. $LONGREAD_UMI_PATH/scripts/dependencies.sh # Path to dependencies script

### Description ----------------------------------------------------------------

USAGE="
-- longread_umi demultiplex_3end: 3'-end dual barcode demultiplexing

    Script for demultiplexing UMI consensus sequences based on 
    custom barcodes. The script demultiplexes raw read data
    and assigns the consensus sequences to a sample by majority vote
    of the raw read assignments. Post processing demultiplxing optimizes 
    consensus yield. The script expects dual barcodes in a barcode file.
    If the same barcode is used in both ends simply repeat barcode. This
    version of the script only looks for barcodes in the 3' end. This demultiplexing
    is for data types which are truncated in the 5' end. Dual barcoding is
    still used.
	
usage: $(basename "$0" .sh) [-h] (-c file -r file -u file -o dir -b file -p string)
(-n range -m value -t value) 

where:
    -h  Show this help text.
    -c  UMI consensus sequences that need demultiplexing.
    -r  Raw read sequences that were used to generate
        the consensus sequences.
    -u  List of raw read names and UMI bin assignments.
    -o  Output directory.
    -b  File containing barcodes.
        Default is "$BARCODES"
    -p  Barcode name prefix. [Default = 'barcode'].
    -n  Barcode numbers used. [Default = '1-120'].
    -m  Minimum number of barcodes found to demultiplex
        sequences. Default 2.
    -t  Number of threads used.

"

### Terminal Arguments ---------------------------------------------------------

# Import user arguments
while getopts ':hzd:c:r:u:o:b:p:n:m:t:' OPTION; do
  case $OPTION in
    h) echo "$USAGE"; exit 1;;
    c) CONSENSUS_IN=$OPTARG;;
	r) RAW_IN=$OPTARG;;
	u) UMI_BIN_RESULTS=$OPTARG;;
	o) OUT_DIR=$OPTARG;;
	b) BARCODE_FILE=$OPTARG;;
	p) BARCODE_PREFIX=$OPTARG;;
	n) BARCODE_RANGES=$OPTARG;;
	m) MIN_BARCODE_COUNT=$OPTARG;;
    t) THREADS=$OPTARG;;
    :) printf "missing argument for -$OPTARG\n" >&2; exit 1;;
    \?) printf "invalid option for -$OPTARG\n" >&2; exit 1;;
  esac
done

# Check missing arguments
MISSING="is missing but required. Exiting."
if [ -z ${CONSENSUS_IN+x} ]; then echo "-c $MISSING"; echo "$USAGE"; exit 1; fi; 
if [ -z ${RAW_IN+x} ]; then echo "-r $MISSING"; echo "$USAGE"; exit 1; fi; 
if [ -z ${UMI_BIN_RESULTS+x} ]; then echo "-u $MISSING"; echo "$USAGE"; exit 1; fi; 
if [ -z ${OUT_DIR+x} ]; then echo "-o $MISSING"; echo "$USAGE"; exit 1; fi; 
if [ -z ${BARCODE_FILE+x} ]; then
  echo "-b missing. Defaulting to $BARCODES."
  BARCODE_FILE="$BARCODES"
fi;
if [ -z ${BARCODE_PREFIX+x} ]; then
  echo "-p missing. Defaulting to 'barcode'."
  BARCODE_PREFIX='barcode'
fi;
if [ -z ${BARCODE_RANGES+x} ]; then
  echo "-n missing. Defaulting to barcode numbers '1-120'."
  BARCODE_RANGES='1-120'
fi
if [ -z ${MIN_BARCODE_COUNT+x} ]; then
  echo "-m missing. Defaulting to minimum barcode count of 2."
  MIN_BARCODE_COUNT=2
fi
if [ -z ${THREADS+x} ]; then echo "-t is missing. Defaulting to 1 thread."; THREADS=1; fi;

# Demultiplexing of UMI consensus sequences sequenced with Nanopore barcodes
mkdir $OUT_DIR

# Raw read demultiplexing

# Format barcodes
$GAWK \
  -v BR="$BARCODE_RANGES" \
  -v BP="$BARCODE_PREFIX" \
  -v OUT="$OUT_DIR" \
  ' 
    # Define function for reverse complement
    function revcomp(s,  i, o) {
      o = ""
      for(i = length; i > 0; i--)
           o = o c[substr(s, i, 1)]
      return(o)
    }
    BEGIN{
      # Define revcomp mapping vector
      c["A"] = "T"
      c["C"] = "G"
      c["G"] = "C"
      c["T"] = "A"
      
      # Format barcode ranges for testing
      split(BR, RANGES, ";")
      for (R in RANGES){
        FW = RANGES[R]
        sub("-.*", "", FW)
        RV = RANGES[R]
        sub(".*-", "", RV)
        FW_L[++j]=FW+0
        RV_L[j]=RV+0
      }
    }
    # Store barcodes matching ranges
    NR>=2{
      BNR=$1
      sub(BP, "", BNR)
      for (R in FW_L){
        if ( BNR+0 >= FW_L[R] && BNR+0 <= RV_L[R]){
          # Store barcode sequences
          BL[$2]=$3 $4
          BL[$5]=$6 $7
      # Store barcode combinations
      DBL[$1]=$2 "\t" $5
        }
      }
    }
    END{
      for (B in BL){
        print ">" B "\n" BL[B] > OUT "/barcode_sequences.fa"
      }
    for (DB in DBL){
      print DB "\t" DBL[DB] > OUT "/barcode_combinations.fa"
    }
    }
  ' $BARCODE_FILE

demultiplex_cutadapt_5(){
  # Input
  local JOBNR=$1
  local OUTDIR=$2
  local BC=$3

  # Demultiplex
  # Input from stdin
  $GAWK '
    NR%4==1{print}
    NR%4==2{print substr($0,1,150)}
    NR%4==3{print}
    NR%4==0{print substr($0,1,150)}
  ' |\
  $CUTADAPT \
    -e 0.15 \
    -O 17 \
    -g file:"$BC" \
    -o $OUTDIR/${JOBNR}_{name}.5.fastq \
    - 
}

demultiplex_cutadapt_3(){
  # Input
  local JOBNR=$1
  local OUTDIR=$2
  local BC=$3

  # Demultiplex
  # Input from stdin
  $GAWK '
    NR%4==1{print}
    NR%4==2{print substr($0,length($0)-149,150)}
    NR%4==3{print}
    NR%4==0{print substr($0,length($0)-149,150)}
  ' |\
  $CUTADAPT \
    -e 0.15 \
    -O 17 \
    -a file:<(cat "$BC" | $SEQTK seq -r -) \
    -o $OUTDIR/${JOBNR}_{name}.3.fastq \
    - 
}    

export -f demultiplex_cutadapt_5
export -f demultiplex_cutadapt_3

# Find 5' barcodes
cat $RAW_IN |\
$GNUPARALLEL \
  --env demultiplex_cutadapt_5 \
  --progress \
  -j $THREADS \
  -L4 \
  -N100 \
  --roundrobin \
  --pipe \
  "demultiplex_cutadapt_5 \
     {#} \
     $OUT_DIR \
     $OUT_DIR/barcode_sequences.fa
  "
  
# Find 3' barcodes
cat $RAW_IN |\
$GNUPARALLEL \
  --env demultiplex_cutadapt_3 \
  --progress \
  -j $THREADS \
  -L4 \
  -N100 \
  --roundrobin \
  --pipe \
  "demultiplex_cutadapt_3 \
     {#} \
     $OUT_DIR \
     $OUT_DIR/barcode_sequences.fa
  "


# Assign consensus to samples based on raw data demultiplexing

$GAWK \
  -v UMI_MAP_FILE="$UMI_BIN_RESULTS" \
  -v BARCODE_COMBINATION_FILE="$OUT_DIR/barcode_combinations.fa" \
'
  # Import UMI -> read connection
  FILENAME==UMI_MAP_FILE{
    UMI=$1
    sub(";.*", "", UMI)
    READ=$2 # Reads split by porechop are disregarded
    u2r[READ]=UMI
    next
  }
  # Import Barcode combinations
  FILENAME==BARCODE_COMBINATION_FILE{
    BARCODE_LIST[$2"|"$3]=$1
    next
  }
  # Import read -> barcode connection
  FNR%4==1 && FILENAME !~ /unknown/{
  # Format readname
    READ=substr($1,2)
  # Format barcode name
  BARCODE=FILENAME
  sub("^.*/[0-9]*_","", BARCODE)
  sub("\\..*fastq", "", BARCODE)
    u2s[u2r[READ]][BARCODE]++
  utr[u2r[READ]]++
  }
  END{
    # Loop over UMIs
    for (UMI in u2s){
    if (UMI != ""){
      UMI_MAX1 = 0
      UMI_MAX2 = 0
      UMI_MAX1_BARCODE = ""
      UMI_MAX2_BARCODE = ""
      # Loop over barcodes assigned to UMI
      for (BARCODE in u2s[UMI]){
        # Check if barcode has most assignments
        if (u2s[UMI][BARCODE]+0 >= UMI_MAX1){
          # Move current highest to second highest
          UMI_MAX2 = UMI_MAX1
          UMI_MAX2_BARCODE = UMI_MAX1_BARCODE
          # Assign new highest
          UMI_MAX1 = u2s[UMI][BARCODE]+0
          UMI_MAX1_BARCODE = BARCODE
        } else if (u2s[UMI][BARCODE]+0 >= UMI_MAX2){
          #Check if barcode has second most assignments
          UMI_MAX2 = u2s[UMI][BARCODE]+0
          UMI_MAX2_BARCODE = BARCODE
        }
      }
      
      # Assign to barcode
      BARCODE_COMBO1=UMI_MAX1_BARCODE"|"UMI_MAX2_BARCODE
      BARCODE_COMBO2=UMI_MAX2_BARCODE"|"UMI_MAX1_BARCODE
      if (BARCODE_COMBO1 in BARCODE_LIST){
        BARCODE_COMBO=BARCODE_LIST[BARCODE_COMBO1]
      } else if(BARCODE_COMBO2 in BARCODE_LIST){
        BARCODE_COMBO=BARCODE_LIST[BARCODE_COMBO2]
      } else{BARCODE_COMBO="?"}
      
      # Calculate assignment ratio and determine if passing
      UMI_RATIO_1=UMI_MAX1/utr[UMI]
      UMI_RATIO_2=UMI_MAX2/utr[UMI]
      if(UMI_RATIO_1 > 0.15 && UMI_RATIO_2 > 0.15 && UMI_RATIO_1 + UMI_RATIO_2 > 0.90 && utr[UMI] > 5 && BARCODE_COMBO != "?"){
        UMI_PASS_FILTER="YES"
      } else {UMI_PASS_FILTER="NO"}
      
      # Print results for UMI
      print UMI, BARCODE_COMBO, UMI_MAX1_BARCODE, UMI_MAX2_BARCODE, UMI_MAX1, UMI_MAX2, utr[UMI], UMI_PASS_FILTER
      }
    }
  }
' \
$UMI_BIN_RESULTS \
$OUT_DIR/barcode_combinations.fa \
$OUT_DIR/*fastq \
> $OUT_DIR/demultiplex.txt

rm $OUT_DIR/*.fastq

# Demultiplexing consensus
$GAWK -v OUT_DIR="$OUT_DIR" '
  NR==FNR{
    if($8 == "YES"){
    umi[$1]=$2
  }
  next
  }
  FNR%2==1 {
    UMI=substr($0,2);
  sub("bins.*", "", UMI)
    if (UMI in umi){
    print $0 > OUT_DIR "/" umi[UMI] ".fa"
      getline; print > OUT_DIR "/" umi[UMI] ".fa"
    } else {
    print $0 > OUT_DIR "/undetermined.fa"
    getline; print > OUT_DIR "/undetermined.fa"
  }
  }
' $OUT_DIR/demultiplex.txt $CONSENSUS_IN
