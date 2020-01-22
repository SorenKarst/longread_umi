#!/bin/bash

# DESCRIPTION
#    Script for demultiplexing of UMI consensus sequences based on 
#    custom barcodes. The script demultiplexes raw read data
#    and creates an assignment based on a consensus of the raw read 
#    assignments. This approach is the most effective strategy for 
#    optimizing yield. The script expects dual barcodes in a barcode file.
#    If the same barcode is used in both ends simply repeat barcode. This
#    version of the script only looks for barcodes in the 3' end. This is
#    if the raw data is truncated in the 5' end. Dualo barcoding is still used.
#    
# IMPLEMENTATION
#    author  Søren Karst (sorenkarst@gmail.com)
#            Ryan Ziels (ziels@mail.ubc.ca)
#    license  GNU General Public License
#
# TO DO

### Source commands and subscripts -------------------------------------
. $LONGREAD_UMI_PATH/scripts/dependencies.sh # Path to dependencies script

### Terminal input ------------------------------------------------------------
CONSENSUS_IN=$1
THREADS=$2
RAW_IN=${3:-reads.fq}
UMI_BIN_RESULTS=${4:-umi_binning/read_binning/umi_bin_map.txt}
OUT_DIR=${5:-demultiplexing}
BARCODE_FILE=${6:-$BARCODES}
BARCODE_PREFIX=${7:-"barcode"}
BARCODE_RANGES=${8:-"1-120"}
MIN_BARCODE_COUNT=${9:-2}

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
      print DB, DBL[DB] > OUT "/barcode_combinations.fa"
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

rm ./$OUT_DIR/*.fastq

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
