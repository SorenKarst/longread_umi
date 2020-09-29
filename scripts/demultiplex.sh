#!/bin/bash

# DESCRIPTION
#    Script for demultiplexing of UMI consensus sequences based on 
#    custom barcodes. The script demultiplexes raw read data
#    and creates an assignment based on a consensus of the raw read 
#    assignments. This approach is the most effective strategy for 
#    optimizing yield. The script expects dual barcodes in a barcode file.
#    If the same barcode is used in both ends simply repeat barcode.
#    
# IMPLEMENTATION
#    author	Søren Karst (sorenkarst@gmail.com)
#               Ryan Ziels (ziels@mail.ubc.ca)
#    license	GNU General Public License
#
# TO DO

### Source commands and subscripts -------------------------------------
. $LONGREAD_UMI_PATH/scripts/dependencies.sh # Path to dependencies script

### Description ----------------------------------------------------------------

USAGE="
-- longread_umi demultiplex: Dual barcode demultiplexing

   Script for demultiplexing UMI consensus sequences based on 
   custom barcodes. The script demultiplexes raw read data
   and assigns the consensus sequences to a sample by majority vote
   of the raw read assignments. Post processing demultiplxing optimizes 
   consensus yield. The script expects dual barcodes in a barcode file.
   If the same barcode is used in both ends simply repeat barcode.

usage: $(basename "$0" .sh) [-h] (-c file -r file -u file -o dir -b file)
(-p string -n range -t value) 

where:
    -h  Show this help text.
    -c  UMI consensus sequences that need demultiplexing.
    -r  Raw read sequences that were used to generate
        the consensus sequences.
    -u  List of raw read names and UMI bin assignments.
    -o  Output directory.
    -b  File containing barcodes. 
        [Default = "$BARCODES"].
    -p  Barcode name prefix [Default = 'barcode'].
    -n  Barcode numbers used. [Default  = '1-120'].
    -t  Number of threads used.

"

### Terminal Arguments ---------------------------------------------------------

# Import user arguments
while getopts ':hzd:c:r:u:o:b:p:n:t:' OPTION; do
  case $OPTION in
    h) echo "$USAGE"; exit 1;;
    c) CONSENSUS_IN=$OPTARG;;
	r) RAW_IN=$OPTARG;;
	u) UMI_BIN_RESULTS=$OPTARG;;
	o) OUT_DIR=$OPTARG;;
	b) BARCODE_FILE=$OPTARG;;
	p) BARCODE_PREFIX=$OPTARG;;
	n) BARCODE_RANGES=$OPTARG;;
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
if [ -z ${THREADS+x} ]; then echo "-t is missing. Defaulting to 1 thread."; THREADS=1; fi;

# Demultiplexing of UMI consensus sequences sequenced with Nanopore barcodes
mkdir $OUT_DIR

# Raw read demultiplexing

# Format barcodes
BC=$(
  $GAWK \
    -v BR="$BARCODE_RANGES" \
    -v BP="$BARCODE_PREFIX" \
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
            FW_SEQ=$3 $4
            RV_SEQ=revcomp($6 $7)
            FW_SEQ_RC=revcomp($3 $4)
            RV_SEQ_RC=$6 $7
            print ">" $1 "\n" FW_SEQ "..." RV_SEQ
            print ">" $1 "\n" RV_SEQ_RC "..." FW_SEQ_RC  
          }
        }
      }
    ' $BARCODE_FILE
)

echo "$BC" > $OUT_DIR/barcodes_used.fa
 
demultiplex_cutadapt(){
  # Input
  local JOBNR=$1
  local OUTDIR=$2
  local BC=$3

  # Demultiplex
  # Input from stdin
  $GAWK '
    NR%4==1{print}
    NR%4==2{print substr($0,1,150)"XXXXXXXXXXX" substr($0,length($0)-149,150)}
    NR%4==3{print}
    NR%4==0{print substr($0,1,150)"XXXXXXXXXXX" substr($0,length($0)-149,150)}
  ' |\
  $CUTADAPT \
    -e 0.15 \
    -O 17 \
    -g file:<(echo "$BC") \
    -o $OUTDIR/${JOBNR}tmp_{name}.fastq \
    -
}

export -f demultiplex_cutadapt

cat $RAW_IN |\
$GNUPARALLEL \
  --env demultiplex_cutadapt \
  --progress \
  -j $THREADS \
  -L 4 \
  -N 100 \
  --roundrobin \
  --pipe \
  "demultiplex_cutadapt \
     {#} \
     $OUT_DIR \
     \"$BC\"
  "

# Assign consensus to samples based on raw data demultiplexing

$GAWK '
  NR==FNR{
    UMI=$1
  	sub(";.*", "", UMI)
	  READ=$2 # Reads split by porechop are disregarded
    u2r[READ]=UMI
	  next
  }
  FNR%4==1{
    READ=substr($1,2)
	  SAMPLE=FILENAME
	  sub(".*tmp_", "", SAMPLE)
	  sub(".fastq", "", SAMPLE)
    u2s[u2r[READ]][SAMPLE]++
	  utr[u2r[READ]]++
  }
  END{
    for (UMI in u2s){
	  if (UMI != ""){
	    UMI_MAX = 0
	    UMI_MAX_SAMPLE = ""
		UMI_PASS_FILTER="NO"
	    for (SAMPLE in u2s[UMI]){
	      if (u2s[UMI][SAMPLE]+0 >= UMI_MAX){
		    UMI_MAX = u2s[UMI][SAMPLE]+0
		    UMI_MAX_SAMPLE = SAMPLE
			UMI_RATIO=u2s[UMI][SAMPLE]/utr[UMI]
			if(UMI_RATIO > 0.65){UMI_PASS_FILTER="YES"}
		  }
		}
	    print UMI, UMI_MAX_SAMPLE, UMI_MAX, UMI_RATIO, UMI_PASS_FILTER
	  }
	}
  }
' $UMI_BIN_RESULTS $OUT_DIR/*${BARCODE_PREFIX}*fastq > $OUT_DIR/demultiplex.txt
rm $OUT_DIR/*tmp_*fastq

# Demultiplexing consensus
$GAWK -v OUT_DIR="$OUT_DIR" '
  NR==FNR{
    if($5 == "YES"){
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
