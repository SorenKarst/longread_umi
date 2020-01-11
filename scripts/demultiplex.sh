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
#               Ryans Ziels (ziels@mail.ubc.ca)
#    license	GNU General Public License
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
    -o $OUTDIR/${JOBNR}_{name}.fastq \
    -
}

export -f demultiplex_cutadapt

cat $RAW_IN |\
$GNUPARALLEL \
  --progress \
  -j $THREADS \
  -L4 \
  -N100 \
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
	  sub(".*_", "", SAMPLE)
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
rm ./$OUT_DIR/*.fastq

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
