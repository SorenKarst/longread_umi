#!/bin/bash

# DESCRIPTION
#    Script for demultiplexing of UMI consensus sequences based on 
#    Nanopore native barcodes. The script demultiplexes raw read data
#    and creates an assignment based on a consensus of the raw read 
#    assignments. This approach is the most effective strategy for 
#    optimizing yield.
#    
# IMPLEMENTATION
#    author	Søren Karst (sorenkarst@gmail.com)
#               Ryans Ziels (ziels@mail.ubc.ca)
#    license	GNU General Public License
#
# TO DO


### Terminal input ------------------------------------------------------------
CONSENSUS_IN=$1
THREADS=$2
RAW_IN=${3:-reads.fq}
UMI_BIN_RESULTS=${4:-umi_binning/read_binning/umi_bin_map.txt}
OUT_DIR=${5:-demultiplexing}
BARCODE_FILE=${6:-$NANOPORE_BARCODES}

# Demultiplexing of UMI consensus sequences sequenced with Nanopore barcodes
mkdir $OUT_DIR

# Raw read demultiplexing
BC=$($GAWK '
      NR>=2{print ">" $1 "\nAATGTACTTCGTTCAGTTACGTATTGCTAAGGTTAA" $2 "CAGCACCT"}
     ' $BARCODE_FILE)
BCRV=$(echo "$BC" |\
  $SEQTK seq -r -)
  

demultiplex_cutadapt(){
  # Input
  local JOBNR=$1
  local OUTDIR=$2
  local BC=$3
  local BCRV=$4

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
    -a file:<(echo "$BCRV") \
    -o $OUTDIR/${JOBNR}_{name}.fastq \
    -
}

export -f demultiplex_cutadapt

cat $RAW_IN |\
$GNUPARALLEL \
  --progress \
  -j $THREADS \
  -L4 \
  -N1 \
  --roundrobin \
  --pipe \
  "demultiplex_cutadapt \
     {#} \
     $OUT_DIR \
     \"$BC\" \
	 \"$BCRV\"
  "

# Assign consensuus to samples based on raw data demultiplexing

$GAWK '
  NR==FNR{
    UMI=$1
	sub(";.*", "", UMI)
	READ=$2
	sub("_.*", "", READ)
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
' $UMI_BIN_RESULTS $OUT_DIR/*BC*fastq > $OUT_DIR/sample_assignment.txt
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
' $OUT_DIR/sample_assignment.txt $CONSENSUS_IN
