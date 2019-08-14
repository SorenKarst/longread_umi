#!/bin/bash
# DESCRIPTION
#    longread_umi nanopore_pipeline script. 
#    
# IMPLEMENTATION
#    author	Søren Karst (sorenkarst@gmail.com)
#               Ryans Ziels (ziels@mail.ubc.ca)
#    license	GNU General Public License
#
# To-do:
# - Fix logging

### Description ----------------------------------------------------------------

USAGE="$(basename "$0" .sh) [-h] [-d file -n value -c value -o dir -s value -e value 
-m value -M value -f string -F string -r string -R string -t value -T ] 
-- longread_umi nanopore_pipeline: Generates UMI consensus sequences from
   raw Nanopore fastq reads with UMIs in both terminals.

where:
    -h  Show this help text.
    -d  Single file containing raw Nanopore data in fastq format.
    -c  Minimum read coverage for using UMI consensus sequences for 
        variant calling.
    -o  Output directory.
    -s  Check start of read up to s bp for UMIs.
    -e  Check end of read up to f bp for UMIs.
    -m  Minimum read length.
    -M  Maximum read length.
    -f  Forward adaptor sequence. 
    -F  Forward primer sequence.
    -r  Reverse adaptor sequence.
    -R  Reverse primer sequence.
    -n  Process n number of bins. If not defined all bins are processed.
        Pratical for testing large datasets.
    -w  Use predefined workflow with settings for s, e, m, M, f, F, r, R.
        rrna_operon [70, 80, 3500, 6000, CAAGCAGAAGACGGCATACGAGAT,
        AGRGTTYGATYMTGGCTCAG, AATGATACGGCGACCACCGAGATC, CGACATCGAGGTGCCAAAC]
    -t  Number of threads to use.
    -T  Use all available cores for Medaka consensus jobs. One job is started
        pr. threads (-t). With -T all cores are available to all jobs.
        Can drasticly speed up consensus calling but use with care.
"

### Terminal Arguments ---------------------------------------------------------

# Import user arguments
while getopts ':hzd:c:o:s:e:m:M:f:F:r:R:n:w:t:T' OPTION; do
  case $OPTION in
    h) echo "$USAGE"; exit 1;;
    d) INPUT_READS=$OPTARG;;
    c) UMI_COVERAGE_MIN=$OPTARG;;
    o) OUT_DIR=$OPTARG;;
    s) START_READ_CHECK=$OPTARG;;
    e) END_READ_CHECK=$OPTARG;;
    m) MIN_LENGTH=$OPTARG;;
    M) MAX_LENGTH=$OPTARG;;
    f) FW1=$OPTARG;;
    F) FW2=$OPTARG;;
    r) RV1=$OPTARG;;
    R) RV2=$OPTARG;;  
    n) UMI_SUBSET_N=$OPTARG;;
    w) WORKFLOW=$OPTARG;;
    t) THREADS=$OPTARG;;
    T) TURBO="YES";;
    :) printf "missing argument for -$OPTARG\n" >&2; exit 1;;
    \?) printf "invalid option for -$OPTARG\n" >&2; exit 1;;
  esac
done

# Check missing arguments
MISSING="is missing but required. Exiting."
if [ "$WORKFLOW" == rrna_operon ]; then
  START_READ_CHECK=70
  END_READ_CHECK=80
  MIN_LENGTH=3500
  MAX_LENGTH=6000
  FW1=CAAGCAGAAGACGGCATACGAGAT
  FW2=AGRGTTYGATYMTGGCTCAG
  RV1=AATGATACGGCGACCACCGAGATC
  RV2=CGACATCGAGGTGCCAAAC
fi
if [ -z ${INPUT_READS+x} ]; then echo "-d $MISSING"; echo "$USAGE"; exit 1; fi; 
if [ -z ${UMI_COVERAGE_MIN+x} ]; then echo "-c $MISSING"; echo "$USAGE"; exit 1; fi;
if [ -z ${OUT_DIR+x} ]; then echo "-o $MISSING"; echo "$USAGE"; exit 1; fi;
if [ -z ${START_READ_CHECK+x} ]; then echo "-s $MISSING"; echo "$USAGE"; exit 1; fi;
if [ -z ${END_READ_CHECK+x} ]; then echo "-e $MISSING"; echo "$USAGE"; exit 1; fi;
if [ -z ${MIN_LENGTH+x} ]; then echo "-m $MISSING"; echo "$USAGE"; exit 1; fi;
if [ -z ${MAX_LENGTH+x} ]; then echo "-M $MISSING"; echo "$USAGE"; exit 1; fi;
if [ -z ${FW1+x} ]; then echo "-f $MISSING"; echo "$USAGE"; exit 1; fi;
if [ -z ${FW2+x} ]; then echo "-F $MISSING"; echo "$USAGE"; exit 1; fi;
if [ -z ${RV1+x} ]; then echo "-r $MISSING"; echo "$USAGE"; exit 1; fi;
if [ -z ${RV2+x} ]; then echo "-R $MISSING"; echo "$USAGE"; exit 1; fi;
if [ -z ${THREADS+x} ]; then echo "-t is missing. Defaulting to 1 thread."; THREADS=1; fi;
if [ -z ${TURBO+x} ]; then echo "-T is missing. Turbo is disabled."; TURBO=NO; fi;

### Source commands and subscripts -------------------------------------
. $LONGREAD_UMI_PATH/scripts/dependencies.sh # Path to dependencies script

### Pipeline -----------------------------------------------------------
# Logging
LOG_DIR=logs
mkdir $LOG_DIR

LOG_NAME="$LOG_DIR/longread_umi_nanopore_pipeline_log_$(date +"%Y-%m-%d-%T").txt"
echo "longread_umi nanopore_pipeline log" >> $LOG_NAME
longread_umi_version_dump $LOG_NAME
exec &> >(tee -a "$LOG_NAME")
exec 2>&1
echo ""
echo "### Settings:"
echo "Input reads: $INPUT_READS"
echo "Bin size cutoff: $UMI_COVERAGE_MIN"
echo "Output directory: $OUT_DIR"
echo "Check start of read: $START_READ_CHECK"
echo "Check end of read: $END_READ_CHECK"
echo "Minimum read length: $MIN_LENGTH"
echo "Maximum read length: $MAX_LENGTH"
echo "Forward adaptor sequence: $FW1"
echo "Forward primer sequence: $FW2"
echo "Reverse adaptor sequence: $RV1"
echo "Reverse adaptor primer: $RV2" 
echo "UMI subsampling: $UMI_SUBSET_N"
echo "Preset workflow: $WORKFLOW"
echo "Threads: $THREADS"
echo "Turbo: $TURBO"
echo ""

# Read filtering and UMI binning
UMI_DIR=$OUT_DIR/umi_binning

longread_umi umi_binning  \
  -d $INPUT_READS      `# Raw nanopore data in fastq format`\
  -o umi_binning       `# Output folder`\
  -m $MIN_LENGTH       `# Min read length`\
  -M $MAX_LENGTH       `# Max read length` \
  -s $START_READ_CHECK `# Start of read to check` \
  -e $END_READ_CHECK   `# End of read to check` \
  -f $FW1              `# Forward adaptor sequence` \
  -F $FW2              `# Forward primer sequence` \
  -r $RV1              `# Reverse adaptor sequence` \
  -R $RV2              `# Reverse primer sequence` \
  -t $THREADS          `# Number of threads`

# Sample UMI bins for testing
if [ ! -z ${UMI_SUBSET_N+x} ]; then
  find umi_binning/read_binning/bins \
    -name 'umi*bins.fastq' | sed -e 's|^.*/||' -e 's|\..*||' |\
    head -n $UMI_SUBSET_N > sample$UMI_SUBSET_N.txt
fi

# Consensus
CON_DIR=$OUT_DIR/racon
longread_umi consensus_racon \
  $UMI_DIR/read_binning/bins    `# Path to UMI bins`\
  $CON_DIR                      `# Output folder`\
  4                             `# Number of racon polishing times`\
  $THREADS                      `# Number of threads`\
  sample$UMI_SUBSET_N.txt       `# List of bins to process`

# Polishing
POLISH_DIR1=$OUT_DIR/racon_medaka
longread_umi polish_medaka \
  $CON_DIR/consensus_*.fa       `# Path to consensus data`\
  $MAX_LENGTH                   `# Sensible chunk size`\
  $UMI_DIR/read_binning/bins    `# Path to UMI bins`\
  $POLISH_DIR1                  `# Output folder`\
  $THREADS                      `# Number of threads`\
  sample$UMI_SUBSET_N.txt       `# List of bins to process` \
  $TURBO                        `# Uses ALL threads with medaka`

POLISH_DIR2=$OUT_DIR/racon_medaka_medaka
longread_umi polish_medaka \
  $POLISH_DIR1/consensus_*.fa   `# Path to consensus data`\
  $MAX_LENGTH                   `# Sensible chunk size`\
  $UMI_DIR/read_binning/bins    `# Path to UMI bins`\
  $POLISH_DIR2                  `# Output folder`\
  $THREADS                      `# Number of threads`\
  sample$UMI_SUBSET_N.txt       `# List of bins to process` \
  $TURBO                        `# Uses ALL threads with medaka`

# Trim UMI consensus data
longread_umi trim_amplicon \
  $POLISH_DIR2         `# Path to consensus data`\
  '"consensus*fa"'    `# Consensus file pattern. Regex must be flanked by '"..."'`\
  $OUT_DIR             `# Output folder`\
  $FW2                 `# Forward primer sequence`\
  $RV2                 `# Reverse primer sequence`\
  $MIN_LENGTH          `# Min read length`\
  $MAX_LENGTH          `# Max read length` \
  $THREADS             `# Number of threads` \
  $LOG_DIR

# Generate variants

## Subset to UMI consensus sequences with min read coverage
POLISH_DIR2_NAME=${POLISH_DIR2##*/}
POLISH_DIR2_NAME=${POLISH_DIR2_NAME%.*}
  
$GAWK -v UBS="$UMI_COVERAGE_MIN" '
  /^>/{
    match($0,/;ubs=([0-9]+)/, s)
    if (s[1]+0 >= UBS+0){
      print
      getline
      print
    }
  }
' $OUT_DIR/consensus_${POLISH_DIR2_NAME}.fa \
> $OUT_DIR/consensus_${POLISH_DIR2_NAME}_${UMI_COVERAGE_MIN}.fa

## Variant calling of from UMI consensus sequences
longread_umi variants \
  $OUT_DIR/consensus_${POLISH_DIR2_NAME}_${UMI_COVERAGE_MIN}.fa `# Path to consensus data`\
  $OUT_DIR/variants `# Output folder`\
  $THREADS `# Number of threads`

## Copy variants
cp $OUT_DIR/variants/variants.fa $OUT_DIR

## Testing
exit 0
THREADS=60
INPUT_READS=reads.fq
UMI_SUBSET_N=1000000
UMI_COVERAGE_MIN=30
