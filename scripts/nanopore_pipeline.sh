#!/bin/bash
# DESCRIPTION
#    longread_umi nanopore_pipeline script. 
#    
# IMPLEMENTATION
#    author   Søren Karst (sorenkarst@gmail.com)
#             Ryan Ziels (ziels@mail.ubc.ca)
#    license  GNU General Public License
#
# To-do:
# - Fix logging

### Description ----------------------------------------------------------------

USAGE="
-- longread_umi nanopore_pipeline: Generate UMI consensus sequences from Nanopore data
   
usage: $(basename "$0" .sh) [-h] [ -k flag] (-d file -v value -o dir -s value) 
(-e value -m value -M value -f string -F string -r string -R string )
( -c value -p value -n value -u dir -U string -t value -T value ) 

where:
    -h  Show this help text.
    -d  Single file containing raw Nanopore data in fastq format.
    -v  Minimum read coverage for using UMI consensus sequences for 
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
    -c  Number of iterative rounds of consensus calling with Racon.
    -p  Number of iterative rounds of consensus calling with Medaka.
    -q  Medaka model used for polishing. r941_min_high_g360, r103_min_high_g360 etc.
    -k  Flag for keeping failed bins in output.
    -n  Process n number of bins. If not defined all bins are processed.
        Pratical for testing large datasets.
    -u  Directory with UMI binned reads.
    -U  UMI filter settings. Define settings for:
        - UMI match error mean (UMEM): Mean match error between reads in a bin
          and the UMI reference.
        - UMI match error SD (UMESD): Standard deviation for match error between
          reads in a bin and the UMI reference.
        - Bin cluster ratio (BCR): Ratio between UMI bin size and UMI cluster size.
        - Read orientation ratio (ROR): n(+ strand reads)/n(all reads). '0' is the
          means disabled.
        Settings can be provided as a string: 'UMEM/UMESD/BCR/ROR'
        Or as a preset:
        - 'r941_min_high_g360' == '3;2;6;0.3'
        - 'r103_min_high_g360' == '3;2.5;12;0.3'
    -t  Number of threads to use.
    -T  Number of medaka jobs to start. Threads pr. job is threads/jobs.
        [Default = 1].
"

### Terminal Arguments ---------------------------------------------------------

# Import user arguments
while getopts ':hzd:v:o:s:e:m:M:f:F:r:R:c:p:q:kn:u:U:t:T:' OPTION; do
  case $OPTION in
    h) echo "$USAGE"; exit 1;;
    d) INPUT_READS=$OPTARG;;
    v) UMI_COVERAGE_MIN=$OPTARG;;
    o) OUT_DIR=$OPTARG;;
    s) START_READ_CHECK=$OPTARG;;
    e) END_READ_CHECK=$OPTARG;;
    m) MIN_LENGTH=$OPTARG;;
    M) MAX_LENGTH=$OPTARG;;
    f) FW1=$OPTARG;;
    F) FW2=$OPTARG;;
    r) RV1=$OPTARG;;
    R) RV2=$OPTARG;;  
    c) CON_N=$OPTARG;;
    p) POL_N=$OPTARG;;
    q) MEDAKA_MODEL=$OPTARG;;
    k) KEEP="YES";;
    n) UMI_SUBSET_N=$OPTARG;;
    u) UMI_DIR=$OPTARG;;
    U) UMI_FILTER_SETTINGS=$OPTARG;;
    t) THREADS=$OPTARG;;
    T) MEDAKA_JOBS=$OPTARG;;
    :) printf "missing argument for -$OPTARG\n" >&2; exit 1;;
    \?) printf "invalid option for -$OPTARG\n" >&2; exit 1;;
  esac
done

# Check missing arguments
MISSING="is missing but required. Exiting."

if [ -z ${UMI_FILTER_SETTINGS+x} ]; then 
  echo "-U $MISSING"
  echo "$USAGE"
  exit 1
elif [ "$UMI_FILTER_SETTINGS" == "r941_min_high_g360" ]; then
  UMI_MATCH_ERROR=3
  UMI_MATCH_ERROR_SD=2
  BIN_CLUSTER_RATIO=6
  RO_FRAC=0.3
elif [ "$UMI_FILTER_SETTINGS" == "r103_min_high_g360" ]; then
  UMI_MATCH_ERROR=3
  UMI_MATCH_ERROR_SD=2.5
  BIN_CLUSTER_RATIO=12
  RO_FRAC=0.3
else
  ufs=(`echo $UMI_FILTER_SETTINGS | cut -d ";"  --output-delimiter=" " -f 1-`)
  UMI_MATCH_ERROR=${ufs[0]}
  UMI_MATCH_ERROR_SD=${ufs[1]}
  BIN_CLUSTER_RATIO=${ufs[2]}
  RO_FRAC=${ufs[3]}
  if [[ -z $UMI_MATCH_ERROR || -z $UMI_MATCH_ERROR_SD || -z $BIN_CLUSTER_RATIO || -z $RO_FRAC ]]; then
    echo "One or more filter settings is unset:"
    echo "UMI match error mean: $UMI_MATCH_ERROR"
    echo "UMI match error SD: $UMI_MATCH_ERROR_SD"
    echo "Bin cluster ratio: $BIN_CLUSTER_RATIO"
    echo "Read orientation ratio: $RO_FRAC"
    echo "Exiting ..."
    exit 1
  fi
fi

if [ -z ${INPUT_READS+x} ]; then echo "-d $MISSING"; echo "$USAGE"; exit 1; fi; 
if [ -z ${UMI_COVERAGE_MIN+x} ]; then echo "-v $MISSING"; echo "$USAGE"; exit 1; fi;
if [ -z ${OUT_DIR+x} ]; then echo "-o $MISSING"; echo "$USAGE"; exit 1; fi;
if [ -z ${START_READ_CHECK+x} ]; then echo "-s $MISSING"; echo "$USAGE"; exit 1; fi;
if [ -z ${END_READ_CHECK+x} ]; then echo "-e $MISSING"; echo "$USAGE"; exit 1; fi;
if [ -z ${MIN_LENGTH+x} ]; then echo "-m $MISSING"; echo "$USAGE"; exit 1; fi;
if [ -z ${MAX_LENGTH+x} ]; then echo "-M $MISSING"; echo "$USAGE"; exit 1; fi;
if [ -z ${FW1+x} ]; then echo "-f $MISSING"; echo "$USAGE"; exit 1; fi;
if [ -z ${FW2+x} ]; then echo "-F $MISSING"; echo "$USAGE"; exit 1; fi;
if [ -z ${RV1+x} ]; then echo "-r $MISSING"; echo "$USAGE"; exit 1; fi;
if [ -z ${RV2+x} ]; then echo "-R $MISSING"; echo "$USAGE"; exit 1; fi;
if [ -z ${CON_N+x} ]; then echo "-c $MISSING"; echo "$USAGE"; exit 1; fi;
if [ -z ${POL_N+x} ]; then echo "-p $MISSING"; echo "$USAGE"; exit 1; fi;
if [ -z ${MEDAKA_MODEL+x} ]; then echo "-q $MISSING"; echo "$USAGE"; exit 1; fi;
if [ -z ${KEEP+x} ]; then KEEP="NO"; fi;
if [ -z ${THREADS+x} ]; then echo "-t is missing. Defaulting to 1 thread."; THREADS=1; fi;
if [ -z ${MEDAKA_JOBS+x} ]; then echo "-T is missing. Medaka jobs set to 1."; MEDAKA_JOBS=1; fi;

### Source commands and subscripts -------------------------------------
. $LONGREAD_UMI_PATH/scripts/dependencies.sh # Path to dependencies script

if [ -d $OUT_DIR ]; then
  echo ""
  echo "$OUT_DIR exists. Remove existing directory or rename desired output directory."
  echo "Analysis aborted ..."
  echo ""
  exit 1 
else
  mkdir $OUT_DIR
fi

### Pipeline -----------------------------------------------------------
# Logging
LOG_DIR=$OUT_DIR/logs
mkdir $LOG_DIR

LOG_NAME="$LOG_DIR/longread_umi_nanopore_pipeline_log_$(date +"%Y-%m-%d-%T").txt"
echo "longread_umi nanopore_pipeline log" >> $LOG_NAME
longread_umi_version_dump $LOG_NAME
exec &> >(tee -a "$LOG_NAME")
exec 2>&1
echo ""
echo "### Settings:"
echo "Input reads: $INPUT_READS"
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
echo "Racon consensus rounds: $CON_N"
echo "Medaka consensus rounds: $POL_N"
echo "Medaka model: $MEDAKA_MODEL"
echo "UMI filter settings: $UMI_FILTER_SETTINGS"
echo "UMI match error mean: $UMI_MATCH_ERROR"
echo "UMI match error SD: $UMI_MATCH_ERROR_SD"
echo "Bin cluster ratio: $BIN_CLUSTER_RATIO"
echo "Read orientation ratio: $RO_FRAC"
echo "Keep failed bins: $KEEP"
echo "Bin size cutoff: $UMI_COVERAGE_MIN"
echo "UMI binning dir: $UMI_DIR"
echo "Threads: $THREADS"
echo "Medaka jobs: $MEDAKA_JOBS"
echo ""

# Read filtering and UMI binning
if [ -z ${UMI_DIR+x} ]; then
  UMI_DIR=$OUT_DIR/umi_binning
  longread_umi umi_binning  \
    -d $INPUT_READS      `# Raw nanopore data in fastq format`\
    -o $UMI_DIR          `# Output folder`\
    -m $MIN_LENGTH       `# Min read length`\
    -M $MAX_LENGTH       `# Max read length` \
    -s $START_READ_CHECK `# Start of read to check` \
    -e $END_READ_CHECK   `# End of read to check` \
    -f $FW1              `# Forward adaptor sequence` \
    -F $FW2              `# Forward primer sequence` \
    -r $RV1              `# Reverse adaptor sequence` \
    -R $RV2              `# Reverse primer sequence` \
    -u $UMI_MATCH_ERROR  `# UMI match error filter` \
    -U $UMI_MATCH_ERROR_SD`# UMI match error SD filter` \
    -O $RO_FRAC          `# Min read orientation fraction` \
    -S $BIN_CLUSTER_RATIO`# UMI bin/cluster ratio cutoff` \
    -N 10000             `# Maximum number of reads +/-` \
    -t $THREADS          `# Number of threads`
fi

# Sample UMI bins for testing
$GAWK \
  -v SEED="$RANDOM"  \
  -v SUB_N="$UMI_SUBSET_N" \
  -v KEEP="$KEEP" \
  '
  # Handle failed UMI bins
  $0 !~ /fail/ && KEEP == "NO" && NR>1{
    umi[NR]=$1
  }
  KEEP == "YES"  && NR>1 {
    umi[NR]=$1
  }
  END {
    # Handle UMI bin subsetting
    if (UMI_SUBSET_N+0 > 0){
      # Determine bins to subset
      UMI_N = length(umi)
      if (SUB_N > UMI_N){
        SUB_N = UMI_N
      }
      # Sample bins
      srand(SEED)
      while (SAMPLED <= SUB_N){
        # Sample line
        SAMPLED++
        LINE=int((length(umi)+1) * rand())
        UMI_NAME=umi[LINE]
        sub(";.*", "bins.fastq", UMI_NAME)
        print UMI_NAME
        # Remove sampled UMI
        delete umi[LINE]
      }
    } else {
      for (i in umi){
        UMI_NAME=umi[i]
        sub(";.*", "", UMI_NAME)
        print UMI_NAME "bins"
        print UMI_NAME "_"
      }
    }
  }
  ' \
  $UMI_DIR/read_binning/umi_binning_stats.txt |\
  sort -u \
  > $OUT_DIR/processed_bins.txt
  

# Consensus
CON_NAME=raconx${CON_N}
CON_DIR=$OUT_DIR/$CON_NAME
longread_umi consensus_racon \
  -d $UMI_DIR/read_binning/bins           `# Path to UMI bins`\
  -o ${CON_DIR}                           `# Output folder`\
  -p map-ont                              `# Minimap preset`\
  -a "--no-trimming"                      `# Extra args for racon`\
  -r $CON_N                               `# Number of racon polishing times`\
  -t $THREADS                             `# Number of threads`\
  -n $OUT_DIR/processed_bins.txt    `# List of bins to process`

# Polishing
CON=${CON_DIR}/consensus_${CON_NAME}.fa
for j in `seq 1 $POL_N`; do
  POLISH_NAME=medakax${j}
  POLISH_DIR=${CON_DIR}_${POLISH_NAME}
  longread_umi polish_medaka \
    -c $CON                              `# Path to consensus data`\
    -m $MEDAKA_MODEL                     `# Path to consensus data`\
    -l $MAX_LENGTH                       `# Sensible chunk size`\
    -d $UMI_DIR                          `# Path to UMI bins`\
    -o $POLISH_DIR                       `# Output folder`\
    -t $THREADS                          `# Number of threads`\
    -n $OUT_DIR/processed_bins.txt       `# List of bins to process`\
    -T $MEDAKA_JOBS                      `# Uses ALL threads with medaka`
  CON=$POLISH_DIR/consensus_${CON_NAME}_${POLISH_NAME}.fa
done
  

# Trim UMI consensus data
longread_umi trim_amplicon \
  -d $POLISH_DIR          `# Path to consensus data`\
  -p '"consensus*fa"'     `# Consensus file pattern. Regex must be flanked by '"..."'`\
  -o $OUT_DIR             `# Output folder`\
  -F $FW2                 `# Forward primer sequence`\
  -R $RV2                 `# Reverse primer sequence`\
  -m $MIN_LENGTH          `# Min read length`\
  -M $MAX_LENGTH          `# Max read length` \
  -t $THREADS             `# Number of threads` \
  -l $LOG_DIR

# Generate variants

## Subset to UMI consensus sequences with min read coverage
$GAWK -v UBS="$UMI_COVERAGE_MIN" '
  /^>/{
    match($0,/;ubs=([0-9]+)/, s)
    if (s[1]+0 >= UBS+0){
      print
      getline
      print
    }
  }
' $OUT_DIR/consensus_${CON_NAME}_${POLISH_NAME}.fa \
> $OUT_DIR/consensus_${CON_NAME}_${POLISH_NAME}_${UMI_COVERAGE_MIN}.fa

## Variant calling of from UMI consensus sequences
longread_umi variants \
  -c $OUT_DIR/consensus_${CON_NAME}_${POLISH_NAME}_${UMI_COVERAGE_MIN}.fa `# Path to consensus data`\
  -o $OUT_DIR/variants `# Output folder`\
  -t $THREADS `# Number of threads`

## Copy variants
cp $OUT_DIR/variants/variants.fa $OUT_DIR

