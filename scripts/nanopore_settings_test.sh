#!/bin/bash
# DESCRIPTION
#    Testing settings for racon consensus and medaka polishing.
#    
# IMPLEMENTATION
#    author   Søren Karst (sorenkarst@gmail.com)
#             Ryan Ziels (ziels@mail.ubc.ca)
#    license  GNU General Public License
#
# To-do:
# - Fix logging
# - Broken due to missing settings for UMI binning

### Description ----------------------------------------------------------------

USAGE="
-- longread_umi nanopore_settings_test: Test impact of polishing rounds on UMI consensus.

usage: $(basename "$0" .sh) [-h] (-d file -n value -c value -o dir -s value -e value) 
(-m value -M value -f string -F string -r string -R string -t value -T value) 
(-x value -y value -q string -p -u dir ) 

where:
    -h  Show this help text.
    -d  Single file containing raw Nanopore data in fastq format.
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
    -T  Number of medaka jobs to start. Threads pr. job is threads/jobs.
        [Default = 1].
    -x  Test Racon consensus rounds from 1 to <value>.
    -y  Test Medaka polishing rounds from 1 to <value>.
    -q  Medaka model used for polishing. r941_min_high, r10_min_high etc.
    -p  Flag to disable Nanopore trimming and filtering.
    -u  Directory with UMI binned reads.

Test run:
longread_umi nanopore_settings_test 
  -d test_reads.fq 
  -o settings_test 
  -w rrna_operon 
  -t 100 
  -T 20 
  -x 4 
  -y 3 
  -n 1000
"

### Terminal Arguments ---------------------------------------------------------

# Import user arguments
while getopts ':hzd:o:s:e:m:M:f:F:r:R:n:w:t:T:x:y:q:pu:' OPTION; do
  case $OPTION in
    h) echo "$USAGE"; exit 1;;
    d) INPUT_READS=$OPTARG;;
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
    T) MEDAKA_JOBS=$OPTARG;;
    x) RACON_ROUNDS=$OPTARG;;
    y) MEDAKA_ROUNDS=$OPTARG;;
    q) MEDAKA_MODEL=$OPTARG;;
    p) TRIM_FLAG="-p";;
    u) UMI_DIR=$OPTARG;;
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
if [ -z ${MEDAKA_JOBS+x} ]; then echo "-T is missing. Medaka jobs set to 1."; MEDAKA_JOBS=1; fi;
if [ -z ${RACON_ROUNDS+x} ]; then echo "-x $MISSING"; echo "$USAGE"; exit 1; fi;
if [ -z ${MEDAKA_ROUNDS+x} ]; then echo "-y $MISSING"; echo "$USAGE"; exit 1; fi;
if [ -z ${MEDAKA_MODEL+x} ]; then echo "-q $MISSING"; echo "$USAGE"; exit 1; fi;

### Source commands and subscripts -------------------------------------
. $LONGREAD_UMI_PATH/scripts/dependencies.sh # Path to dependencies script
mkdir $OUT_DIR

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
echo "Medaka jobs: $MEDAKA_JOBS"
echo "Racon rounds: $RACON_ROUNDS"
echo "Medaka rounds: $MEDAKA_ROUNDS"
echo "Medaka model: $MEDAKA_MODEL"
echo "UMI directory: $UMI_DIR"
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
    -t $THREADS          `# Number of threads` \
    $TRIM_FLAG
fi

# Sample UMI bins for testing
if [ ! -z ${UMI_SUBSET_N+x} ]; then
  find  $UMI_DIR -name 'umi*bins.fastq' |\
    sed -e 's|^.*/||' -e 's|\..*||' |\
    head -n $UMI_SUBSET_N > $OUT_DIR/sample$UMI_SUBSET_N.txt
fi

# Consensus and polishing
for i in `seq 1 $RACON_ROUNDS`; do
  # Racon consensus
  CON_NAME=raconx$i
  longread_umi consensus_racon \
    -d $UMI_DIR/read_binning/bins          `# Path to UMI bins`\
    -o $OUT_DIR/$CON_NAME                   `# Output folder`\
    -r $i                                  `# Number of racon polishing times`\
    -t $THREADS                            `# Number of threads`\
    -n $OUT_DIR/sample$UMI_SUBSET_N.txt    `# List of bins to process`
  CON=$OUT_DIR/$CON_NAME/consensus_${CON_NAME}.fa
  # Medaka polishing
  for j in `seq 1 $MEDAKA_ROUNDS`; do
    POLISH_NAME=${CON_NAME}_medakax$j
    longread_umi polish_medaka \
      -c $CON                              `# Path to consensus data`\
      -m $MEDAKA_MODEL                     `# Path to consensus data`\
      -l $MAX_LENGTH                       `# Sensible chunk size`\
      -d $UMI_DIR                          `# Path to UMI bins`\
      -o $OUT_DIR/${POLISH_NAME}   `# Output folder`\
      -t $THREADS                          `# Number of threads`\
      -n $OUT_DIR/sample$UMI_SUBSET_N.txt  `# List of bins to process` \
      -T $MEDAKA_JOBS                      `# Uses ALL threads with medaka`
    CON=$OUT_DIR/${POLISH_NAME}/consensus_${POLISH_NAME}.fa
  done
done

# Trim UMI consensus data
longread_umi trim_amplicon \
  -d '$OUT_DIR/racon*'      `# Path to consensus data`\
  -p '"consensus*fa"'     `# Consensus file pattern. Regex must be flanked by '"..."'`\
  -o $OUT_DIR             `# Output folder`\
  -F $FW2                 `# Forward primer sequence`\
  -R $RV2                 `# Reverse primer sequence`\
  -m $MIN_LENGTH          `# Min read length`\
  -M $MAX_LENGTH          `# Max read length` \
  -t $THREADS             `# Number of threads` \
  -l $LOG_DIR

# Perform qc
CON_LIST=$(echo $OUT_DIR/consensus* | sed 's/ /;/g')
longread_umi qc_pipeline \
  -d $UMI_DIR/trim/reads_tf.fq \
  -c $CON_LIST \
  -u $UMI_DIR \
  -o $OUT_DIR/qc \
  -r "zymo_curated" \
  -t $THREADS

