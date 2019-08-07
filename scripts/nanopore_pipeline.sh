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

USAGE="$(basename "$0" .sh) [-h] [-d file -s value -c value -f value -r value -t value] 
-- longread_umi nanopore_pipeline: Generates UMI consensus sequences from
   raw Nanopore fastq reads with UMIs in both terminals.

where:
    -h  Show this help text.
    -d  Single file containing raw Nanopore data in fastq format.
    -s  Process <value> number of bins.
    -c  Minimum read coverage for using UMI consensus sequences for 
        variant calling.
    -f  Check start of read up to f bp for UMIs. Default: 70.
    -r  Check end of read up to r bp for UMIs. Default: 80.
    -t  Number of threads to use. Ps. Medaka bug means medaka uses all available..
"

### Terminal Arguments ---------------------------------------------------------

# Import user arguments
while getopts ':hzd:s:c:f:r:t:' OPTION; do
  case $OPTION in
    h) echo "$USAGE"; exit 1;;
    d) INPUT_READS=$OPTARG;;
    s) UMI_SUBSET_N=$OPTARG;;
    c) UMI_COVERAGE_MIN=$OPTARG;;
    f) START_READ_CHECK=$OPTARG;;
    r) END_READ_CHECK=$OPTARG;;
    t) THREADS=$OPTARG;;
    :) printf "missing argument for -$OPTARG\n" >&2; exit 1;;
    \?) printf "invalid option for -$OPTARG\n" >&2; exit 1;;
  esac
done

# Check missing arguments
MISSING="is missing but required. Exiting."
if [ -z ${INPUT_READS+x} ]; then echo "-d $MISSING"; echo "$USAGE"; exit 1; fi; 
if [ -z ${UMI_SUBSET_N+x} ]; then echo "-s $MISSING"; echo "$USAGE"; exit 1; fi; 
if [ -z ${UMI_COVERAGE_MIN+x} ]; then echo "-c $MISSING"; echo "$USAGE"; exit 1; fi;
if [ -z ${START_READ_CHECK+x} ]; then START_READ_CHECK=70; fi; 
if [ -z ${END_READ_CHECK+x} ]; then END_READ_CHECK=80; fi;
if [ -z ${THREADS+x} ]; then echo "-t is missing. Defaulting to 1 thread."; THREADS=1; fi;

### Source commands and subscripts -------------------------------------
. $LONGREAD_UMI_PATH/scripts/dependencies.sh # Path to dependencies script

### Pipeline -----------------------------------------------------------
# Logging
LOG_NAME="longread_umi_nanopore_pipeline_log_$(date +"%Y-%m-%d-%T").txt"
echo "longread_umi nanopore_pipeline log" >> $LOG_NAME
longread_umi_version_dump $LOG_NAME
exec &> >(tee -a "$LOG_NAME")
exec 2>&1
echo ""
echo "### Settings:"
echo "Threads: $THREADS"
echo "UMI subsampling: $UMI_SUBSET_N"
echo "Bin size cutoff: $UMI_COVERAGE_MIN"
echo ""

# Read filtering and UMI binning
longread_umi umi_binning  \
  $INPUT_READS `# Raw nanopore data in fastq format`\
  umi_binning  `# Output folder`\
  $THREADS     `# Number of threads`\
  3500         `# Min read length`\
  6000         `# Max read length` \
  $START_READ_CHECK `# Start of read to check` \
  $END_READ_CHECK `# End of read to check`

# Sample UMI bins for testing
find umi_binning/read_binning/bins \
   -name 'umi*bins.fastq' | sed -e 's|^.*/||' -e 's|\..*||' |\
   head -n $UMI_SUBSET_N> sample$UMI_SUBSET_N.txt

# Consensus
longread_umi consensus_racon \
  umi_binning/read_binning/bins `# Path to UMI bins`\
  racon                         `# Output folder`\
  $THREADS                      `# Number of threads`\
  sample$UMI_SUBSET_N.txt       `# List of bins to process`

# Polishing
longread_umi polish_medaka \
  racon/consensus_*.fa          `# Path to consensus data`\
  umi_binning/read_binning/bins `# Path to UMI bins`\
  racon_medaka                  `# Output folder`\
  $THREADS                      `# Number of threads`\
  sample$UMI_SUBSET_N.txt       `# List of bins to process`

longread_umi polish_medaka \
  racon_medaka/consensus_*.fa   `# Path to consensus data`\
  umi_binning/read_binning/bins `# Path to UMI bins`\
  racon_medaka_medaka           `# Output folder`\
  $THREADS                      `# Number of threads`\
  sample$UMI_SUBSET_N.txt       `# List of bins to process`

# Trim UMI consensus data
longread_umi trim_amplicon \
  racon_medaka_medaka  `# Path to consensus data`\
  consensus*fa         `# Consensus file pattern`\
  .                    `# Output folder`\
  rrna_8f2490r         `# Primers used`\
  $THREADS             `# Number of threads`

# Generate variants

## Subset to UMI consensus sequences with min read coverage
$GAWK -v bsco="$UMI_COVERAGE_MIN" '
  /^>/{
    s=$0
    gsub(".*size=", "", s)
    if (s+0 >= bsco+0){
      print
      getline
      print
    }
  }
' consensus_racon_medaka_medaka.fa \
> consensus_racon_medaka_medaka_${UMI_COVERAGE_MIN}.fa

## Variant calling of from UMI consensus sequences
longread_umi variants \
  consensus_racon_medaka_medaka_${UMI_COVERAGE_MIN}.fa `# Path to consensus data`\
  variants `# Output folder`\
  $THREADS `# Number of threads`

## Copy phase
cp variants/variants*fa .


## Testing
exit 0
THREADS=60
INPUT_READS=reads.fq
UMI_SUBSET_N=1000000
UMI_COVERAGE_MIN=30
