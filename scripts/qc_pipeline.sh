#!/bin/bash
# DESCRIPTION
#    Script for performing quality control of longread UMI
#    consensus data.
#    
# IMPLEMENTATION
#    author   Søren Karst (sorenkarst@gmail.com)
#             Ryan Ziels (ziels@mail.ubc.ca)
#    license  GNU General Public License

### Source commands and subscripts -------------------------------------
. $LONGREAD_UMI_PATH/scripts/dependencies.sh # Path to dependencies script

### Description ----------------------------------------------------------------

USAGE="
-- longread_umi qc_pipeline: UMI consensus data statistics and compare to references
   Calculates data statistics (UMI bins size, UMI cluster size, yield, read length etc.).
   Mapping of read data and UMI consensus sequences to reference sequences to allow for 
   error profiling. Detects chimeras using uchime2_ref. Detects contamination by
   comparing mapping results to known references and the SILVA database - only works
   if reference database contains all expected sequences. Alternatively, use variants.fa
   as reference database.
   
usage: $(basename "$0" .sh) [-h] (-d files -c files -r files -s file -u dir -o dir -t value) 

where:
    -h  Show this help text.
    -d  List of read files seperated by ';'
        i.e. 'reads.fq;trim/reads_tf.fq'
        First read file used for read classification. 'reads_tf.fq' recommended.
    -c  List of consensus files seperated by ';'
        i.e. 'consensus_medaka_medaka.fa;racon/consensus_racon.fa'
        First consensus file used for comparison to alternative refs
        and for chimera checking. Subsequent consensus sequences only mapped to
        first reference.
    -r  List of reference files seperated by ';'.
        First reference is used for all mappings. Subsequent references
        only used for mapping first consensus file.
        'zymo_curated' refers to:
        "$REF_CURATED"
        'zymo_vendor' refers to:
        "$REF_VENDOR"
    -s  SILVA reference database in fasta format used for detecting contamination.
    -u  UMI consensus output folder.
    -o  Output folder. Default 'qc'.
    -t  Number of threads to use.

Example of SILVA database download:
wget https://www.arb-silva.de/fileadmin/silva_databases/
release_132/Exports/SILVA_132_SSURef_Nr99_tax_silva.fasta.gz
gunzip SILVA_132_SSURef_Nr99_tax_silva.fasta.gz
"

### Terminal Arguments ---------------------------------------------------------

# Import user arguments
while getopts ':hzd:c:r:s:u:o:t:' OPTION; do
  case $OPTION in
    h) echo "$USAGE"; exit 1;;
    d) READ_LIST=$OPTARG;;
    c) CON_LIST=$OPTARG;;
    r) REF_LIST=$OPTARG;;
    s) SILVA=$OPTARG;;
    u) UMI_DIR=$OPTARG;;
    o) OUT=$OPTARG;;
    t) THREADS=$OPTARG;;
    :) printf "missing argument for -$OPTARG\n" >&2; exit 1;;
    \?) printf "invalid option for -$OPTARG\n" >&2; exit 1;;
  esac
done

# Check missing arguments
MISSING="is missing but required. Exiting."
if [ -z ${CON_LIST+x} ]; then echo "-c $MISSING"; echo ""; echo "$USAGE"; exit 1; fi; 
if [ -z ${REF_LIST+x} ]; then echo "-r $MISSING"; echo ""; echo "$USAGE"; exit 1; fi;
if [ -z ${OUT+x} ]; then
  echo "-o is missing. Using default output dir 'qc'."
  OUT=qc
fi
if [ -z ${THREADS+x} ]; then echo "-t is missing. Defaulting to 1 thread."; THREADS=1; fi;

### Data processing -----------------------------------------------------

# Prepare folders and vars
mkdir $OUT

REF_LIST=$(echo $REF_LIST |\
  sed -e 's/[,;\t]/ /g' \
      -e "s|zymo_curated|$REF_CURATED|g" \
      -e "s|zymo_vendor|$REF_VENDOR|g")
REF1=$(echo "$REF_LIST " | cut -d' ' -f1)
REF_=$(echo "$REF_LIST " | cut -d' ' -f2-)

CON_LIST=$(echo $CON_LIST | sed -e 's/[,;\t]/ /g')
CON1=$(echo "$CON_LIST " | cut -d' ' -f1)
CON_=$(echo "$CON_LIST " | cut -d' ' -f2-)
CON1_NAME=${CON1%.*}
CON1_NAME=${CON1_NAME##*/}

# Prepare binning statistics
if [ ! -z ${UMI_DIR+x} ]; then
 cp $UMI_DIR/umi_binning/read_binning/umi_bin_map.txt $OUT/
 cp $UMI_DIR/umi_binning/read_binning/umi_binning_stats.txt $OUT/
fi

# Process read data
if [ ! -z ${READ_LIST+x} ]; then
  # Format read list
  READ_LIST=$(echo $READ_LIST | sed -e 's|[,;\t]| |g')
  READ1=$(echo "$READ_LIST " | cut -d' ' -f1)
  READ_=$(echo "$READ_LIST " | cut -d' ' -f2-)
  echo "data_type,read_count,bp_total,bp_average" > $OUT/data_stats.txt
  
  for READ_FILE in $READ_LIST; do
    # Format name
    READ_NAME=${READ_FILE%.*};
    READ_NAME=${READ_NAME##*/};
  
    # Subset data
    $SEQTK sample \
      -s1334 \
      $READ_FILE \
      5000 |\
      $SEQTK seq -a - \
      > $OUT/${READ_NAME}.fa
  
    # Read stats
    fastq_stats(){
      $GAWK -v sample="$2" '
        NR%4==2{
          rc++
          bp+=length
        } END {
          print sample","rc","bp","bp/rc
        }
      ' $1
    }
    fastq_stats $READ_FILE $READ_NAME >> $OUT/data_stats.txt
 done
 
 $MINIMAP2 \
   -x map-ont \
   $REF1 \
   $READ1 \
   -t $THREADS |\
   $GAWK '$13 ~ /tp:A:P/{split($6,tax,"_"); print $1, tax[1]"_"tax[2]}'\
   > $OUT/read_classification.txt
fi


# Prepare consensus data
cp -t $OUT $CON_LIST

# Mapping
for DATA_FILE in $OUT/*.fa; do
  DATA_NAME=${DATA_FILE%.*};
  DATA_NAME=${DATA_NAME##*/};
  $MINIMAP2 \
    -ax map-ont \
    $REF1 \
    $DATA_FILE \
    -t $THREADS --cs |\
    $SAMTOOLS view -F 2308 - |\
    cut -f1-9,12,21 > $OUT/$DATA_NAME.sam
done

for REF in $REF_; do
  REF_NAME=${REF%.*};
  REF_NAME=${REF_NAME##*/};
  $MINIMAP2 -ax map-ont \
    $REF \
    $CON1 \
    -t $THREADS --cs |\
    $SAMTOOLS view -F 2308 - |\
    cut -f1-9,12,21 \
    > $OUT/${CON1_NAME}_${REF_NAME}.sam
done

# Copy refs
cp -t $OUT $REF_LIST

# Detect chimeras
$USEARCH \
  -uchime2_ref $CON1 \
  -uchimeout $OUT/${CON1_NAME}_chimera.txt \
  -db $REF1 \
  -strand plus \
  -mode sensitive

# Detect contamination
if [ ! -z ${SILVA+x} ]; then
  SILVA_NAME=${SILVA%.*};
  SILVA_NAME=${SILVA_NAME##*/};
  REF1_NAME=${REF1%.*};
  REF1_NAME=${REF1_NAME##*/};

  # Extract SSU from data
  $CUTADAPT \
    --discard-untrimmed \
    -m 1200 \
    -M 2000 \
    -a TGYACWCACCGCCCGTC \
    $CON1 > $OUT/${CON1_NAME}_ssu.fa

  # Map to SILVA database
  $MINIMAP2 \
    -ax map-ont \
    $SILVA \
    $OUT/${CON1_NAME}_ssu.fa \
    -t $THREADS --cs |\
    $SAMTOOLS view -F 2308 - |\
    cut -f1-9,12,21 \
    > $OUT/${CON1_NAME}_ssu_${SILVA_NAME}.sam
	

  # Map to ref database
  $MINIMAP2 \
    -ax map-ont \
    $REF1 \
    $OUT/${CON1_NAME}_ssu.fa \
    -t $THREADS --cs |\
    $SAMTOOLS view -F 2308 - |\
    cut -f1-9,12,21 \
    > $OUT/${CON1_NAME}_ssu_${REF1_NAME}.sam
  
  # Export target taxonomy
  $GAWK \
    -v SAM="$OUT/${CON1_NAME}_ssu_${SILVA_NAME}.sam" \
    -v SILVA="$SILVA" \
    '
      (FILENAME == SAM){
        TARGETS[$3]++ 
      }
      (FILENAME == SILVA && $0 ~ /^>/){
        TARGET=substr($1,2)
        if (TARGET in TARGETS){
          print TARGET, $2
        }
      }
    ' \
    $OUT/${CON1_NAME}_ssu_${SILVA_NAME}.sam \
    $SILVA > $OUT/${CON1_NAME}_ssu_${SILVA_NAME}_tax.txt
fi

# Testing
exit 0
longread_umi qc_pipeline \
  -d "test_reads.fq;umi_binning/trim/reads_tf.fq" \
  -c "consensus_racon_medaka_medaka.fa;variants_all.fa" \
  -r "zymo_curated;zymo_vendor;variants_all.fa" \
  -s SILVA_132_SSURef_Nr99_tax_silva.fasta
