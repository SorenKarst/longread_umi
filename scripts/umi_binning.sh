#!/bin/bash

# DESCRIPTION
#    Script for binning long reads based on UMIs. Part of 
#    longread_umi.
#    
# IMPLEMENTATION
#    author   Søren Karst (sorenkarst@gmail.com)
#             Ryan Ziels (ziels@mail.ubc.ca)
#    license  GNU General Public License
#
# TO DO
#    Add terminal messages.
#    Optimize trimming and filtering for speed.
#    Add bin size limit
#    Add mapping against adaptors to remove UMI artifacts

### Source commands and subscripts -------------------------------------
. $LONGREAD_UMI_PATH/scripts/dependencies.sh # Path to dependencies script

### info ---------------------------------------------------------------
USAGE="
-- longread_umi umi_binning: Longread UMI detection and read binning.
   Tool requires UMIs in both ends of the read flanked by defined
   adaptor regions.

usage: $(basename "$0" .sh) [-h] (-d file -o dir -m value -M value )
(-s value -e value -f string -F string -r string -R string -p )
(-u value -U value -O value -S value -t value -T value) 

where:
    -h  Show this help text.
    -d  Reads in fastq format.
    -o  Output directory.
    -m  Minimum read length.
    -M  Maximum read length.
    -s  Check start of read up to s bp for UMIs.
    -e  Check end of read up to f bp for UMIs.
    -f  Forward adaptor sequence. 
    -F  Forward primer sequence.
    -r  Reverse adaptor sequence.
    -R  Reverse primer sequence.
    -p  Flag to disable Nanopore trimming and filtering.
        Use with PacBio reads.
    -u  Discard bins with a mean UMI match error above u.
    -U  Discard bins with a UMI match error standard
        deviation above U.
    -O  Normalize read orientation fraction to 'O' if n(-/+ strand)/n(total)
        is less than 'O'. Disable filter by setting -O to 0.
    -N  Max number of reads with +/- orientation. [Default = 10000]
    -S  UMI bin size/UMI cluster size cutoff.
    -t  Number of threads to use.
    -T  Number of threads to use for splitting reads into UMI bins. [Default is same as -t]
"

### Terminal Arguments ---------------------------------------------------------

# Import user arguments
while getopts ':hzd:o:m:M:s:e:f:F:r:R:pt:u:U:O:N:S:T:' OPTION; do
  case $OPTION in
    h) echo "$USAGE"; exit 1;;
    d) READ_IN=$OPTARG;;
    o) OUT_DIR=$OPTARG;;
    m) MIN_LENGTH=$OPTARG;;
    M) MAX_LENGTH=$OPTARG;;
    s) START_READ_CHECK=$OPTARG;;
    e) END_READ_CHECK=$OPTARG;;
    f) FW1=$OPTARG;;
    F) FW2=$OPTARG;;
    r) RV1=$OPTARG;;
    R) RV2=$OPTARG;;
    p) TRIM_FLAG=YES;;
    u) UMI_MATCH_ERROR=$OPTARG;;
    U) UMI_MATCH_ERROR_SD=$OPTARG;;
    O) RO_FRAC=$OPTARG;;
    N) MAX_BIN_SIZE=$OPTARG;;
    S) BIN_CLUSTER_RATIO=$OPTARG;;
    t) THREADS=$OPTARG;;
    T) BIN_THREADS=$OPTARG;;
    :) printf "missing argument for -$OPTARG\n" >&2; exit 1;;
    \?) printf "invalid option for -$OPTARG\n" >&2; exit 1;;
  esac
done

# Check missing arguments
MISSING="is missing but required. Exiting."
if [ -z ${READ_IN+x} ]; then echo "-d $MISSING"; echo ""; echo "$USAGE"; exit 1; fi;
if [ -z ${OUT_DIR+x} ]; then echo "-o $MISSING"; echo ""; echo "$USAGE"; exit 1; fi;
if [ -z ${MIN_LENGTH+x} ]; then echo "-m $MISSING"; echo ""; echo "$USAGE"; exit 1; fi; 
if [ -z ${MAX_LENGTH+x} ]; then echo "-M $MISSING"; echo ""; echo "$USAGE"; exit 1; fi;
if [ -z ${START_READ_CHECK+x} ]; then echo "-s $MISSING"; echo ""; echo "$USAGE"; exit 1; fi;
if [ -z ${END_READ_CHECK+x} ]; then echo "-e $MISSING"; echo ""; echo "$USAGE"; exit 1; fi;
if [ -z ${FW1+x} ]; then echo "-f $MISSING"; echo ""; echo "$USAGE"; exit 1; fi;
if [ -z ${FW2+x} ]; then echo "-F $MISSING"; echo ""; echo "$USAGE"; exit 1; fi;
if [ -z ${RV1+x} ]; then echo "-r $MISSING"; echo ""; echo "$USAGE"; exit 1; fi;
if [ -z ${RV2+x} ]; then echo "-R $MISSING"; echo ""; echo "$USAGE"; exit 1; fi;
if [ -z ${UMI_MATCH_ERROR+x} ]; then echo "-u $MISSING"; echo ""; echo "$USAGE"; exit 1; fi;
if [ -z ${UMI_MATCH_ERROR_SD+x} ]; then echo "-U $MISSING"; echo ""; echo "$USAGE"; exit 1; fi;
if [ -z ${RO_FRAC+x} ]; then echo "-O $MISSING"; echo ""; echo "$USAGE"; exit 1; fi;
if [ -z ${MAX_BIN_SIZE+x} ]; then echo "-N is missing. Defaulting to 10000 +/- reads ."; MAX_BIN_SIZE=10000; fi;
if [ -z ${BIN_CLUSTER_RATIO+x} ]; then echo "-S $MISSING"; echo ""; echo "$USAGE"; exit 1; fi;
if [ -z ${THREADS+x} ]; then echo "-t is missing. Defaulting to 1 thread."; THREADS=1; fi;
if [ -z ${BIN_THREADS+x} ]; then BIN_THREADS=$THREADS; fi;

### Primer formating
revcom() {
  echo $1 |\
  $GAWK '{print ">dummy\n" $0}' |\
  $SEQTK seq -r - |\
  $GAWK '!/^>/'  
}
FW1R=$(revcom "$FW1")
FW2R=$(revcom "$FW2")
RV1R=$(revcom "$RV1")
RV2R=$(revcom "$RV2")

### Read trimming and filtering -----------------------------------------------
mkdir $OUT_DIR
TRIM_DIR=$OUT_DIR/trim
mkdir $TRIM_DIR

# Trim data
if [ -z ${TRIM_FLAG+x} ]; then
  # Prepare porechop adapters

  # Define part of adapters.py to modify
  LEAD='^                    end_sequence=('\''SQK-NSK007_Y_Bottom'\'', '\''GCAATACGTAACTGAACGAAGT'\'')),$'
  TAIL='^def make_full_native_barcode_adapter(barcode_num):'
  
  # Format adapters - 12 bp substrings
  FW1_12=$(echo ${FW1:0:12})
  RV1_12=$(echo ${RV1:0:12})
  FW1R_12=$(revcom "$FW1_12")
  RV1R_12=$(revcom "$RV1_12")
  
  # Generate custom adapters.py
  ADAPTER_FMT="
            Adapter('LU_ADP_FWRV',
                    start_sequence=('lu_adpfwrv', '${FW1R_12}${RV1_12}'),
                    end_sequence=('lu_adpfwrv_rv', '${RV1R_12}${FW1_12}')),
            Adapter('LU_ADP_FWFW',
                    start_sequence=('lu_adpfwfw', '${FW1R_12}${FW1_12}'),
                    end_sequence=('lu_adpfwfw_rv', '${FW1R_12}${FW1_12}')),
            Adapter('LU_ADP_RVRV',
                    start_sequence=('lu_adprvrv', '${RV1R_12}${RV1_12}'),
                    end_sequence=('lu_adprvrv_rv', '${RV1R_12}${RV1_12}'))]"

  echo -e "$ADAPTER_FMT\n\n\n" > $TRIM_DIR/adapters.tmp

  sed \
  -e "/$LEAD/,/$TAIL/{ /$LEAD/{p; r $TRIM_DIR/adapters.tmp
      }; /$TAIL/p; d }"  $LONGREAD_UMI_PATH/scripts/adapters.py \
  > $TRIM_DIR/adapters.py
  
  # Add working folder to python path
  export PYTHONPATH=$PYTHONPATH:$TRIM_DIR
  
  # Perform porechop and filtlong in parallel
  FT_THREADS=$(( $THREADS/10 ))
  if (( FT_THREADS < 1 )); then
    FT_THREADS=1
  elif (( FT_THREADS > THREADS )); then
    FT_THREADS=1
  fi

  cat $READ_IN | $GNUPARALLEL --progress -j 10 -L 4 --round-robin --pipe \
    "cat > $TRIM_DIR/{#}.tmp;\
    $PORECHOP_UMI \
      -i $TRIM_DIR/{#}.tmp \
      -o $TRIM_DIR/{#}_trim.tmp \
      --threads $FT_THREADS \
      --min_split_read_size $MIN_LENGTH \
      --adapter_threshold  80 \
	  --min_trim_size 20 \
      --extra_end_trim 0 \
      --extra_middle_trim_good_side 0 \
      --extra_middle_trim_bad_side 0 \
      --middle_threshold 80 \
      --check_reads 5000; \
    $FILTLONG --min_length $MIN_LENGTH --min_mean_q 70 $TRIM_DIR/{#}_trim.tmp |\
      $CUTADAPT -j $FT_THREADS -m $MIN_LENGTH -M $MAX_LENGTH - \
        -o $TRIM_DIR/{#}_filt.tmp;"

  # Concatenate temp files
  cat $TRIM_DIR/*_filt.tmp > $TRIM_DIR/reads_tf.fq
  rm $TRIM_DIR/*.tmp
  rm -rf $TRIM_DIR/__pycache__
else
# Create symlink if already trimmed.
  ln -s $(readlink -f $READ_IN) $(readlink -f $TRIM_DIR)/reads_tf.fq  
fi

### Extract UMI references sequences ------------------------------------------- 
mkdir $OUT_DIR/umi_ref
UMI_DIR=$OUT_DIR/umi_ref

# Extract UMI terminal region
$GAWK -v UD="$UMI_DIR" 'NR%4==1{
       print $0 > UD"/reads_tf_start.fq";
       print $0 > UD"/reads_tf_end.fq";  
     }
     NR%4==2{
       print substr($0, 1, 200) > UD"/reads_tf_start.fq";
       print substr($0, length($0) - 199, 200)  > UD"/reads_tf_end.fq";  
     }
     NR%4==3{
       print $0 > UD"/reads_tf_start.fq";
       print $0 > UD"/reads_tf_end.fq";   
     }
     NR%4==0{
       print substr($0, 1, 200) > UD"/reads_tf_start.fq";
       print substr($0, length($0) - 199, 200)  > UD"/reads_tf_end.fq";  
     }
' $TRIM_DIR/reads_tf.fq


# Extract UMI pairs with correct lengths
$CUTADAPT -j $THREADS -e 0.2 -O 11 -m 18 -M 18 \
  --discard-untrimmed \
  -g $FW1...$FW2 -g $RV1...$RV2 \
  -G $RV2R...$RV1R -G $FW2R...$FW1R \
  -o $UMI_DIR/umi1.fq -p $UMI_DIR/umi2.fq \
  $UMI_DIR/reads_tf_start.fq $UMI_DIR/reads_tf_end.fq \
  > $UMI_DIR/perfect_trim.log

paste -d "" <( sed -n '1~4s/^@/>/p;2~4p' $UMI_DIR/umi1.fq ) \
            <( sed -n '1~4s/^@/>/p;2~4p' $UMI_DIR/umi2.fq ) |\
  cut -d " " -f1 > $UMI_DIR/umi12.fa

# Extract UMI pairs with correct patterns 

# Pattern: (NNNYRNNNYRNNNYRNNN NNNYRNNNYRNNNYRNNN)
PATTERN="[ATCG]{3}[CT][AG][ATCG]{3}[CT][AG][ATCG]{3}[CT][AG][ATCG]{6}\
[CT][AG][ATCG]{3}[CT][AG][ATCG]{3}[CT][AG][ATCG]{3}"
grep -B1 -E "$PATTERN" $UMI_DIR/umi12.fa |\
  sed '/^--$/d' > $UMI_DIR/umi12f.fa

# Cluster UMI pairs
$USEARCH \
  -fastx_uniques \
  $UMI_DIR/umi12f.fa \
  -fastaout $UMI_DIR/umi12u.fa \
  -sizeout \
  -minuniquesize 1 \
  -relabel umi \
  -strand both

$USEARCH \
  -cluster_fast $UMI_DIR/umi12u.fa \
  -id 0.90 \
  -centroids $UMI_DIR/umi12c.fa \
  -uc $UMI_DIR/umi12c.txt \
  -sizein \
  -sizeout \
  -strand both \
  -sort size \
  -maxaccepts 0 \
  -maxrejects 0 \
  -mincols 34 # Doesn't work...

$GAWK \
  '
  /^>/{
    SIZE=$0
    gsub(".*size=|;", "", SIZE)
    if (SIZE+0 >= 2){
      print
      getline
      print
    }
  }
  ' \
  $UMI_DIR/umi12c.fa \
  > $UMI_DIR/umi_ref.fa 
  
# Chimera screening

# Split UMIs into sub UMIs
$GAWK \
  '
    /^>/{
      HEAD=$0
      getline
      print HEAD "_1\n" substr($0,1,18) "\n" HEAD "_2\n" substr($0,19,18)
    }
  ' $UMI_DIR/umi_ref.fa \
  > $UMI_DIR/umi_ref_sub.fa

# Cluster sub UMIs to detect chimeras
$USEARCH \
  -cluster_fast $UMI_DIR/umi_ref_sub.fa \
  -id 0.94 \
  -uc $UMI_DIR/umi_ref_chimera.txt \
  -sizein \
  -sizeout \
  -strand both \
  -sort size \
  -maxaccepts 0 \
  -maxrejects 0 \
  -mincols 17

# Derivate screening

# Cluster UMIs to detect potential derivates
$USEARCH \
  -cluster_fast $UMI_DIR/umi_ref.fa \
  -id 0.83 \
  -uc $UMI_DIR/umi_ref_derivates.txt \
  -sizein \
  -sizeout \
  -strand both \
  -sort size \
  -maxaccepts 0 \
  -maxrejects 0 \
  -mincols 32 


### Bin reads based on UMIs ----------------------------------------------------
mkdir $OUT_DIR/read_binning
mkdir $OUT_DIR/read_binning/bins
BINNING_DIR=$OUT_DIR/read_binning

# Extract UMI region
$GAWK -v BD="$BINNING_DIR" -v TL="$START_READ_CHECK" '
  NR%4==1{
    print ">" substr($1,2) > BD"/reads_tf_umi1.fa";
  }
  NR%4==2{
    print substr($0, 1, TL) > BD"/reads_tf_umi1.fa";
  }
' $UMI_DIR/reads_tf_start.fq

$GAWK -v BD="$BINNING_DIR" -v TL="$END_READ_CHECK" '
  NR%4==1{
     print ">" substr($1,2) > BD"/reads_tf_umi2.fa";  
   }
   NR%4==2{
     print substr($0, length($0) - TL + 1, TL)  > BD"/reads_tf_umi2.fa";  
   }
' $UMI_DIR/reads_tf_end.fq


# Divide in barcode1 and barcode2 files
cat $UMI_DIR/umi_ref.fa <($SEQTK seq -r $UMI_DIR/umi_ref.fa |\
  $GAWK 'NR%2==1{print $0 "_rc"; getline; print};') |\
  $GAWK -v BD="$BINNING_DIR" 'NR%2==1{
       print $0 > BD"/umi_ref_b1.fa";
       print $0 > BD"/umi_ref_b2.fa";  
     }
     NR%2==0{
       print substr($0, 1, 18) > BD"/umi_ref_b1.fa";
       print substr($0, 19, 18)  > BD"/umi_ref_b2.fa";  
     }'

# Map UMIs to UMI references
## Important settings:
## -N : diasble iterative search. All possible hits are found.
## -F 20 : Removes unmapped and reverse read matches. Keeps UMIs
##         in correct orientations.

$BWA index \
  $BINNING_DIR/reads_tf_umi1.fa
$BWA \
  index \
  $BINNING_DIR/reads_tf_umi2.fa

$BWA aln \
  $BINNING_DIR/reads_tf_umi1.fa \
  $BINNING_DIR/umi_ref_b1.fa \
  -n 3 \
  -t $THREADS \
  -N \
  > $BINNING_DIR/umi1_map.sai

$BWA samse \
  -n 10000000 \
  $BINNING_DIR/reads_tf_umi1.fa \
  $BINNING_DIR/umi1_map.sai \
  $BINNING_DIR/umi_ref_b1.fa |\
$SAMTOOLS view \
  -F 20 \
  - \
  > $BINNING_DIR/umi1_map.sam

$BWA aln \
  $BINNING_DIR/reads_tf_umi2.fa \
  $BINNING_DIR/umi_ref_b2.fa \
  -n 3 \
  -t $THREADS \
  -N \
  > $BINNING_DIR/umi2_map.sai
 
$BWA samse \
  -n 10000000 \
  $BINNING_DIR/reads_tf_umi2.fa \
  $BINNING_DIR/umi2_map.sai \
  $BINNING_DIR/umi_ref_b2.fa |\
$SAMTOOLS view \
  -F 20 \
  - \
  > $BINNING_DIR/umi2_map.sam

# UMI binning and filtering

$GAWK \
  -v BD="$BINNING_DIR" \
  -v UM1="$BINNING_DIR/umi1_map.sam" \
  -v UM2="$BINNING_DIR/umi2_map.sam" \
  -v URC="$UMI_DIR/umi_ref_chimera.txt" \
  -v URD="$UMI_DIR/umi_ref_derivates.txt" \
  -v UME_MATCH_ERROR="$UMI_MATCH_ERROR" \
  -v UME_MATCH_ERROR_SD="$UMI_MATCH_ERROR_SD" \
  -v RO_FRAC="$RO_FRAC" \
  -v MAX_BIN_SIZE="$MAX_BIN_SIZE"  \
  -v BIN_CLUSTER_RATIO="$BIN_CLUSTER_RATIO" \
  '
  # Read UMI match file
  FILENAME == UM1 && FNR == 1 {
    print "[" strftime("%T") "] ### Read-UMI match filtering ###" > "/dev/stderr"
    print "[" strftime("%T") "] Reading UMI1 match file..." > "/dev/stderr"
  }
  FILENAME == UM1 {
    # Extract data from optional fields
    for (i=12; i <= NF; i++){
      # Find NM field and remove prefix (primary hit err)
      if($i ~ /^NM:i:/){sub("NM:i:", "", $i); perr = $i}
      # Find secondary hit field, remove prefix and split hits
      if($i ~ /^XA:Z:/){sub("XA:Z:", "", $i); split($i, shits, ";")}
    }
    # Add primary mapping to hit list
    err1[$1][$3]=perr
    # Add secondary mapping to hit list
    #Iterate over each hit
    for (i in shits){
      # Split hit in subdata (read, pos, cigar, err)  
      split(shits[i], tmp, ",")
      # Add hit if non-empty, not seen before and not target reverse strand
      if (tmp[1] != "" && !(tmp[1] in err1[$1]) && tmp[2] ~ "+"){
        err1[$1][tmp[1]] = tmp[4]
      }
    }
  }
  FILENAME == UM2 && FNR == 1 {
   print "[" strftime("%T") "] Reading UMI2 match file..." > "/dev/stderr"
  }
  # Read UMI match file
  FILENAME == UM2 {
    # Extract data from optional fields
    for (i=12; i <= NF; i++){
      # Find NM field and remove prefix (primary hit err)
      if($i ~ /^NM:i:/){sub("NM:i:", "", $i); perr = $i}
      # Find secondary hit field and remove prefix
      if($i ~ /^XA:Z:/){sub("XA:Z:", "", $i); split($i, shits, ";")}
    }
    # Add primary mapping to hit list
    err2[$1][$3]=perr
    # Add secondary mapping to hit list
    # Split list of hits 
    #Iterate over each hit
    for (i in shits){
      # Split hit in subdata (read, pos, cigar, err)
      split(shits[i], tmp, ",")
      # Add hit if non-empty, not seen before and not target reverse strand
      if (tmp[1] != "" && !(tmp[1] in err2[$1]) && tmp[2] ~ "+"){
        err2[$1][tmp[1]] = tmp[4]
      }
    }
  #--> Output is err1 and err2 2d arrays (umi x reads) where values are match errors
  }
  # Read chimera check file
  FILENAME == URC && $1 != "C" {
    CQUERY=$9
    sub("_.*", "", CQUERY)
    CREF=$10
    sub("_.*", "", CREF)
    if($10 == "*" && !(CQUERY in chimera_check)){
      chimera_match[CQUERY]="NA"
      chimera_check[CQUERY]="chimera_ok"
    } else if (CQUERY == CREF){
      chimera_match[CQUERY]="tandem"
      chimera_check[CQUERY]="chimera_fail"
    } else if ($10 != "*"){
      chimera_match[CQUERY]=$10
      chimera_check[CQUERY]="chimera_fail"
    }
  }
  
  # Read derivate check file
  FILENAME == URD && $1 != "C" {
    if($10 == "*"){
      derivate_match[$9]="NA"
      derivate_check[$9]="derivate_ok"
    } else {
      derivate_match[$9]=$10
      derivate_check[$9]="derivate_fail"
    }
  }
  END {
    print "[" strftime("%T") "] UMI match filtering..." > "/dev/stderr"
    # Filter reads based on UMI match error
    for (umi in err1){    
      for (read in err1[umi]){
        # Define vars
        e1 = err1[umi][read]
        e2 = err2[umi][read]
        # Filter reads not matching both UMIs
        if (e1 != "" && e2 != ""){
          # Filter based on mapping error 
          if (e1 + e2 <= 6 && e1 <= 3 && e2 <= 3){
            # Add read to bin list or replace bin assignment if error is lower
            if (!(read in match_err)){
              match_umi[read] = umi
              match_err[read] = e1 + e2
            } else if (match_err[read] > e1 + e2 ){
              match_umi[read] = umi
              match_err[read] = e1 + e2
            } 
          }
        }
      }
    }
    #--> Output is match_umi 1d array with reads [key] are linked to umi [value] 
    #--> Output is match_err 1d array with reads [key] are linked to total match err [value]
    
    # Extract read strandedness (+/-) from UMI names and count raw UMI bin assignments
    for (r in match_umi){
      UMI=match_umi[r]
      # read orientation and clean UMI name
      if (match(UMI, /_rc/) != 0){
         match_ro[r]="-"
         sub("_rc", "", UMI)
         umi_ro_neg[UMI]++
         match_umi[r] = UMI
      } else {
         match_ro[r]="+"
         umi_ro_plus[UMI]++
      }
      # Count reads pr UMI
      umi_n_raw[UMI]++
    }
    
    # Read orientation filtering 
    if (RO_FRAC != 0){
      print "[" strftime("%T") "] Read orientation filtering..." > "/dev/stderr"
    
      # Calculate read orientation fraction
      for (u in umi_ro_plus){
        # Check read orientation fraction
        if (umi_ro_plus[u] >= 1 && umi_ro_neg[u] >= 1){
          if (umi_ro_plus[u]/(umi_ro_neg[u]+umi_ro_plus[u]) < RO_FRAC ){
            rof_check[u]="rof_subset"
            ROF_N = umi_ro_plus[u]*(1/RO_FRAC-1)
            rof_sub_target[u] = ROF_N
            rof_sub_neg_n[u] = ROF_N
            rof_sub_pos_n[u] = ROF_N
          } else if (umi_ro_neg[u]/(umi_ro_neg[u]+umi_ro_plus[u]) < RO_FRAC ){
            rof_check[u]="rof_subset"
            ROF_N = umi_ro_neg[u]*(1/RO_FRAC-1)
            rof_sub_target[u] = ROF_N
            rof_sub_neg_n[u] = ROF_N
            rof_sub_pos_n[u]= ROF_N
          } else {
            rof_check[u]="rof_ok"
            rof_sub_target[u] = "NA"
          }
        } else {
          rof_check[u]="rof_fail"
          rof_sub_target[u] = "NA"
        }
      }
      
      # Subset reads
      for (r in match_umi){
        UMI=match_umi[r]
        if (rof_sub_target[UMI] != "NA"){
          if(match_ro[r] == "+"){
            if(rof_sub_pos_n[UMI]-- <= 0){
              # Remove unused reads from match_umi/match_err arrays
              delete match_umi[r]
              delete match_err[r]
            }
          } else if (match_ro[r] == "-"){
            if(rof_sub_neg_n[UMI]-- <= 0){
              # Remove unused reads from match_umi/match_err arrays
              delete match_umi[r]
              delete match_err[r]
            }
          }
        }
      }
    } else {
      for (u in umi_n_raw){
        rof_check[u]="rof_disabled"
        rof_sub_target[u]="NA"
      }
    }

    print "[" strftime("%T") "] UMI match error filtering..." > "/dev/stderr"

    # Calculate UME stats
    for (r in match_umi){
      UMI=match_umi[r]
      # UMI match error stats
      umi_me_sum[UMI] += match_err[r]
      umi_me_sq[UMI] += (match_err[r])^2
      # Create list of UMIs
      umi_n[UMI]++ 
    }

    # Check UMI match error
    for (u in umi_n){
      ume_mean[u] = umi_me_sum[u]/umi_n[u]
      ume_sd[u] = sqrt((umi_me_sq[u]-umi_me_sum[u]^2/umi_n[u])/umi_n[u])
      if (ume_mean[u] > UME_MATCH_ERROR || ume_sd[u] > UME_MATCH_ERROR_SD){
        ume_check[u] = "ume_fail"
      } else {
        ume_check[u] = "ume_ok"
      }
    }

    print "[" strftime("%T") "] UMI bin/cluster size ratio filtering..." > "/dev/stderr"
    for (u in umi_n){
      CLUSTER_SIZE=u
      gsub(".*;size=|;", "", CLUSTER_SIZE)
      bcr[u]=umi_n_raw[u]/CLUSTER_SIZE
      if (bcr[u] > BIN_CLUSTER_RATIO){
        bcr_check[u] = "bcr_fail"
      } else {
        bcr_check[u] = "bcr_ok"
      }
    }

    # Print filtering stats
    print \
      "umi_name",\
      "read_n",\
      "read_plus_n",\
      "read_neg_n",\
      "read_or_ratio",\
      "read_or_filter",\
      "read_or_sub_n",\
      "umi_match_err_mean",\
      "umi_match_err_sd",\
      "umi_match_err_filter",\
      "bin_cluster_ratio",\
      "bin_cluster_ratio_filter",\
      "chimera_match",\
      "chimera_filter",\
      "derivate_match",\
      "derivate_filter",\
      "filter_result" \
      > BD"/umi_binning_stats.txt"
    for (u in umi_n_raw){
      if (rof_check[u] ume_check[u] bcr_check[u] chimera_check[u] derivate_check[u] ~ /fail/){
        filter_result[u] = "fail"
      } else {
        filter_result[u] = "ok"
      }
      print \
        u,\
        umi_n_raw[u]+0,\
        umi_ro_plus[u]+0,\
        umi_ro_neg[u]+0,\
        umi_ro_plus[u]/(umi_ro_neg[u]+umi_ro_plus[u])+0, \
        rof_check[u],\
        rof_sub_target[u],\
        ume_mean[u],\
        ume_sd[u],\
        ume_check[u],\
        bcr[u],\
        bcr_check[u],\
        chimera_match[u],\
        chimera_check[u],\
        derivate_match[u],\
        derivate_check[u],\
        filter_result[u] \
        > BD"/umi_binning_stats.txt"
    }

    print "[" strftime("%T") "] Print UMI matches..." > "/dev/stderr" 
    for (r in match_umi){
      UMI=match_umi[r]
      print UMI, r, match_err[r]
    }
    
    # Print to terminal
    print "[" strftime("%T") "] Done." > "/dev/stderr" 
  }
' \
$BINNING_DIR/umi1_map.sam \
$BINNING_DIR/umi2_map.sam \
$UMI_DIR/umi_ref_derivates.txt \
$UMI_DIR/umi_ref_chimera.txt \
> $BINNING_DIR/umi_bin_map.txt

# Extract binned reads

umi_binning() {
  # Input
  local READS=$1
  local UMIMAP=$2
  local OUT=$3

  # Binning
  $GAWK '
    NR==FNR{
      UMI_SUBSET[$0]=""
      next
    }
    {
      if ($1 in UMI_SUBSET){
        print $0
      }
    }
  ' - $UMIMAP |\
  $GAWK -v out="$OUT" '
    NR==FNR {
      # Get read name
      sub(";.*", "", $1);
      # Associate read name and umi match
      bin[$2]=$1;
      next;
    }
    FNR%4==1 {
      read=substr($1,2);
      bin_tmp=bin[read]
      if ( bin_tmp != "" ){
        binfile=out"/"bin_tmp"bins.fastq";
        READ_RECORD=$0
        getline; READ_RECORD=READ_RECORD"\n"$0
        getline; READ_RECORD=READ_RECORD"\n"$0
        getline; print READ_RECORD"\n"$0 > binfile;
      }
    }
  ' - $READS
}

export -f umi_binning

cut -d " " -f1 $BINNING_DIR/umi_bin_map.txt |\
  sort -u |\
  $GNUPARALLEL \
    --env umi_binning \
    -N 4000 \
	-j $THREADS \
	--pipe \
  "
    mkdir $BINNING_DIR/bins/{#}
    cat |\
      umi_binning \
        $TRIM_DIR/reads_tf.fq \
        $BINNING_DIR/umi_bin_map.txt \
        $BINNING_DIR/bins/{#}
  "
