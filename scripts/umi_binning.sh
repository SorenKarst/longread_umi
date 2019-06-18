#!/bin/bash

# DESCRIPTION
#    Script for binning Nanopore reads based on UMIs. Part of the 
#    longread-UMI-pipeline.
#    
# IMPLEMENTATION
#    author	Søren Karst (sorenkarst@gmail.com)
#               Ryans Ziels (ziels@mail.ubc.ca)
#    license	GNU General Public License
#
# TO DO
#    Add terminal messages.
#    Optimize trimming and filtering for speed.
#    Add bin size limit to 200 x
#    Add mapping against adaptors to remove UMI artifacts

### Terminal input ------------------------------------------------------------
READ_IN=$1
OUT_DIR=$2
THREADS=$3
MIN_LENGTH=$4
MAX_LENGTH=$5
TERMINAL1_CHECK_RANGE=${6:-70}
TERMINAL2_CHECK_RANGE=${7:-80}
FW1=${8:-CAAGCAGAAGACGGCATACGAGAT} #RC: ATCTCGTATGCCGTCTTCTGCTTG
FW2=${9:-AGRGTTYGATYMTGGCTCAG} #RC: CTGAGCCAKRATCRAACYCT
RV1=${10:-AATGATACGGCGACCACCGAGATC} #RC: GATCTCGGTGGTCGCCGTATCATT
RV2=${11:-CGACATCGAGGTGCCAAAC} #RC: GTTTGGCACCTCGATGTCG


### Primer formating
revcom() {
  echo $1 |\
  awk '{print ">dummy\n" $0}' |\
  $SEQTK seq -r - |\
  awk '!/^>/'  
}
FW1R=$(revcom "$FW1")
FW2R=$(revcom "$FW2")
RV1R=$(revcom "$RV1")
RV2R=$(revcom "$RV2")

### Read trimming and filtering -----------------------------------------------
mkdir $OUT_DIR
TRIM_DIR=$OUT_DIR/trim

# Check if data is trimmed
if [ ! -f "$TRIM_DIR/reads_tf.fq" ]; then
  # Create trim dir
  mkdir $TRIM_DIR

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
      --extra_end_trim 0 \
      --extra_middle_trim_good_side 0 \
      --extra_middle_trim_bad_side 0 \
      --middle_threshold 80 \
      --check_reads 1000; \
    $FILTLONG --min_length $MIN_LENGTH --min_mean_q 70 $TRIM_DIR/{#}_trim.tmp |\
      $CUTADAPT -j $FT_THREADS -m $MIN_LENGTH -M $MAX_LENGTH - \
        -o $TRIM_DIR/{#}_filt.tmp;"

  # Concatenate temp files
  cat $TRIM_DIR/*_filt.tmp > $TRIM_DIR/reads_tf.fq
  rm $TRIM_DIR/*.tmp
fi

### Extract UMI references sequences ------------------------------------------- 
mkdir $OUT_DIR/umi_ref
UMI_DIR=$OUT_DIR/umi_ref

# Extract UMIs and 5' flanking region
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

# Extract UMIs with correct length in both ends
$CUTADAPT -j $THREADS -e 0.2 -O 11 -m 18 -M 18 \
  --discard-untrimmed \
  -g $FW1...$FW2 -g $RV1...$RV2 \
  -G $RV2R...$RV1R -G $FW2R...$FW1R \
  -o $UMI_DIR/umi1.fq -p $UMI_DIR/umi2.fq \
  $UMI_DIR/reads_tf_start.fq $UMI_DIR/reads_tf_end.fq \
  > $UMI_DIR/trim.log

paste -d "" <( sed -n '1~4s/^@/>/p;2~4p' $UMI_DIR/umi1.fq ) \
            <( sed -n '1~4s/^@/>/p;2~4p' $UMI_DIR/umi2.fq ) |\
  cut -d " " -f1 > $UMI_DIR/umi12.fa

# Extract UMIs with correct patterns 

# Pattern: (NNNYRNNNYRNNNYRNNN NNNYRNNNYRNNNYRNNN)
PATTERN="[ATCG]{3}[CT][AG][ATCG]{3}[CT][AG][ATCG]{3}[CT][AG][ATCG]{6}\
[CT][AG][ATCG]{3}[CT][AG][ATCG]{3}[CT][AG][ATCG]{3}"
grep -B1 -E "$PATTERN" $UMI_DIR/umi12.fa |\
  sed '/^--$/d' > $UMI_DIR/umi12f.fa

# Extract UMIs with abundance >= 2 
$USEARCH -fastx_uniques $UMI_DIR/umi12f.fa \
  -fastaout $UMI_DIR/umi12u.fa -sizeout -minuniquesize 2 \
  -relabel umi -strand both

$USEARCH -cluster_fast $UMI_DIR/umi12u.fa -id 0.85 \
  -centroids $UMI_DIR/umi12c.fa -uc $UMI_DIR/umi12c.txt \
  -sizein -sizeout -strand both 


# Remove potential chimeras
paste <(cat $UMI_DIR/umi12c.fa | paste - - ) \
  <($GAWK '!/^>/{print}' $UMI_DIR/umi12c.fa | rev | tr ATCG TAGC) |\
  $GAWK -v UD="$UMI_DIR" 'NR==FNR {
      #Format columns
      split($1, a, /[>;]/);
      sub("size=", "", a[3]);
      # Extract UMI1 and UMI2 in both orientations
      s1 = substr($2, 1, 18);
      s2 = substr($2, 19, 36);
      s1rc= substr($3, 1, 18);
      s2rc= substr($3, 19, 36);
      # Register UMI1 size if larger than current highest or if empty
      if ((g1n[s1]+0) <= (a[3]+0) || g1n[s1] == ""){
        g1n[s1] = a[3];
        g1[s1] = a[2];
      }
      # Register UMI2 size if larger than current highest or if empty
      if ((g2n[s2]+0) <= (a[3]+0) || g2n[s2] == ""){
        g2n[s2] = a[3];
        g2[s2] = a[2];
      }
      # Register UMI1rc size if larger than current highest or if empty
      if ((g1n[s1rc]+0) <= (a[3]+0) || g1n[s1rc] == ""){
        g1n[s1rc] = a[3];
        g1[s1rc] = a[2];
      }
      # Register UMI2rc size if larger than current highest or if empty
      if ((g2n[s2rc]+0) <= (a[3]+0) || g2n[s2rc] == ""){
        g2n[s2rc] = a[3];
        g2[s2rc] = a[2];
      }
      # Register UMI1 and UMI matches for current UMI
      u[a[2]] = a[3];
      s1a[a[2]] = s1;
      s2a[a[2]] = s2;
      s1arc[a[2]] = s1rc;
      s2arc[a[2]] = s2rc;
    } END {
      for (i in u){
        keep="no";
        if (g1[s1a[i]] == i && g2[s2a[i]] == i && g1[s1arc[i]] == i && g2[s2arc[i]] == i && s1a[i] != s1arc[i]){
          keep="yes";
          print ">"i";"u[i]"\n"s1a[i]s2a[i] > UD"/umi_ref.fa";
        } else if (s1a[i] == s1arc[i]){
          keep="tandem"
          print ">"i";"u[i]"\n"s1a[i]s2a[i] > UD"/umi_ref.fa";
        }
        print i, n[i], s1a[i], s2a[i], keep, g1[s1a[i]]"/"g2[s2a[i]]"/"g1[s1arc[i]]"/"g2[s2arc[i]], u[i]
      }  
    }' > $UMI_DIR/umi_ref.txt

# Checks: grep -c "no" umi_ref/umi_ref.txt; grep -c "yes" umi_ref/umi_ref.txt;
# $GAWK -F";" 'FNR==NR{T+=$2}END{print T}' umi_ref/umi_ref.fa

### Bin reads based on UMIs ----------------------------------------------------
mkdir $OUT_DIR/read_binning
mkdir $OUT_DIR/read_binning/bins
BINNING_DIR=$OUT_DIR/read_binning

# Extract UMI region
$GAWK -v BD="$BINNING_DIR" -v TL="$TERMINAL1_CHECK_RANGE" '
  NR%4==1{
    print ">" substr($1,2) > BD"/reads_tf_umi1.fa";
  }
  NR%4==2{
    print substr($0, 1, TL) > BD"/reads_tf_umi1.fa";
  }
' $UMI_DIR/reads_tf_start.fq

$GAWK -v BD="$BINNING_DIR" -v TL="$TERMINAL2_CHECK_RANGE" '
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

$BWA index $BINNING_DIR/reads_tf_umi1.fa
$BWA index $BINNING_DIR/reads_tf_umi2.fa

$BWA aln $BINNING_DIR/reads_tf_umi1.fa $BINNING_DIR/umi_ref_b1.fa \
  -n 3 -t $THREADS -N > $BINNING_DIR/umi1_map.sai
$BWA samse -n 10000000 $BINNING_DIR/reads_tf_umi1.fa $BINNING_DIR/umi1_map.sai\
  $BINNING_DIR/umi_ref_b1.fa | $SAMTOOLS view -F 20 - > $BINNING_DIR/umi1_map.sam

$BWA aln $BINNING_DIR/reads_tf_umi2.fa $BINNING_DIR/umi_ref_b2.fa \
  -n 3 -t $THREADS -N > $BINNING_DIR/umi2_map.sai
$BWA samse -n 10000000 $BINNING_DIR/reads_tf_umi2.fa $BINNING_DIR/umi2_map.sai\
  $BINNING_DIR/umi_ref_b2.fa | $SAMTOOLS view -F 20 - > $BINNING_DIR/umi2_map.sam

# Merge umi1 and 2 mappings

$GAWK -v BD="$BINNING_DIR" '
  # Print to terminal
  NR==1 {
    print "[" strftime("%T") "] ### Read-UMI match filtering ###" > "/dev/stderr";
    print "[" strftime("%T") "] Reading UMI1 match file..." > "/dev/stderr";
  }
  # Read UMI match file
  NR==FNR{
    # Extract data from optional fields
    for (i=12; i <= NF; i++){
      # Find NM field and remove prefix (primary hit err)
      if($i ~ /^NM:i:/){sub("NM:i:", "", $i); perr = $i};
      # Find secondary hit field, remove prefix and split hits
      if($i ~ /^XA:Z:/){sub("XA:Z:", "", $i); split($i, shits, ";")};
    }
    # Add primary mapping to hit list
    err1[$1][$3]=perr;
    # Add secondary mapping to hit list
    #Iterate over each hit
    for (i in shits){
      # Split hit in subdata (read, pos, cigar, err)  
      split(shits[i], tmp, ",");
      # Add hit if non-empty, not seen before and not target reverse strand
      if (tmp[1] != "" && !(tmp[1] in err1[$1]) && tmp[2] ~ "+"){
        err1[$1][tmp[1]] = tmp[4];
      }
    }
    next;
  }
  # Print to terminal (NB: Require next in above statement)
  FNR==1 {
   print "[" strftime("%T") "] Reading UMI2 match file..." > "/dev/stderr";
  }
  # Read UMI match file
  {
    # Extract data from optional fields
    for (i=12; i <= NF; i++){
      # Find NM field and remove prefix (primary hit err)
      if($i ~ /^NM:i:/){sub("NM:i:", "", $i); perr = $i};
      # Find secondary hit field and remove prefix
      if($i ~ /^XA:Z:/){sub("XA:Z:", "", $i); split($i, shits, ";")};
    }
    # Add primary mapping to hit list
    err2[$1][$3]=perr;
    # Add secondary mapping to hit list
    # Split list of hits 
    #Iterate over each hit
    for (i in shits){
      # Split hit in subdata (read, pos, cigar, err)
      split(shits[i], tmp, ",");
      # Add hit if non-empty, not seen before and not target reverse strand
      if (tmp[1] != "" && !(tmp[1] in err2[$1]) && tmp[2] ~ "+"){
        err2[$1][tmp[1]] = tmp[4];
      }
    }
  } END {
    # Print to terminal
    print "[" strftime("%T") "] UMI match filtering..." > "/dev/stderr"; 
    # Loop over UMIs and reads
    for (umi in err1){    
      for (read in err1[umi]){
        # Define vars
        e1 = err1[umi][read];
        e2 = err2[umi][read];
        # Filter reads not matching both UMIs
        if (e1 != "" && e2 != ""){
          # Filter based on mapping error 
          if (e1 + e2 <= 6 && e1 <= 3 && e2 <= 3){
            # Add read to bin list or replace bin assignment if error is lower
            if (!(read in match_err)){
              match_umi[read] = umi;
              match_err[read] = e1 + e2;
            } else if (match_err[read] > e1 + e2 ){
              match_umi[read] = umi;
              match_err[read] = e1 + e2;
            } 
          }
        }
      }
    }
    # Print to terminal 
    print "[" strftime("%T") "] Print UMI matches..." > "/dev/stderr"; 
    # Print binned reads
    for (s in match_umi){
      sub("_rc", "", match_umi[s]);
      umi_n[match_umi[s]]++;
      print match_umi[s], s, match_err[s]
    }
    # Print stats
    for (u in umi_n) print u, umi_n[u] > BD"/umi_binning_stats.txt";    
    # Print to terminal
    print "[" strftime("%T") "] Done." > "/dev/stderr"; 
  }
' $BINNING_DIR/umi1_map.sam $BINNING_DIR/umi2_map.sam > $BINNING_DIR/umi_bin_map.txt


# Extract binned reads

umi_binning() {
  # Input
  local UMIMAP=$1
  local OUT=$2

  # Binning
  $GAWK -v out="$OUT" '
    BEGIN {g=1; outsub="./"out"/"g; system("mkdir \047" outsub "\047");}
    NR==FNR {
      # Get read name
      sub(";.*", "", $1);
      # Associate read name and umi match
      bin[$2]=$1;
      # Assign umi to a folder group if it has none
      if (foldergrp[$1] == ""){
        j++;
        if (j <= 4000){
          foldergrp[$1]=g;
        } else {
          j = 0;
          g++;
          foldergrp[$1]=g;
          outsub="./"out"/"g;
          system("mkdir \047" outsub "\047");
        }
      }
      next;
    }
    FNR%4==1 {
      read=substr($1,2);
      bin_tmp=bin[read]
      if ( bin_tmp != "" ){
        binfile=out"/"foldergrp[bin_tmp]"/"bin_tmp"bins.fastq";
        print > binfile;
        getline; print > binfile;
        getline; print > binfile;
        getline; print > binfile;
      }
    }
  ' $UMIMAP -
}

export -f umi_binning

cat $TRIM_DIR/reads_tf.fq | $GNUPARALLEL -L4 -j $THREADS --block 300M --pipe \
  "mkdir $BINNING_DIR/bins/job{#}; cat | umi_binning $BINNING_DIR/umi_bin_map.txt\
  $BINNING_DIR/bins/job{#}"

aggregate_bins() {
  # Input
  local IN=$1
  local OUTDIR=$2
  local OUTNAME=$3
  local JOB=$4

  # Determine output folder
  local BIN=$(( ($JOB - 1)/4000 ))
  mkdir -p $OUTDIR/$BIN

  # Aggregate data
  cat $IN > $OUTDIR/$BIN/$OUTNAME  
}

export -f aggregate_bins

find $BINNING_DIR/bins/*/*/ -name "*bins.fastq" -printf "%f\n" |\
  sort | uniq | $GNUPARALLEL -j $THREADS "aggregate_bins '$BINNING_DIR/bins/*/*/'{/} \
  $BINNING_DIR/bins {/} {#}"

rm -r $BINNING_DIR/bins/job*


