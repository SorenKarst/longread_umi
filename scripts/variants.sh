#!/bin/bash
# DESCRIPTION
#    Phasing and variant calling of consensus sequences to obtain
#    lowest possible error rate
#    
# IMPLEMENTATION
#    author	Søren Karst (sorenkarst@gmail.com)
#               Ryans Ziels (ziels@mail.ubc.ca)
#    version	0.1.0
#    license	GNU General Public License
# TODO
# - extract_vars and bin_vars can greately simplified by proper use
#   of samtools mpileup which references the positions in the reads
#   already.

### Terminal input ------------------------------------------------------------
CONSENSUS_FILE=$1 #Consensus file path name
OUT_DIR=$2 #output folder name
THREADS=$3 #Number of threads

### Custom functions ----------------------------------------------------------
bam_read_split() {
  ### Description
  # Input is a bam file containing reads mapped to a reference.
  # The function splits the bam file per read. The output is per
  # read mapping information in sam format.
 
  # Input 
  IN_BAM=$1

  # Split
  $SAMTOOLS view -h $IN_BAM |\
    $GAWK 'NR<=3{h[++l]=$0; next}
      NR>3{s[++n]=$0; next}
      END {
        for (i in s){
          print h[1]"\n"h[2]"\n"h[3]"\n"s[i]
      }
    }'
}
export -f bam_read_split

extract_vars(){
  ### Description
  # Input is a sam file of a single read mapped to a reference.
  # The function extracts base information at diagnostic
  # coordinates based on a reference positions. The output is a 
  # read name and diagnostic bases, which are used as a bin id.
   
  # Input
  local IN=$(cat)
  local REG=$1
  local REF=$2

  # Format name
  local UMI_NAME=$(echo "$IN" |\
    $GAWK 'NR==4{print $1}')

  # mpileup
  if [ -z "$REG" ]; then
    echo "${UMI_NAME} all"
  elif [ ! -z "$REG" ]; then
    echo "$IN" |\
      $SAMTOOLS view -b - |\
      $BCFTOOLS mpileup \
        -Ov \
        -a "FORMAT/AD,FORMAT/DP" \
        -f $REF \
        -t $REG \
        - |\
      $GAWK -v un="$UMI_NAME" '!/^#/{
        gsub(",|<\\*>", "", $5)
        df=length($4)-length($5)
        if ($5 == ""){p = p $2 $4}
        else if (df > 0){
          for(c=1;c<=df;c++){d=d"-"}
          p = p $2 $5 d
        }
        else if (df <= 0){p = p $2 $5}
      } END {
        print un, p
      }'
  fi
  # Testing
  # samtools view -h Cluster25.bam | head -n 4 | extract_vars "$REG" Cluster25_con.fa > test.txt
}
export -f extract_vars

var_bin() {
  ### Description
  # Input is file of reads in fasta format and a list of read names
  # and diagnostic bases (bin id). The extracts reads and divides them
  # into thier assigned bins. 

  # Input
  local IN=$1
  local READS=$2
  local CN=$3

  # Binning
  cat $IN |\
  $GAWK -v out="$CN" 'NR==FNR{
    seq[$1] = $2
    cl[$2]++
    next
  }
  FNR%2==1{
    # Read consensus read name
    sn = substr($0, 2);
    # Check read has been binned
    if (sn in seq){
      # Write read to varians bin if bin is big enough
      if (cl[seq[sn]] >= 3){
        print ">" sn > out"_"seq[sn]"_bin.fa"
        getline
        print > out"_"seq[sn]"_bin.fa"
      }
      # Write read to none bin if bin is too small
      else if (cl[seq[i]] < 3){
        print ">" sn > out"_none.fa"
        getline
        print > out"_none.fa"
      }
    sn = ""
    }    
  }' - <($SEQTK seq -l0 $READS)
  
  # Testing
  #samtools view -h Cluster0.bam | head -n 4 |\
  # extract_vars "$REG" Cluster0_con.fa | var_bin \
  #../../../final/consensus_sracon_medaka_medaka_b30.fa test
  # seqtk seq -l0 ../../clusters/Cluster0 > test.fa
}
export -f var_bin

phased_consensus(){
  ### Description
  # Input is bin file of reads in fasta format and an output
  # name. The binned reads are used to create a raw consensus
  # using usearch, which is polished using bcftools.

  # Input
  local PB=$1
  local PB_NAME=$2

  # Format
  local OUT=${PB%/*}

  # Usearch raw consensus
  $USEARCH \
    -cluster_fast $PB \
    -id 0.99 \
    -strand both \
    -consout $OUT/${PB_NAME}_bincon.fa \
    -sizeout \
    -relabel ${PB_NAME}

  # Polish with $BCFTOOLS
  $MINIMAP2 \
    -t 1 \
    -ax map-ont \
    $OUT/${PB_NAME}_bincon.fa \
    $PB |\
    $SAMTOOLS view -b - \
    > $OUT/${PB_NAME}_temp.bam

  $SAMTOOLS sort \
    -o $OUT/${PB_NAME}.bam \
    $OUT/${PB_NAME}_temp.bam
  $SAMTOOLS index $OUT/${PB_NAME}.bam
  rm $OUT/${PB_NAME}_temp.bam   

  $BCFTOOLS mpileup \
    -Ov \
    -d 1000000 \
    -L 1000000 \
    -a "FORMAT/AD,FORMAT/DP" \
    -f $OUT/${PB_NAME}_bincon.fa \
    $OUT/${PB_NAME}.bam \
    > $OUT/${PB_NAME}_pileup.vcf

  $BCFTOOLS norm \
    -f $OUT/${PB_NAME}_bincon.fa \
    $OUT/${PB_NAME}_pileup.vcf \
    -Ov \
    -o $OUT/${PB_NAME}_norm.vcf

  $BCFTOOLS view -i 'AD[0:1]/FORMAT/DP>0.5' \
    -Oz \
    $OUT/${PB_NAME}_norm.vcf \
    > $OUT/${PB_NAME}_polish.vcf.gz
  $BCFTOOLS index $OUT/${PB_NAME}_polish.vcf.gz
    
  $BCFTOOLS consensus $OUT/${PB_NAME}_polish.vcf.gz \
    -f $OUT/${PB_NAME}_bincon.fa |\
    $SEQTK seq -l0 - |\
      $GAWK -v name="$PB_NAME" -v out="$OUT" '
        /^>/{
          ++n;
          if (n == 1){
            split($0,size,";")
            gsub("size=", ";", size[2])
            print ">" name size[2] > out "/" name "_phased.fa"
            getline
            print > out "/" name "_phased.fa"
          }
        } END {
          if (n > 1){msg = "warning"} else {msg="ok"}
          print name, n, msg >> out "/bincon_log.txt"
        }'
}
export -f phased_consensus

phasing () {
  # Input
  local IN=$(cat) # Single cluster sequence
  local VARIANT_OUT=$1 # Output dir for cluster
  local CLUSTER_READ_DIR=$2 # Path to cluster reads
  local CONSENSUS_FILE=$3 # Path to consensus reads
  local CLUSTER_THREADS=$4 # Number of threads allocated

  # Name format
  local CLUSTER_NAME=$(echo "$IN" |\
    $GAWK 'NR==1{gsub(">|;.*","",$0); print $0}')

  # Prepare output folder
  local CLUSTER_OUT=$VARIANT_OUT/${CLUSTER_NAME}
  mkdir $CLUSTER_OUT
  echo "$IN" > $CLUSTER_OUT/${CLUSTER_NAME}_con.fa

  # Map cluster reads to cluster sequences
  $MINIMAP2 \
    -t $CLUSTER_THREADS \
    -ax asm5\
    $CLUSTER_OUT/${CLUSTER_NAME}_con.fa \
    $CLUSTER_READ_DIR/${CLUSTER_NAME} |\
    $SAMTOOLS view -b - \
    > $CLUSTER_OUT/temp.bam

  $SAMTOOLS sort \
    -o $CLUSTER_OUT/${CLUSTER_NAME}.bam \
    $CLUSTER_OUT/temp.bam
  $SAMTOOLS index $CLUSTER_OUT/${CLUSTER_NAME}.bam
  rm $CLUSTER_OUT/temp.bam   

  # Call variants
  $BCFTOOLS mpileup \
    -Ov \
    -d 1000000 \
    -L 1000000 \
    -a "FORMAT/AD,FORMAT/DP" \
    -f $CLUSTER_OUT/${CLUSTER_NAME}_con.fa \
    $CLUSTER_OUT/${CLUSTER_NAME}.bam \
    > $CLUSTER_OUT/variants.vcf
  $BCFTOOLS view -i 'AD[0:1-]>2' \
    $CLUSTER_OUT/variants.vcf \
    > $CLUSTER_OUT/calls.vcf

  # Convert variants to list
  REG=$(
    $GAWK '
      # Skip comment lines
      !/^#/ {
        reg[$1":"$2]++
       } END {
         for (i in reg){
           if(++n == 1) regout=i
           else if(++n > 1) regout=regout","i
         }
         print regout
       }' $CLUSTER_OUT/calls.vcf

      # Homopolymers effects should have been removed.
      # If not test: && /tolower($5) !~ /a{3,}|t{3,}|c{3,}|g{3,}/
      # but remember this code does not take SNPs close to hps in to account
      # Exampel: GAAA -> AAAA
  )

  # Phasing

  ## Export extract_vars and dependencies
  export -f extract_vars
  export SAMTOOLS=$SAMTOOLS
  export BCFTOOLS=$BCFTOOLS
  export GAWK=$GAWK
  export SEQTK=$SEQTK 

  bam_read_split $CLUSTER_OUT/${CLUSTER_NAME}.bam |\
    $GNUPARALLEL -j $CLUSTER_THREADS -L 4 -N 1 --pipe \
    "cat | extract_vars \"$REG\" \"$CLUSTER_OUT/${CLUSTER_NAME}_con.fa\"" \
    > $CLUSTER_OUT/read_variants.txt

  var_bin \
    $CLUSTER_OUT/read_variants.txt \
    $CONSENSUS_FILE \
    $CLUSTER_OUT/${CLUSTER_NAME}

  # Consensus from phased reads

  ## Export phased_consensus and dependencies
  export -f phased_consensus
  export SAMTOOLS=$SAMTOOLS
  export BCFTOOLS=$BCFTOOLS
  export GAWK=$GAWK
  export MINIMAP2=$MINIMAP2
  export USEARCH=$USEARCH

  find $CLUSTER_OUT -type f -name "*_bin.fa" |\
    $GNUPARALLEL -j $CLUSTER_THREADS \
    --rpl '{name} s:.*/::; s/_bin.fa$//' \
    "phased_consensus {} {name}"
}
export -f phasing

subsample_phasing(){
  # Input
  local IN_DIR="$1"
  local OUT_DIR=$2
  local SAMPLE_N=$3
  local THREADS=$4

  # Subsample to n
  $GAWK -v out="$OUT_DIR" -v sn="$SAMPLE_N" '
    /^>/{ 
      h=$0
      getline
      s=$0
      a[++n]=h"\n"s
    }
    ENDFILE{
      srand();
      fn=FILENAME
      gsub(".*/|_bin.fa", "", fn)
      while ( length(a) >= sn){
        ++bf
        for (i = 1; i <= sn; i++){
          sp=int(1+rand()*length(a))
          print a[sp] > out "/" fn "_" bf "_sample.fa"
          delete a[sp]
          asort(a)
        }
      }
      n=0
      bf=0  
    }
  ' $IN_DIR/*_bin.fa 

  ## Export phased_consensus and dependencies
  export -f phased_consensus
  export SAMTOOLS=$SAMTOOLS
  export BCFTOOLS=$BCFTOOLS
  export GAWK=$GAWK
  export MINIMAP2=$MINIMAP2
  export USEARCH=$USEARCH

  find $OUT_DIR -type f -name "*_sample.fa" |\
    $GNUPARALLEL -j $THREADS \
      --rpl '{name} s:.*/::; s/_sample.fa$//' \
      "phased_consensus {} {name}"
}
export -f subsample_phasing

### Consensus phasing pipeline -------------------------------------------------

# Threads handling
CLUSTER_JOBS=10
CLUSTER_THREADS=$(( $THREADS/$CLUSTER_JOBS ))

# Format names
CONSENSUS_NAME=${CONSENSUS_FILE##*/}
CONSENSUS_NAME=${CONSENSUS_NAME%.*}

# Prepare output folders
if [ -d "$OUT_DIR" ]; then
  echo "Output folder exists. Exiting..."
  exit 0
fi
mkdir $OUT_DIR
mkdir $OUT_DIR/clusters

# Hard mask homopolymers
$SEQTK seq -l0 $CONSENSUS_FILE |\
$GAWK '
  /^>/{print}
  !/^>/{
    gsub(/A{3,}/, "AA", $0)
    gsub(/T{3,}/, "TT", $0)
    gsub(/C{3,}/, "CC", $0)
    gsub(/G{3,}/, "GG", $0)
    print $0
  }' > $OUT_DIR/${CONSENSUS_NAME}_m_temp.fa

# Usearch cluster round 1
$USEARCH \
  -cluster_fast $OUT_DIR/${CONSENSUS_NAME}_m_temp.fa \
  -id 0.995 \
  -strand both \
  -consout $OUT_DIR/${CONSENSUS_NAME}_c1_temp.fa \
  -clusters $OUT_DIR/clusters/Cluster \
  -sort length \
  -sizeout

# Usearch cluster round 2
$USEARCH \
  -cluster_fast $OUT_DIR/${CONSENSUS_NAME}_c1_temp.fa \
  -id 0.995 \
  -strand both \
  -centroids $OUT_DIR/${CONSENSUS_NAME}_c2_temp.fa\
  -uc $OUT_DIR/merge_c2_temp.uc \
  -sort length \
  -sizein \
  -sizeout

# Merge clusters
$GAWK -v cd="$OUT_DIR/clusters/" '
  $10 !~ /\*/{
    sub(";.*", "", $9)
    sub(";.*", "", $10)
    system("cat " cd $9 " >> " cd $10)
    system("rm " cd $9)
  }' $OUT_DIR/merge_c2_temp.uc 

# Filter clusters
$USEARCH \
  -sortbysize \
  $OUT_DIR/${CONSENSUS_NAME}_c2_temp.fa \
  -fastaout $OUT_DIR/${CONSENSUS_NAME}_clusters.fa \
  -minsize 3
rm $OUT_DIR/*temp*

# Phasing and consensus
VARIANT_OUT=$OUT_DIR/phasing_consensus
mkdir -p $VARIANT_OUT

cat $OUT_DIR/${CONSENSUS_NAME}_clusters.fa | $SEQTK seq -l0 - |\
  $GNUPARALLEL --progress -j$CLUSTER_JOBS --recstart ">" -N 1 --pipe \
  "cat | phasing $VARIANT_OUT $OUT_DIR/clusters \
   $CONSENSUS_FILE $CLUSTER_THREADS"
cat $VARIANT_OUT/*/*variants.fa > $OUT_DIR/variants_all.fa

### Testing
exit 0

### Sample phasing pipeline -------------------------------------------------
NR=(3 5 7)
for i in "${NR[@]}"; do
  SUBSAMPLES_OUT=$OUT_DIR/subsampling_n$i
  mkdir $SUBSAMPLES_OUT
  subsample_phasing \
    $VARIANT_OUT'/*/' \
    $SUBSAMPLES_OUT \
    $i \
    $THREADS
  cat $SUBSAMPLES_OUT/*variants.fa > $OUT_DIR/variants_n${i}.fa
done

# Soft mask homopolymers 
#seqtk seq -l0 $CONSENSUS_FILE |\
#awk '
#  /^>/{print}
#  !/^>/{
#    rec=""
#    while(match($0,/(A{3,}|T{3,}|C{3,}|G{3,})/)){
#      rec = rec substr($0,1,RSTART+1) tolower(substr($0,RSTART+2,RLENGTH-2))
#      $0=substr($0,RSTART+RLENGTH)
#    }
#    print rec $0
#  }' > $OUT_DIR/${CONSENSUS_NAME}_masked.fa 
