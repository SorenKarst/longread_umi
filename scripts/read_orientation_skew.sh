#!/bin/bash

# DESCRIPTION
#    Script for generating UMI bins with skewed +/- read orientation. Part of 
#    longread_umi.
#    
# IMPLEMENTATION
#    author   Søren Karst (sorenkarst@gmail.com)
#             Ryan Ziels (ziels@mail.ubc.ca)
#    license  GNU General Public License

### Source commands and subscripts -------------------------------------
. $LONGREAD_UMI_PATH/scripts/dependencies.sh # Path to dependencies script

### info ---------------------------------------------------------------
USAGE="
-- longread_umi read_orientation_skew: Script for simulating UMI bins with 
   skewed +/- read orientation ratios.

usage: $(basename "$0" .sh) [-h] (-b dir -d file -o dir -r integer range -R double range )
(-u value -U value -O value -S value -t value -T value -s value) 

where:
    -h  Show this help text.
    -b  UMI binning directory [Default: umi_binning] containing output from longread_umi umi_binning.
    -d  Reads in fastq format.
    -o  Output directory.
    -r  Bin size range to simulate +/- skew in. Only bins with 
        that contain n(negative strand reads) >= range_max and 
        n(positive strand reads) >= range_max will be used for
        simulation. This ensures the same bins are used for all
        simulated bin sizes. The range should be provided as a string
        of integers separated by ';' eg. '10;20;30'
    -R  Range of fractions of simulated bins that are made up by 
        positive strand reads. The range should be provided as a string
        of fractions separated by ';' eg. '0.0;0.1;0.2;0.3;0.4;0.5;0.6;0.7;0.8;0.9;1.0'
    -u  Discard bins with a mean UMI match error above u.
    -U  Discard bins with a UMI match error standard
        deviation above U.
    -O  Normalize read orientation fraction to 'O' if < 'O' reads are
        either +/- strand orientation. [Default = 0] which is disabled.
    -N  Max number of reads with +/- orientation. [Default = 10000]
    -S  UMI bin size/UMI cluster size cutoff. [Default = 10]
    -t  Number of threads to use.
    -T  Number of threads to use for splitting reads into UMI bins. [Default is same as -t]
    -s  Subsample number of UMI bins. [Default = 1000]
"

### Terminal Arguments ---------------------------------------------------------

# Import user arguments
while getopts ':hzb:d:o:R:r:t:u:U:O:N:S:T:s:' OPTION; do
  case $OPTION in
    h) echo "$USAGE"; exit 1;;
    b) IN_DIR=$OPTARG;;
    d) READ_IN=$OPTARG;;
    o) OUT_DIR=$OPTARG;;
    r) SCREEN_UBS_RANGE=$OPTARG;;
    R) SCRREN_ROR_RANGE=$OPTARG;;
    u) UMI_MATCH_ERROR=$OPTARG;;
    U) UMI_MATCH_ERROR_SD=$OPTARG;;
    O) RO_FRAC=$OPTARG;;
    N) MAX_BIN_SIZE=$OPTARG;;
    S) BIN_CLUSTER_RATIO=$OPTARG;;
    t) THREADS=$OPTARG;;
    T) BIN_THREADS=$OPTARG;;
    s) SUBSAMPLE_BINS=$OPTARG;;
    :) printf "missing argument for -$OPTARG\n" >&2; exit 1;;
    \?) printf "invalid option for -$OPTARG\n" >&2; exit 1;;
  esac
done

# Check missing arguments
MISSING="is missing but required. Exiting."
if [ -z ${IN_DIR+x} ]; then echo "-b $MISSING"; echo ""; echo "$USAGE"; exit 1; fi;
if [ -z ${READ_IN+x} ]; then echo "-d $MISSING"; echo ""; echo "$USAGE"; exit 1; fi;
if [ -z ${OUT_DIR+x} ]; then echo "-o $MISSING"; echo ""; echo "$USAGE"; exit 1; fi;
if [ -z ${SCREEN_UBS_RANGE+x} ]; then echo "-r $MISSING"; echo ""; echo "$USAGE"; exit 1; fi;
if [ -z ${SCRREN_ROR_RANGE+x} ]; then echo "-R $MISSING"; echo ""; echo "$USAGE"; exit 1; fi;
if [ -z ${UMI_MATCH_ERROR+x} ]; then echo "-u $MISSING"; echo ""; echo "$USAGE"; exit 1; fi;
if [ -z ${UMI_MATCH_ERROR_SD+x} ]; then echo "-U $MISSING"; echo ""; echo "$USAGE"; exit 1; fi;
if [ -z ${RO_FRAC+x} ]; then echo "-O is missing. Read orientation filter disabled."; RO_FRAC=0; fi;
if [ -z ${MAX_BIN_SIZE+x} ]; then echo "-N is missing. Defaulting to 10000 +/- reads ."; MAX_BIN_SIZE=10000; fi;
if [ -z ${BIN_CLUSTER_RATIO+x} ]; then echo "-S is missing. Defaulting to 10 ."; BIN_CLUSTER_RATIO=10; fi;
if [ -z ${THREADS+x} ]; then echo "-t is missing. Defaulting to 1 thread."; THREADS=1; fi;
if [ -z ${BIN_THREADS+x} ]; then BIN_THREADS=$THREADS; fi;
if [ -z ${SUBSAMPLE_BINS+x} ]; then echo "-s is missing. Defaulting to 1000 bins."; SUBSAMPLE_BINS=1000; fi;

# UMI binning and filtering
mkdir -p $OUT_DIR/bins

UMI_DIR=$IN_DIR/umi_ref
BINNING_DIR=$IN_DIR/read_binning

$GAWK \
  -v BD="$OUT_DIR" \
  -v UM1="$BINNING_DIR/umi1_map.sam" \
  -v UM2="$BINNING_DIR/umi2_map.sam" \
  -v URC="$UMI_DIR/umi_ref_chimera.txt" \
  -v URD="$UMI_DIR/umi_ref_derivates.txt" \
  -v UME_MATCH_ERROR="$UMI_MATCH_ERROR" \
  -v UME_MATCH_ERROR_SD="$UMI_MATCH_ERROR_SD" \
  -v RO_FRAC="$RO_FRAC" \
  -v MAX_BIN_SIZE="$MAX_BIN_SIZE"  \
  -v BIN_CLUSTER_RATIO="$BIN_CLUSTER_RATIO" \
  -v SCREEN_UBS_RANGE="$SCREEN_UBS_RANGE" \
  -v SCREEN_ROR_RANGE="$SCRREN_ROR_RANGE" \
  -v SUBSAMPLE_BINS="$SUBSAMPLE_BINS" \
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
      "chimera_filter", \
      "derivate_match", \
      "derivate_filter" \
      > BD"/umi_binning_stats.txt"
    for (u in umi_n_raw){
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
        derivate_check[u] \
        > BD"/umi_binning_stats.txt"
    }

    print "[" strftime("%T") "] Print UMI bins with range of read orientation ratios" > "/dev/stderr" 
    # Preparation
    split(SCREEN_ROR_RANGE, srr, ";")
    split(SCREEN_UBS_RANGE, sur, ";")
    for (i in sur){
      if (sur[i]+0 > SUR_MAX+0){
        SUR_MAX=sur[i]+0
      }
    }
    srand() # Update random seed
    # Scan UMIs
    for (u in umi_n){
      if (SS <= SUBSAMPLE_BINS){
        # Format UMI name
        UMI_NAME=u
        gsub(";.*", "", UMI_NAME)
        # Use UMI if enough reads of both orientations exist
        if (umi_ro_plus[u]+0 >= SUR_MAX && umi_ro_neg[u]+0 >= SUR_MAX){
          # Count UMI bins generated
          SS++
          # Reset read array
          delete umi_read_pos
          delete umi_read_neg
          # Extract all reads from each UMI
          for (r in match_umi){
            if (match_umi[r] == u){
              if (match_ro[r] == "+"){
                umi_read_pos[length(umi_read_pos)+1]=r
              } else if (match_ro[r] == "-"){
                umi_read_neg[length(umi_read_neg)+1]=r
              }
            }
          }

          # Create subsampling arrays 
          SS_POS=""
          for (j in umi_read_pos){
            SS_POS=SS_POS";"j
          }
          sub("^;", "", SS_POS)
          split(SS_POS, ss_pos, ";")
          SS_NEG=""
          for (j in umi_read_neg){
            SS_NEG=SS_NEG";"j
          }
          sub("^;", "", SS_NEG)
          split(SS_NEG, ss_neg, ";")

          # Loop over screen UBS range
          for (b in sur){
            # Loop over screen ror range
            for (s in srr){
              # Target number of reads
              N_READS_POS=int(sur[b]*srr[s])
              N_READS_NEG=sur[b]-N_READS_POS

              # Load subsampling array
              for (j in ss_pos){
                ss_pos_tmp[j]=ss_pos[j]
              }
              for (j in ss_neg){
                ss_neg_tmp[j]=ss_neg[j]
              }

              # Subsample reads
              for (i = 1; i <= N_READS_POS; i++){
                # Randomly pick index position and extract read
                SS_POS_I = int(length(ss_pos_tmp)*rand()) + 1
                r = umi_read_pos[ss_pos_tmp[SS_POS_I]]
                # Print read to bin
                print UMI_NAME"_"sur[b]"_"substr(N_READS_POS/sur[b],1,3), r, match_err[r], "strand_pos"
                # Update subsampling array
                delete ss_pos_tmp[SS_POS_I]
                SS_POS_TMP=""
                for (j in ss_pos_tmp){
                  SS_POS_TMP=SS_POS_TMP";"ss_pos_tmp[j]
                }
                sub("^;", "", SS_POS_TMP)
                split(SS_POS_TMP, ss_pos_tmp, ";")
              }
              for (i = 1; i <= N_READS_NEG; i++){
                # Randomly pick index position and extract read
                SS_NEG_I = int(length(ss_neg_tmp)*rand()) + 1
                r = umi_read_neg[ss_neg_tmp[SS_NEG_I]]
                # Print read to bin
                print UMI_NAME"_"sur[b]"_"substr(N_READS_POS/sur[b],1,3), r, match_err[r], "strand_neg"
                # Update subsampling array
                delete ss_neg_tmp[SS_NEG_I]
                SS_NEG_TMP=""
                for (j in ss_neg_tmp){
                  SS_NEG_TMP=SS_NEG_TMP";"ss_neg_tmp[j]
                }
                sub("^;", "", SS_NEG_TMP)
                split(SS_NEG_TMP, ss_neg_tmp, ";")
              }
            }
          }
        }
      }
    }

    # Print to terminal
    print "[" strftime("%T") "] Done." > "/dev/stderr" 
  }
' \
$BINNING_DIR/umi1_map.sam \
$BINNING_DIR/umi2_map.sam \
$UMI_DIR/umi_ref_derivates.txt \
$UMI_DIR/umi_ref_chimera.txt \
> $OUT_DIR/umi_bin_map.txt

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
      bin[$2]=bin[$2]";"$1
      strand[$2]=$4
      next;
    }
    # Stream fastq records and write to bins
    FNR%4==1 {
      # Format read name
      read=substr($1,2)

      # Store read in memory
      READ_RECORD=$0" "strand[read]
      getline; READ_RECORD=READ_RECORD"\n"$0
      getline; READ_RECORD=READ_RECORD"\n"$0
      getline; READ_RECORD=READ_RECORD"\n"$0

      # Allocate reads to bins
      BIN_TMP=bin[read]
      sub("^;", "", BIN_TMP)
      split(BIN_TMP,bin_tmp,";")
      for (i in bin_tmp){
        if ( bin_tmp[i] != "" ){
          print READ_RECORD > out"/"bin_tmp[i]"bins.fastq"
        }
      }
    }
  ' - $READS
}

export -f umi_binning

cut -d " " -f1 $OUT_DIR/umi_bin_map.txt |\
  sort -u |\
  $GNUPARALLEL \
    --env umi_binning \
    -N 4000 \
	-j $THREADS \
	--pipe \
  "
    mkdir $OUT_DIR/bins/{#}
    cat |\
      umi_binning \
        $READ_IN \
        $OUT_DIR/umi_bin_map.txt \
        $OUT_DIR/bins/{#}
  "