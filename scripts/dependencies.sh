
#!/bin/bash
# DESCRIPTION
#    Paths to dependencies for longread-UMI-pipeline 
#
# IMPLEMENTATION
#    author	Søren Karst (sorenkarst@gmail.com)
#               Ryans Ziels (ziels@mail.ubc.ca)
#    license	GNU General Public License

# Program paths

export SEQTK=seqtk
export GNUPARALLEL=parallel
export RACON=racon
export MINIMAP2=minimap2
export GAWK=gawk
export SAMTOOLS=samtools
export BCFTOOLS=bcftools
#export MEDAKA_ENV_START='. /space/users/smk/software/medaka/bin/activate'
#export MEDAKA_ENV_STOP='deactivate'
export MEDAKA_MODEL=r941_min_high
export CUTADAPT=cutadapt
export PORECHOP_UMI=porechop
export FILTLONG=filtlong
export BWA=bwa
export USEARCH=usearch_path

# Scripts paths
export UMI_BINNING=$PIPELINE_PATH/scripts/umi_binning.sh
export CONSENSUS_SRACON=$PIPELINE_PATH/scripts/consensus_seed-racon.sh
export POLISH_MEDAKA=$PIPELINE_PATH/scripts/polish_medaka.sh
export TRIM_AMPLICON=$PIPELINE_PATH/scripts/trim_amplicon.sh
export VARIANTS=$PIPELINE_PATH/scripts/variants.sh
export REF=$PIPELINE_PATH/scripts/zymo-ref-uniq_2019-03-15.fa
export REF_VENDOR=$PIPELINE_PATH/scripts/zymo-ref-uniq_vendor.fa
export NANOPORE_BARCODES=$PIPELINE_PATH/scripts/nanopore_barcodes.csv

# Version dump
ncec_version_dump (){
  OUT=${1:-./longread-UMI-pipeline_version_dump.txt}

  echo "Script start: $(date +%Y-%m-%d-%T)"  >> $OUT
  echo "Software Version:" >> $OUT
  echo "seqtk - $($SEQTK 2>&1 >/dev/null | grep 'Version')" >> $OUT 
  echo "Parallel - $($GNUPARALLEL --version | head -n 1)" >> $OUT 
  echo "Usearch - $($USEARCH --version)" >> $OUT 
  echo "Racon - $($RACON --version)" >> $OUT
  echo "Minimap2 - $($MINIMAP2 --version)" >> $OUT
  echo "medaka - $(eval $MEDAKA_ENV_START; medaka --version | cut -d" " -f2; eval $MEDAKA_ENV_STOP)"  >> $OUT
  echo "medaka model - ${MEDAKA_MODEL##*/}"  >> $OUT
  echo "Gawk - $($GAWK --version | head -n 1)" >> $OUT 
  echo "Cutadapt - $($CUTADAPT --version | head -n 1)" >> $OUT 
  echo "Porechop - $($PORECHOP_UMI --version) + add UMI adaptors to adaptors.py" >> $OUT 
  echo "Filtlong - $($FILTLONG --version)" >> $OUT
  echo "BWA - $($BWA 2>&1 >/dev/null | grep 'Version')" >> $OUT
  echo "Samtools - $($SAMTOOLS 2>&1 >/dev/null | grep 'Version')" >> $OUT
  echo "Bcftools - $($BCFTOOLS --version | head -n 1)" >> $OUT
}

### Version dump
# source dependencies.sh
# ncec_version_dump
