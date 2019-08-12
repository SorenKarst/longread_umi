#!/bin/bash
# DESCRIPTION
#    Paths to dependencies for longread-UMI-pipeline 
#
# IMPLEMENTATION
#    author	Søren Karst (sorenkarst@gmail.com)
#               Ryans Ziels (ziels@mail.ubc.ca)
#    license	GNU General Public License

# Program paths

export SEQTK=/space/sharedbin/bin/seqtk
export GNUPARALLEL=/space/users/smk/bin/parallel
export RACON=/space/users/smk/bin/racon
export MINIMAP2=/space/users/smk/bin/minimap2-2.5
export GAWK=/usr/bin/gawk
export SAMTOOLS=/space/users/smk/Software/samtools_1.9/bin/samtools
export BCFTOOLS=/space/users/smk/Software/bcftools_1.9/bin/bcftools
export MEDAKA_ENV_START='module load Python; . /space/users/smk/Software/medaka/bin/activate'
export MEDAKA_ENV_STOP='deactivate; module purge'
export MEDAKA_MODEL=r941_min_high
export CUTADAPT=/space/users/smk/.local/bin/cutadapt
export PORECHOP_UMI=/space/users/smk/bin/porechop_umi
export FILTLONG=/space/users/smk/Software/Filtlong-0.2.0/bin/filtlong
export BWA=/usr/bin/bwa
export USEARCH=/space/sharedbin/bin/usearch11

# longread_umi paths
export REF_CURATED=$LONGREAD_UMI_PATH/scripts/zymo-ref-uniq_2019-03-15.fa
export REF_VENDOR=$LONGREAD_UMI_PATH/scripts/zymo-ref-uniq_vendor.fa
export NANOPORE_BARCODES=$LONGREAD_UMI_PATH/scripts/nanopore_barcodes.csv

# Version dump
longread_umi_version_dump (){
  local OUT=${1:-./longread_umi_version_dump.txt}

  echo "Script start: $(date +%Y-%m-%d-%T)"  >> $OUT
  echo "Software Version:" >> $OUT
  echo "longread_umi - $(git --git-dir ${LONGREAD_UMI_PATH}/.git describe --tag)" >> $OUT
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
