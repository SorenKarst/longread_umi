#!/bin/bash
# DESCRIPTION
#    Paths to dependencies for longread-UMI-pipeline 
#
# IMPLEMENTATION
#    author	Søren Karst (sorenkarst@gmail.com)
#               Ryans Ziels (ziels@mail.ubc.ca)
#    version	0.1.0
#    license	GNU General Public License

# Program paths
export SEQTK=/space/sharedbin/bin/seqtk # https://github.com/lh3/seqtk
export GNUPARALLEL=/space/users/smk/bin/parallel #(wget pi.dk/3 -qO - ||  curl pi.dk/3/) | bash
export USEARCH=/space/sharedbin/bin/usearch11 #https://drive5.com/usearch/download.html
export RACON=/space/users/smk/bin/racon #https://github.com/isovic/racon
export MINIMAP2=/space/users/smk/Software/minimap2-2.16/minimap2 #https://github.com/lh3/minimap2
export GAWK=/usr/bin/gawk 
export MEDAKA_ENV_START='. /space/users/smk/Software/medaka-0.7.0/bin/activate' # https://github.com/nanoporetech/medaka.
# Our implementation of Medaka requires installation in virtual environment as recommended by the developers.
export MEDAKA_ENV_STOP='deactivate;'
export MEDAKA_MODEL=r941_min_high # Specifies which medaka model to use. 
export CUTADAPT=/space/users/smk/bin/cutadapt #https://cutadapt.readthedocs.io/en/stable/installation.html
export PORECHOP_UMI=/space/users/smk/bin/porechop_umi #https://github.com/rrwick/Porechop. We have multiple installations of porechop and renamed one to porechop_umi for convenience.
export FILTLONG=/space/users/smk/Software/Filtlong-0.2.0/bin/filtlong #https://github.com/rrwick/Filtlong
export BWA=/space/users/smk/Software/bwa-0.7.17/bwa #https://github.com/lh3/bwa
export SAMTOOLS=/space/users/smk/Software/samtools_1.9/bin/samtools #http://www.htslib.org/download/
export BCFTOOLS=/space/users/smk/Software/bcftools_1.9/bin/bcftools #http://www.htslib.org/download/

# Scripts paths
export UMI_BINNING=$PIPELINE_PATH/scripts/umi_binning.sh
export CONSENSUS_SRACON=$PIPELINE_PATH/scripts/consensus_seed-racon.sh
export POLISH_MEDAKA=$PIPELINE_PATH/scripts/polish_medaka.sh
export TRIM_AMPLICON=$PIPELINE_PATH/scripts/trim_amplicon.sh
export VARIANTS=$PIPELINE_PATH/scripts/variants.sh
export REF=$PIPELINE_PATH/scripts/zymo-ref-uniq_2019-03-15.fa
export REF_VENDOR=$PIPELINE_PATH/scripts/zymo-ref-uniq_vendor.fa

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
