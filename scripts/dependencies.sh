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
export SEQTK=/space/sharedbin/bin/seqtk
export GNUPARALLEL=/space/users/smk/bin/parallel
export USEARCH=/space/sharedbin/bin/usearch11
export RACON=/space/users/smk/bin/racon
export MINIMAP2=/space/users/smk/Software/minimap2-2.16/minimap2
export GAWK=/usr/bin/gawk
export MEDAKA_ENV_START='module load Python/3.6.4-foss-2018a;. /space/users/smk/Software/medaka-0.7.0/bin/activate'
export MEDAKA_ENV_STOP='deactivate;module purge'
export MEDAKA_MODEL=/space/users/smk/Software/medaka-0.7.0/lib/python3.6/site-packages/medaka/data/r941_min_high_model.hdf5
export CUTADAPT=/space/users/smk/bin/cutadapt
export PORECHOP_UMI=/space/users/smk/bin/porechop_umi
export FILTLONG=/space/users/smk/Software/Filtlong-0.2.0/bin/filtlong
export BWA=/space/users/smk/Software/bwa-0.7.17/bwa
export SAMTOOLS=/space/users/smk/Software/samtools_1.9/bin/samtools
export REFORMAT=/space/users/smk/Software/bbmap/reformat.sh
export FASTQSTATS=/usr/bin/fastq-stats
export BCFTOOLS=/space/users/smk/Software/bcftools_1.9/bin/bcftools

# Scripts paths
export UMI_BINNING=./scripts/umi_binning.sh
export CONSENSUS_SPADES=./scripts/consensus_spades.sh
export CONSENSUS_SRACON=./scripts/consensus_seed-racon.sh
export POLISH_MEDAKA=./scripts/polish_medaka.sh
export POLISH_RACON=./scripts/polish_racon.sh
export TRIM_AMPLICON=./scripts/trim_amplicon.sh
export VALIDATE=./scripts/validation.sh
export VARIANTS=./scripts/variants.sh

# Version dump
ncec_version_dump (){
  OUT=${1:-./ncec_version_dump_$(date +%Y-%m-%d-%T).txt}

  echo "Software Version:" >> $OUT
  echo "seqtk - $($SEQTK 2>&1 >/dev/null | grep 'Version')" >> $OUT 
  echo "Parallel - $($GNUPARALLEL --version | head -n 1)" >> $OUT 
  echo "Usearch - $($USEARCH --version)" >> $OUT 
  echo "Racon - $($RACON --version)" >> $OUT
  echo "Minimap2 - $($MINIMAP2 --version)" >> $OUT
  echo "medaka - ?"  >> $OUT
  echo "medaka model - ${MEDAKA_MODEL##*/}"  >> $OUT
  echo "Gawk - $($GAWK --version | head -n 1)" >> $OUT 
  echo "Cutadapt - $($CUTADAPT --version | head -n 1)" >> $OUT 
  echo "Porechop - $($PORECHOP_UMI --version) + add UMI adaptors to adaptors.py" >> $OUT 
  echo "Filtlong - $($FILTLONG --version)" >> $OUT
  echo "BWA - $($BWA 2>&1 >/dev/null | grep 'Version')" >> $OUT
  echo "Samtools - $($SAMTOOLS 2>&1 >/dev/null | grep 'Version')" >> $OUT
  echo "bbmap reformat.sh - $($REFORMAT --version 2>&1 | awk 'NR==2')" >> $OUT
  echo "Bcftools - $($BCFTOOLS --version | head -n 1)" >> $OUT
}

### Version dump
# source ncec_dependencies.sh
# ncec_version_dump
