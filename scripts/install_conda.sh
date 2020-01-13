#!/bin/bash
# DESCRIPTION
#    Install longread_umi as conda environment.
#
# IMPLEMENTATION
#    author	Søren Karst (sorenkarst@gmail.com)
#           Ryan Ziels (ziels@mail.ubc.ca)
#    license	GNU General Public License


# Check conda installation ----------------------------------------------------
if [[ -z $(which conda) ]]; then
  # Ask to install
  read \
    -n 1 \
    -p "Conda not found. Install miniconda3 (y/n)? " \
    ASK_CONDA_INSTALL    
  
  if [ "$ASK_CONDA_INSTALL" != "y" ]; then
    echo ""
    echo "Quitting installation script..."
    echo ""
    exit 1 
  else
    # Install conda
    [ -f Miniconda3-latest-Linux-x86_64.sh ] ||\
      wget "https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh"
    bash ./Miniconda3-latest-Linux-x86_64.sh   
  fi
else
  echo ""
  echo "Conda found"
  echo "version: $(conda -V)"
  echo ""
fi

# Install longread-UMI conda env ----------------------------------------------
echo ""
echo "Installing longread_umi conda environment.."
echo ""

# Define conda env yml
echo "name: longread_umi
channels:
- conda-forge
- bioconda
- defaults
dependencies:
- seqtk=1.3
- parallel=20191122
- racon=1.4.3
- minimap2=2.17
- medaka=0.11.2
- gawk=4.1.3
- cutadapt=2.7
- filtlong=0.2.0
- bwa=0.7.17
- samtools=1.9
- bcftools=1.9
- git
" > ./longread_umi.yml

# Install conda env
conda env create -f ./longread_umi.yml

eval "$(conda shell.bash hook)"
conda activate longread_umi

# Install porechop
WDIR=$PWD
git clone \
  https://github.com/rrwick/Porechop.git \
  $CONDA_PREFIX/Porechop

cd $CONDA_PREFIX/Porechop
python3 setup.py install
cd $WDIR  

# Download longread-UMI from git
git clone \
  https://github.com/SorenKarst/longread-UMI-pipeline.git \
  $CONDA_PREFIX/longread_umi

# Modify adapters.py
\cp \
  $CONDA_PREFIX/longread_umi/scripts/adapters.py \
  $CONDA_PREFIX/lib/python3.6/site-packages/porechop/adapters.py

# Create links to pipeline
find \
  $CONDA_PREFIX/longread_umi/ \
  -name "*.sh" \
  -exec chmod +x {} \;
  
ln -s \
  $CONDA_PREFIX/longread_umi/longread_umi.sh \
  $CONDA_PREFIX/bin/longread_umi
  
  
# Create link to usearch installation
read \
  -p "Type path to usearch excutable and press enter:  " \
  USEARCH_PATH
  
chmod +x $USEARCH_PATH

ln -s \
  $USEARCH_PATH \
  $CONDA_PREFIX/bin/usearch  
  
# Check installation
if [[ -z $(which longread_umi) ]]; then
  echo ""
  echo "Can't locate longread_umi"
  echo "longread_umi installation failed..."
  echo ""
else
  echo ""
  echo "longread_umi installation success..."
  echo ""
  echo "Path to conda environment: $CONDA_PREFIX"
  echo "Path to pipeline files: $CONDA_PREFIX/longread_umi"
  echo ""
  echo "Run pipeline test:"
  echo ""
  echo "conda activate longread_umi" 
  echo "longread_umi nanopore_pipeline \
-d $CONDA_PREFIX/longread_umi/test_data/test_reads.fq \
-o test \
-v 30 \
-w rrna_operon \
-t 1 \
-q r941_min_high_g330"
  echo "conda deactivate"  
  echo ""  
fi

conda deactivate

# Cleanup
if [ -f Miniconda3-latest-Linux-x86_64.sh  ]; then 
  rm -f ./Miniconda3-latest-Linux-x86_64.sh
fi
if [ -f install_conda.sh  ]; then 
  rm -f ./install_conda.sh
fi
if [ -f longread_umi.yml  ]; then 
  rm -f ./longread_umi.yml
fi
