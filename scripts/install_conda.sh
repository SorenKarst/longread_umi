#!/bin/bash
# DESCRIPTION
#    Install longread-UMI-pipeline as conda environment.
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
    -p "Conda not found. Install miniconda (y/n)? " \
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
echo "Installing longread-UMI conda environment.."
echo ""

# Define conda env yml
echo "name: longread-UMI
channels:
- conda-forge
- bioconda
- defaults
dependencies:
- seqtk=1.3
- parallel=20190522
- racon=1.3.3
- minimap2=2.11
- medaka=0.7.1
- gawk=4.1.3
- cutadapt=2.3
- porechop=0.2.3_seqan2.1.1
- filtlong=0.2.0
- bwa=0.7.17
- samtools=1.9
- bcftools=1.9
- git
" > ./longread-UMI.yml

# Install conda env
conda env create -f ./longread-UMI.yml

# Download longread-UMI from git
eval "$(conda shell.bash hook)"
conda activate longread-UMI
git clone \
  -b v0.2 \
  https://github.com/SorenKarst/longread-UMI-pipeline.git \
  $CONDA_PREFIX/longread-UMI

# Modify adapters.py
\cp \
  $CONDA_PREFIX/longread-UMI/scripts/adapters.py \
  $CONDA_PREFIX/lib/python3.6/site-packages/porechop/adapters.py

# Create links to pipeline
find \
  $CONDA_PREFIX/longread-UMI/ \
  -name "*.sh" \
  -exec chmod +x {} \;
  
ln -s \
  $CONDA_PREFIX/longread-UMI/longread_UMI_pipeline.sh \
  $CONDA_PREFIX/bin/longread-UMI-pipeline
  
ln -s \
  $CONDA_PREFIX/longread-UMI/longread_UMI_mockanalysis.sh \
  $CONDA_PREFIX/bin/longread-UMI-mockanalysis
  
# Create link to usearch installation
read \
  -p "Type path to usearch excutable and press enter:  " \
  USEARCH_PATH
  
chmod +x $USEARCH_PATH

ln -s \
  $USEARCH_PATH \
  $CONDA_PREFIX/bin/usearch  
  
# Check installation
if [[ -z $(which longread-UMI-pipeline) ]]; then
  echo ""
  echo "Can't locate longread-UMI-pipeline."
  echo "longread-UMI installation failed..."
  echo ""
else
  echo ""
  echo "longread-UMI installation success..."
  echo ""
  echo "Path to conda environment: $CONDA_PREFIX"
  echo "Path to pipeline files: $CONDA_PREFIX/longread-UMI"
  echo ""
  echo "Run pipeline test:"
  echo ""
  echo "mkdir test; cd test"
  echo "conda activate longread-UMI" 
  echo "longread-UMI-pipeline -d $CONDA_PREFIX/longread-UMI/test_data/test_reads.fq -s 10 -c 30 -t 1"
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
if [ -f longread-UMI.yml  ]; then 
  rm -f ./longread-UMI.yml
fi
