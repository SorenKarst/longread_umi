#!/bin/bash
# DESCRIPTION
#    Instructions for installing dependencies for  
#    longread-UMI-pipeline.
#
#    USE AT OWN RISK!!!
#    Will likely require system specific modifications to work.
#    Preferably run manually one dependency at a time.
#
#    Move script to directory where you want dependencies installed.
#    Run in terminal bash install_dependencies.sh
#
# IMPLEMENTATION
#    author   Søren Karst (sorenkarst@gmail.com)
#             Ryans Ziels (ziels@mail.ubc.ca)
#    license  GNU General Public License

### Terminal input
BRANCH=${1:-master}

# Store software dir path
SOFTWARE_DIR=$PWD

### Check for installation depenencies
if command -v conda >/dev/null 2>&1 ; then
    echo "conda found"
    echo "version: $(conda -V)"
else
    echo "conda not found. Install conda and re-run install script:
	# Miniconda3 install example
	wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
	sh Miniconda3-latest-Linux-x86_64.sh -p $PWD/miniconda3
	source ~/.bashrc
	conda config --add channels defaults
	conda config --add channels bioconda
	conda config --add channels conda-forge
	conda config --set auto_activate_base false
	source ~/.bashrc"
	exit 0
fi


### Create file with paths

echo '' > ./longread-UMI-pipeline_paths.txt

# Make ~/bin if it doesn't exist
mkdir -p ~/bin

# pip
curl "https://bootstrap.pypa.io/get-pip.py" -o "get-pip.py"
python3 get-pip.py --user
rm ./get-pip.py

# Python virtual environment
# python3 -m pip install virtualenv --user

# Cmake
git clone https://github.com/scivision/cmake-utils.git -b v1.4.0.0;
cd cmake-utils
python3 cmake_setup.py 3.15.2 \
  --install_path $SOFTWARE_DIR/cmake
cd ..
rm -rf ./cmake-utils

### Install longread-UMI-pipeline
git clone https://github.com/SorenKarst/longread-UMI-pipeline -b $BRANCH
cd ./longread-UMI-pipeline
find . -name "*.sh" -exec chmod +x {} \;
cd ..
ln -s $SOFTWARE_DIR/longread-UMI-pipeline/longread_umi.sh ~/bin/longread_umi

### Install dependencies automaticly

# Seqtk
git clone https://github.com/lh3/seqtk.git;
cd seqtk; make
cd ..
echo "export SEQTK=$SOFTWARE_DIR/seqtk/seqtk" >> ./longread-UMI-pipeline_paths.txt

# GNU Parallel
(wget pi.dk/3 -qO - ||  curl pi.dk/3/) | bash
echo "export GNUPARALLEL=$(which parallel)" >> ./longread-UMI-pipeline_paths.txt
rm -rf ./parallel*

# Racon
git clone --recursive https://github.com/isovic/racon.git racon
cd racon
mkdir build
cd build
$SOFTWARE_DIR/cmake/cmake*/bin/cmake -DCMAKE_BUILD_TYPE=Release ..
make
cd ../..
echo "export RACON=$SOFTWARE_DIR/racon/build/bin/racon" >> ./longread-UMI-pipeline_paths.txt

# Minimap2
git clone https://github.com/lh3/minimap2
cd minimap2 && make
cd ..
ln -s $SOFTWARE_DIR/minimap2/minimap2 ~/bin/minimap2
echo "export MINIMAP2=$(which minimap2)" >> ./longread-UMI-pipeline_paths.txt

# Gawk
# Check presence by:
# which gawk
# If not present install
echo "export GAWK=$(which gawk)" >> ./longread-UMI-pipeline_paths.txt

#Samtools
wget "https://github.com/samtools/samtools/releases/download/1.9/samtools-1.9.tar.bz2"
tar -xjf ./samtools-1.9.tar.bz2
cd samtools-1.9
./configure \
  --prefix=$SOFTWARE_DIR/samtools_1.9 \
  --disable-bz2 \
  --disable-lzma
make
make install
cd ..
rm -rf ./samtools-1.9*
ln -s $SOFTWARE_DIR/samtools_1.9/bin/samtools ~/bin/samtools
echo "export SAMTOOLS=$(which samtools)" >> ./longread-UMI-pipeline_paths.txt

#Bcftools
wget "https://github.com/samtools/bcftools/releases/download/1.9/bcftools-1.9.tar.bz2"
tar -xjf ./bcftools-1.9.tar.bz2
cd bcftools-1.9 
./configure \
  --prefix=$SOFTWARE_DIR/bcftools_1.9 \
  --disable-bz2 \
  --disable-lzma
make
make install
cd ..
rm -r ./bcftools-1.9*
ln -s $SOFTWARE_DIR/bcftools_1.9/bin/bcftools ~/bin/bcftools
echo "export BCFTOOLS=$(which bcftools)" >> ./longread-UMI-pipeline_paths.txt

#Htslib
wget "https://github.com/samtools/htslib/releases/download/1.9/htslib-1.9.tar.bz2"
tar -xjf ./htslib-1.9.tar.bz2
cd htslib-1.9 
./configure \
  --prefix=$SOFTWARE_DIR/htslib_1.9 \
  --disable-lzma \
  --disable-bz2 
make
make install
cd ..
rm -r ./htslib-1.9*
ln -s $SOFTWARE_DIR/htslib_1.9/bin/tabix ~/bin/tabix
ln -s $SOFTWARE_DIR/htslib_1.9/bin/bgzip ~/bin/bgzip

# Medaka
conda create -c bioconda -n medaka medaka
echo "export MEDAKA_ENV_START='eval \"\$(conda shell.bash hook)\"; conda activate medaka'" >> ./longread-UMI-pipeline_paths.txt
echo "export MEDAKA_ENV_STOP='conda deactivate'" >> ./longread-UMI-pipeline_paths.txt

# cutadapt
pip3 install --user --upgrade cutadapt
echo "export CUTADAPT=$(which cutadapt)" >> ./longread-UMI-pipeline_paths.txt

# Porechop
git clone https://github.com/rrwick/Porechop.git
cd Porechop
make
cd ..
echo "export PORECHOP_UMI=$SOFTWARE_DIR/Porechop/porechop-runner.py" >> ./longread-UMI-pipeline_paths.txt
mv $SOFTWARE_DIR/Porechop/porechop/adapters.py $SOFTWARE_DIR/Porechop/porechop/adapters_original.py
cp $SOFTWARE_DIR/longread-UMI-pipeline/scripts/adapters.py $SOFTWARE_DIR/Porechop/porechop/

# Filtlong
git clone https://github.com/rrwick/Filtlong.git
cd Filtlong
make -j
cd ..
echo "export FILTLONG=$SOFTWARE_DIR/Filtlong/bin/filtlong" >> ./longread-UMI-pipeline_paths.txt

#BWA
git clone https://github.com/lh3/bwa.git
cd bwa; make
cd ..
echo "export BWA=$SOFTWARE_DIR/bwa/bwa" >> ./longread-UMI-pipeline_paths.txt

### Install dependencies manually

# Usearch
mkdir usearch
echo ""
echo "Download usearch from https://drive5.com/usearch/download.html and place in usearch folder"
echo ""
read -rsp $'Press any key to continue...\n' -n1 key
chmod +x $SOFTWARE_DIR/usearch/usearch*
echo "export USEARCH=$(find $SOFTWARE_DIR/usearch/ -type f -name "usearch*")" >> ./longread-UMI-pipeline_paths.txt

### Add depency paths to dependency.sh

echo "" >> ./longread-UMI-pipeline_paths.txt
LEAD='^# Program paths$'
TAIL='^# longread_umi paths'

sed -i \
  -e "/$LEAD/,/$TAIL/{ /$LEAD/{p; r ./longread-UMI-pipeline_paths.txt
        }; /$TAIL/p; d }"  $SOFTWARE_DIR/longread-UMI-pipeline/scripts/dependencies.sh


### Test longread-UMI-pipeline
cd longread-UMI-pipeline/test_data
longread_umi nanopore_pipeline -d test_reads.fq -o . -v 30 -w rrna_operon -t 1 -q r941_min_high_g303
longread_umi qc_pipeline -d test_reads.fq -c consensus_raconx3_medakax1.fa -r zymo_curated -t 1
