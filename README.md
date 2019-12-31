# longread-UMI-pipeline
A collection of scripts for processing longread UMI data.

## Requirements
1. Tested on Ubuntu 14.04 (Yeah, we know...)
2. Dependencies: See scripts/dependencies.sh and scripts/longread-UMI-pipeline_version_dump.txt

## Automatic installation
1. Naive semi-automatic installation of pipeline and dependencies. The script will overwrite
   folders and files in the working directory. USE AT OWN RISK!!
2. Go to a folder where you want the longread UMI pipeline and dependencies installed.
3. Open a terminal in the installation folder and download the installation script:  
  `wget https://raw.githubusercontent.com/SorenKarst/longread-UMI-pipeline/master/scripts/install_dependencies.sh`
4. Option A (Recommended): Open script in a text editor and copy installation commands for missing dependencies to
   terminal one by one.
4. Option B: Install pipeline and dependencies automatically by excuting script `bash install_dependencies.sh`

## Manual installation

### Clone from github
1. Go to desired installation directory, open a terminal and run:  
   `git clone https://github.com/SorenKarst/longread-UMI-pipeline`

### Make bash scripts executable
1. Go to longread-UMI-pipeline directory, open a terminal and run:  
   `find . -name "*.sh" -exec chmod +x {} \;`

### Create symlink to run longread-UMI-pipeline and longread-UMI-mockanalysis from terminal
1. Create symlink in ~/bin by opening a terminal and run:  
   `mkdir -p ~/bin`  
   `ln -s /path/to/longread-UMI-pipeline/longread_UMI_pipeline.sh ~/bin/longread-UMI-pipeline`  
   `ln -s /path/to/longread-UMI-pipeline/longread_UMI_mockanalysis.sh ~/bin/longread-UMI-mockanalysis`

### Change paths to dependencies
1. Open /path/to/longread-UMI-pipeline/scripts/dependencies.sh in a texteditor.
2. Change all paths under "Program paths" to reflect installation paths on your system.
3. If unsure of the paths try to type `which <function>` in the terminal. I.e. `which racon`.
4. Install any missing dependencies.

### Customize porechop adaptors.py to be able to detect custom primers
1. We recommend to make a seperate installation of porechop to use with the longread-UMI-pipeline.
2. Go to path/to/porechop/porechop/
3. Backup current adapters.py.
4. Replace current adapters.py with path/to/longread-UMI-pipeline/scripts/adapters.py.

### Test scripts
1. Open a terminal anywhere and run:  
  `longread-UMI-pipeline -h` or `/path/to/longread-UMI-pipeline -h`
2. Test longread-UMI-pipeline on test data:  
   Go to /path/to/longread-UMI-pipeline/test_data  
   Open a terminal in the directory and run `longread-UMI-pipeline -d test_reads.fq -s 10 -c 30 -t 1`

### Run pipeline on Zymo mock data
1. Create a working directory, open a terminal, download the Zymo mock fastq data and decompress:  
   `wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR333/003/ERR3336963/ERR3336963_1.fastq.gz; gunzip -c ERR3336963_1.fastq.gz > reads.fq`  
2. Open a terminal in the directory and run:  
  `longread-UMI-pipeline -d reads.fq -s 1000000 -c 30 -t <Number-of-threads>`
3. Open a terminal in the directory and run:  
  `longread-UMI-mockanalysis <Number-of-threads>` 

# American Gut sequence comparison
A collection of scripts and notebooks for comparing operon data to 16S V4 data from the American Gut

## Requirements

1. Tested on Centos 7
2. Dependencies: QIIME 2 2019.10, redbiom >= 0.3.0, BLAST 2.7.1, RNAmmer 1.2, VSEARCH 2.7.0, Web of Life 86k reference.

## Installation

Please see the install guides for [QIIME 2](https://docs.qiime2.org/2019.10/install/) and [redbiom](https://github.com/biocore/redbiom#client). Both can be installed through Conda, however the specific installation instructions for QIIME 2 depend on whether the environment is Linux or OSX.

BLAST and VSEARCH are installed with QIIME 2.

For RNAmmer, please see the installation [instructions](http://www.cbs.dtu.dk/services/doc/rnammer-1.2.readme).

Instructions for obtaining the Web of Life database are [here](https://biocore.github.io/wol/data/genomes/). 

## Taxonomic assignments

RNAmmer (see `scripts/rnammer.sh`) must first be run on per-sample. The per-gene sequences are then extracted (see `scripts/extract_rrna.py`). The individual gene sequences must then be  blast'd (see `scripts/blastn.sh`) against the Web of Life database. Taxonomy can be assigned (see `scripts/assign_taxonomy.py`). The final comparison is performed in a Jupyter Notebook (see `notebooks/analysis.ipynb`)

## Genus comparison

The genus comparison between full length 16S and 16S V4 data assumes the 16S sequences were etracted already (see "Taxonomic assignments"). This repository already contains the extracted 16S sequences (see `notebooks/annotations/q2_seqs.qza`) and the subsequent QIIME 2 processing steps (see `notebooks/annotations`). However, those steps can be rerun with `notebooks/annotations/q2-over-16S.sh`. The methods for obtaining the AGP 16S V4 data, and subsequent comparison, are performed in the Jupyter Notebook `notebooks/genus_16Sv4_vs_16S_comparison.ipynb`).
