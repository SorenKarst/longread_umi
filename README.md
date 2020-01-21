# longread_umi 

A collection of scripts for processing longread UMI data.

**Table of contents**
- [Installation](#installation)
- [Quick start](#quick-start)
- [Usage](#usage)
- [Data and examples](#data-and-examples)

**Citation**  
SM Karst, RM Ziels, RH Kirkegaard, EA Sørensen, D. McDonald, Q Zhu, R Knight, & M Albertsen. (2019). Enabling high-accuracy long-read amplicon sequences using unique molecular identifiers with Nanopore or PacBio sequencing. bioRxiv, 645903.
https://www.biorxiv.org/content/10.1101/645903v3

## Installation

### Conda

1. Requirements/Dependencies \
   OS tested (Linux 3.10.0, Ubuntu 14.04, Ubuntu 16.04) \
  `usearch` >=10
2. Download installer script from terminal \
   `wget https://raw.githubusercontent.com/SorenKarst/longread_umi/master/scripts/install_conda.sh`
3. Run installation script from terminal and follow instructions \
   `bash ./install_conda.sh` 
4. Initiate conda and refresh terminal before using pipeline. \
   `conda init; source ~/.bashrc`  
5. Activate and deactivate conda environment
   ```
   conda activate longread_umi
   ...
   conda deactivat
   
   ```

### Manual

1. Requirements/Dependencies \
   OS tested (Linux 3.10.0, Ubuntu 14.04, Ubuntu 16.04) \
   See `scripts/longread-UMI-pipeline_version_dump.txt`
2. Clone from github in terminal \
   `git clone https://github.com/SorenKarst/longread_umi.git`
3. Make bash scripts executable \
   `find ./longread_umi -name "*.sh" -exec chmod +x {} \;`
4. Install dependencies \
   See `./longread_umi/scripts/install_dependencies.sh` for inspiration.
5. Change paths to dependencies \
   Modify `./longread_umi/scripts/dependencies.sh` in a texteditor.
6. Customize porechop adaptors.py to be able to detect custom primers \
   Replace current `adapters.py` with `./longread_umi/scripts/adapters.py`

## Quick start

### Test data
1. Test the initialization command in terminal  
  `longread_umi -h` or `/path/to/longread_umi.sh -h`
2. Test the nanopore_pipeline in terminal  
  `longread_umi nanopore_pipeline -h` or `/path/to/longread_umi.sh nanopore_pipeline -h`
3. Test longread_umi nanopore_pipeline and qc_pipeline on test data:  
   Go to /path/to/longread-UMI-pipeline/test_data and open a terminal in the directory.
4. Run nanopore pipeline
   ```
   longread_umi nanopore_pipeline \
     -d test_reads.fq \
     -v 30 \
     -o test \
     -s 90 \
     -e 90 \
     -m 3500 \
     -M 6000 \
     -f CAAGCAGAAGACGGCATACGAGAT \
     -F AGRGTTYGATYMTGGCTCAG \
     -r AATGATACGGCGACCACCGAGATC \
     -R CGACATCGAGGTGCCAAAC \
     -c 3 \
     -p 1 \
     -q r941_min_high_g330 \
     -t 1
   ```
5. Run qc pipeline
   ```
   longread_umi qc_pipeline \
     -d test_reads.fq \
     -c test/consensus_raconx3_medakax1.fa \
     -r zymo_curated \
     -t 1 \
     -o test/qc
   ```

### Zymomock rRNA operon data
1. Download the Zymo mock Nanopore R9.4.1 fastq data and decompress
   ```
   wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR333/003/ERR3336963/ERR3336963_1.fastq.gz 
   gunzip -c ERR3336963_1.fastq.gz > reads.fq
   ```
2. Run nanopore pipeline
   ```
   longread_umi nanopore_pipeline \  
     -d reads.fq \
     -o analysis \
     -v 30  \
     -s 90 \
     -e 90 \
     -m 3500 \
     -M 6000 \
     -f CAAGCAGAAGACGGCATACGAGAT \
     -F AGRGTTYGATYMTGGCTCAG \
     -r AATGATACGGCGACCACCGAGATC \
     -R CGACATCGAGGTGCCAAAC \
     -c 3 \
     -p 1 \
     -q r941_min_high_330  \
     -t <Number-of-threads>
   ```
5. Open a terminal in the directory and run
   ```
   longread_umi qc_pipeline \
     -d "umi_binning/trim/reads_tf.fq;reads.fq" \
     -c "consensus_raconx3_medakax1.fa;variants.fa" \
     -r "zymo_curated;zymo_vendor;variants.fa" \
     -t <Number-of-threads> 
   ```

## Usage

...

## Data and examples

- [Example data](docs/DATA.md)
