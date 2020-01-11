# longread-UMI-pipeline
A collection of scripts for processing longread UMI data.

<b> Pipeline version reference: </b> \
Karst, S. M., Ziels, R. M., Kirkegaard, R. H., & Albertsen, M. (2019). Enabling high-accuracy long-read amplicon sequences using unique molecular identifiers and Nanopore sequencing. bioRxiv, 645903.
https://www.biorxiv.org/content/10.1101/645903v2

## Conda installation

1. Requirements/Dependencies
   OS tested (Linux 3.10.0, Ubuntu 14.04, Ubuntu 16.04)
  `usearch` >=10
2. Download installer script from terminal \
   `wget https://raw.githubusercontent.com/SorenKarst/longread-UMI-pipeline/v0.2/scripts/install_conda.sh`
3. Run installation script from terminal and follow instructions \
   `bash ./install_conda.sh` 

## Manual install

1. Requirements/Dependencies \
   OS tested (Linux 3.10.0, Ubuntu 14.04, Ubuntu 16.04) \
   See `scripts/longread-UMI-pipeline_version_dump.txt`
2. Clone from github in terminal \
   `git clone https://github.com/SorenKarst/longread-UMI-pipeline/tree/v0.2`
3. Make bash scripts executable \
   `find ./longread-UMI-pipeline -name "*.sh" -exec chmod +x {} \;`
4. Install dependencies \
   See `./longread-UMI-pipeline/scripts/install_dependencies.sh` for inspiration.
5. Change paths to dependencies \
   Modify `./longread-UMI-pipeline/scripts/dependencies.sh` in a texteditor.
6. Customize porechop adaptors.py to be able to detect custom primers \
   Replace current `adapters.py` with `./longread-UMI-pipeline/scripts/adapters.py`

## Usage

### Test data
1. Open a terminal anywhere and run:  
  `longread-UMI-pipeline -h` or `/path/to/longread-UMI-pipeline -h`
2. Test longread-UMI-pipeline on test data:  
   Go to /path/to/longread-UMI-pipeline/test_data  
   Open a terminal in the directory and run `longread-UMI-pipeline -d test_reads.fq -s 10 -c 30 -t 1`

### Nanopore R9.4.1 data from Zymo mock
1. Create a working directory and download the Zymo mock fastq data from terminal:  
   `wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR333/003/ERR3336963/ERR3336963_1.fastq.gz; gunzip -c ERR3336963_1.fastq.gz > reads.fq`  
2. Run pipeline:  
  `longread-UMI-pipeline -d reads.fq -s 1000000 -c 30 -t <Number-of-threads>`
3. Run data QC:  
  `longread-UMI-mockanalysis <Number-of-threads>`
