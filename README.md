# longread-UMI-pipeline
A collection of scripts for processing longread UMI data.

## Installation

### Requirements
1. Tested on Ubuntu 14.04 (Yeah, we know...)
2. Dependencies: See scripts/dependencies.sh and scripts/longread-UMI-pipeline_version_dump*.txt

### Clone from github
1. Go to desired installation directory, open a terminal and run:
2. `git clone https://github.com/SorenKarst/longread-UMI-pipeline`

### Make bash scripts executable
1. Go to longread-UMI-pipeline directory, open a terminal and run:
   `find . -name "*.sh" -exec chmod +x {} \;`

### Test scripts
1. Open a terminal anywhere and run: `/path/to/longread-UMI-pipeline -h`

### (Optional) Create symlink to run longread-UMI-pipeline from terminal
1. Check if ~/bin exists. If not create it by opening a terminal and running: 
   `mkdir -p ~/bin`
2. Check if ~/bin is in path `echo $PATH`. If not abort here and just use full path to run the script.
1. Create symlink in bin by openening a terminal and run:
   `ln -s /path/to/longread-UMI-pipeline/longread_UMI_pipeline.sh ~/bin/longread-UMI-pipeline`
2. Test symlink by opening a terminal anywhere and run:
   `longread-UMI-pipeline -h`

### Change paths to dependencies
1. Go to /path/to/longread-UMI-pipeline/scripts and open dependencies.sh in a texteditor.
2. Change all paths under "Program paths" to reflect installation paths on your system.
3. If unsure of the paths try to type `which <function>` in the terminal. I.e. `which racon`.
4. Install any missing dependencies.

### Customize porechop adaptors.py to be able to detect custom primers
1. We recommend to make a seperate installation of porechop to use with the longread-UMI-pipeline.
2. Go to path/to/porechop/porechop/adapters.py
3. Backup current adapters.py.
4. Replace current adapters.py with adapters.py from path/to/longread-UMI-pipeline/scripts.

## Running longread-UMI-pipeline

### Running longread-UMI-pipeline on test data
1. Go to /path/to/longread-UMI-pipeline/test_data
2. Open a terminal in the directory and run `longread-UMI-pipeline -d test_reads.fq -s 10 -c 30 -t 1`
   Or run `./test_longread_umi_pipeline.sh`.

### Run pipeline on Zymo mock data
1. Create a working directory, open a terminal, download the Zymo mock fastq data and decompress:
   ENA: From terminal `wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR333/003/ERR3336963/ERR3336963.fastq.gz`
   Figshare: From browser https://doi.org/10.6084/m9.figshare.8175014 (Concatenate file a and b after decompressing)
2. Open a terminal in the directory and run `longread-UMI-pipeline -d <name-of-fastq-file> -s 1000000 -c 30 -t <Number-of-threads>`. 

### Generate data from paper
1. Create a working directory, open a terminal, download the Zymo mock fastq data and decompress:
   ENA: From terminal `wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR333/003/ERR3336963/ERR3336963.fastq.gz`
   Figshare: From browser https://doi.org/10.6084/m9.figshare.8175014 (Concatenate file a and b after decompressing)
2. Open a terminal in the directory and run `/path/to/longread-UMI-pipeline/longread_UMI_mockanalysis.sh <Number-of-threads>`. 
