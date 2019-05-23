# longread-UMI-pipeline
A collection of scripts for processing longread UMI data.

## Installation

### Clone from github
1. Go to desired installation directory, open a terminal and run commands:
2. git clone https://github.com/SorenKarst/longread-UMI-pipeline.git

### Make bash scripts executable
1. Go to longread-UMI-pipeline directory, open a terminal and run commands:
2. find . -name "*.sh" -exec chmod +x {} \;

### Create symlink to run longread-UMI-pipeline from terminal
1. Open terminal and run command:
   ln -s /path/to/longread-UMI-pipeline/longread_UMI_pipeline.sh ~/bin/longread-UMI-pipeline
2. Test symlink by opening a terminal window anywhere
3. Type longread-UMI-pipeline -h

### Change paths to dependencies
1. Go to /path/to/longread-UMI-pipeline/scripts and open dependencies.sh in a texteditor.
2. Change all paths under "Paths to dependencies" to reflect installations on you system.
3. If unsure of the paths try to type `which <function>` in the terminal. I.e. `which racon`.

## Running longread-UMI-pipeline

### Running longread-UMI-pipeline on test data
1. Go to /path/to/longread-UMI-pipeline/test_data
2. Open a terminal in the directory and run `longread-UMI-pipeline -d test_reads.fq -s 10 -c 30 -t 1`
   Or run ./test_longread_umi_pipeline.sh.

### Run pipeline on Zymo mock data
1. Download Zymo mock fastq data from ENA and place it in a working directory.
2. Open a terminal in the directory and run `longread-UMI-pipeline -d <name-of-fastq-file> -s 1000000 -c 30 -t <Number-of-threads>`. 

### Generate data from paper
1. Download Zymo mock fastq data from ENA and place it in a working directory.
2. Open a terminal in the directory and run `/path/to/longread-UMI-pipeline/longread_UMI_mockanalysis.sh <Number-of-threads>`. 
