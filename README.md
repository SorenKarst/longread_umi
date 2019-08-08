# longread-UMI-pipeline
A collection of scripts for processing longread UMI data.

## Requirements
1. Tested on Ubuntu 14.04 (Yeah, we know...)
2. Dependencies: See scripts/dependencies.sh and scripts/longread_umi_version_dump.txt

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

### Create symlink to run longread-UMI-pipeline commands from terminal
1. Create symlink in ~/bin by opening a terminal and run:  
   `mkdir -p ~/bin`  
   `ln -s /path/to/longread-UMI-pipeline/longread_umi.sh ~/bin/longread_umi`  

### Change paths to dependencies
1. Open /path/to/longread-UMI-pipeline/scripts/dependencies.sh in a texteditor.
2. Change all paths under "Program paths" to reflect installation paths on your system.
3. If unsure of the paths try to type `which <function>` in the terminal i.e. `which racon`.
4. Install any missing dependencies.

### Customize porechop adaptors.py to be able to detect custom primers
1. We recommend to make a seperate installation of porechop to use with the longread-UMI-pipeline.
2. Go to path/to/porechop/porechop/
3. Backup current adapters.py.
4. Replace current adapters.py with path/to/longread-UMI-pipeline/scripts/adapters.py.

### Test scripts
1. Open a terminal anywhere and test the initialization command:  
  `longread_umi -h` or `/path/to/longread_umi.sh -h`
2. Open a terminal anywhere and test the nanopore_pipeline command:  
  `longread_umi nanopore_pipeline -h` or `/path/to/longread_umi.sh nanopore_pipeline -h`
3. Test longread_umi nanopore_pipeline and qc_pipeline on test data:  
   Go to /path/to/longread-UMI-pipeline/test_data and open a terminal in the directory.
4. Run `longread_umi nanopore_pipeline -d test_reads.fq -s 10 -c 30 -t 1`
5. Run `longread_umi qc_pipeline -d test_reads.fq -c consensus_racon_medaka_medaka.fa -r zymo_curated -t 1`

### Run pipeline on Zymo mock Nanopore data and perform qc
1. Create a working directory, open a terminal, download the Zymo mock Nanopore fastq data and decompress:  
   `wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR333/003/ERR3336963/ERR3336963_1.fastq.gz; gunzip -c ERR3336963_1.fastq.gz > reads.fq`  
2. Open a terminal in the directory and run:  
  `longread_umi nanopore_pipeline -d reads.fq -s 1000000 -c 30 -t <Number-of-threads>`
3. Open a terminal in the directory and run:  
  `longread_umi qc_pipeline`  
     `-d "umi_binning/trim/reads_tf.fq;reads.fq"`  
     `-c "consensus_racon_medaka_medaka.fa;variants_all.fa"`  
     `-r "zymo_curated;zymo_vendor;variants_all.fa"`  
     ` -t <Number-of-threads>`  
