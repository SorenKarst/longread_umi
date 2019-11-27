# longread-UMI-pipeline
A collection of scripts for processing longread UMI data.
Tested on Linux 3.10.0

<b> Citation: </b> \
Karst, S. M., Ziels, R. M., Kirkegaard, R. H., & Albertsen, M. (2019). Enabling high-accuracy long-read amplicon sequences using unique molecular identifiers and Nanopore sequencing. bioRxiv, 645903.
https://www.biorxiv.org/content/10.1101/645903v2

## Installing the Easy Way -- via Conda
### Requirements/Dependencies 
Conda or Miniconda3 installed  
`usearch` version 10 or higher

### Clone from github
1. Go to desired installation directory, open a terminal and run:  \
   `git clone https://github.com/ziels/longread-UMI-pipeline`

2. Go to scripts directory: \
   `cd longread-UMI-pipeline/scripts`  

3. Modify `dependencies.sh` with path to `usearch`  
Change line `export USEARCH=usearch_path` to give your exact file path to the `usearch` executable file (instead of `usearch_path`). 

### Install conda environment 
   `conda env create -f environment.yaml`

1. Check that Conda env is installed \
   `conda info —-envs` 
   
Make sure you get something like: 

   > `# conda environments:`\
   > `#`\
   >  ` longread-UMI       <path to conda envs>/longread-UMI ` 
   
Note the installation path of the longread-UMI environment (for next steps)

2. Activate conda environment \
   `conda activate longread-UMI` \
   Or, depending on your conda version: `source activate longread-UMI`

### Editing porechop adapters with custom primers
1. Find path of conda environments from command \
   `conda info —-envs` 

2. Check the `porechop` path works: \
   `ls < path to conda environments >/longread-UMI/lib/python3.6/site-packages/porechop` 

Make sure you see an `adapters.py` returned from the above command.

3. Back-up and replace `adapters.py`\
   `mv <path_to_conda_environments>/longread-UMI/lib/python3.6/site-packages/porechop/adapters.py <path_to_conda_environments>/longread-UMI/lib/python3.6/site-packages/porechop/adapters_old.py`

(From within `longread-UMI-pipeline/scripts` directory):\
   `cp ./adapters.py <path_to_conda_environments>/longread-UMI/lib/python3.6/site-packages/porechop/adapters.py`

## Test longread-UMI-pipeline on test data:  
   Go to `/path/to/longread-UMI-pipeline/test_data`   
   Open a terminal in the directory and run \
   `../longread-UMI-pipeline.sh -d test_reads.fq -s 10 -c 30 -t 1`


## Installing the Hard Way -- Manual installation
### Requirements/Dependencies


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

## Run pipeline on Zymo mock data
1. Create a working directory, open a terminal, download the Zymo mock fastq data and decompress:  
   `wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR333/003/ERR3336963/ERR3336963_1.fastq.gz; gunzip -c ERR3336963_1.fastq.gz > reads.fq`  
2. Open a terminal in the directory and run:  
  `longread-UMI-pipeline -d reads.fq -s 1000000 -c 30 -t <Number-of-threads>`
3. Open a terminal in the directory and run:  
  `longread-UMI-mockanalysis <Number-of-threads>`
