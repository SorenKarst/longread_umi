# longread_umi 
A collection of scripts for processing longread UMI data.

<b> Pipeline version reference: </b> \
SM Karst, RM Ziels, RH Kirkegaard, EA Sørensen, D. McDonald, Q Zhu, R Knight, & M Albertsen. (2019). Enabling high-accuracy long-read amplicon sequences using unique molecular identifiers with Nanopore or PacBio sequencing. bioRxiv, 645903.
https://www.biorxiv.org/content/10.1101/645903v3

## Installation

### Conda

1. Requirements/Dependencies \
   OS tested (Linux 3.10.0, Ubuntu 14.04, Ubuntu 16.04) \
  `usearch` >=10
2. Download installer script from terminal \
   `wget https://raw.githubusercontent.com/SorenKarst/longread-UMI-pipeline/master/scripts/install_conda.sh`
3. Run installation script from terminal and follow instructions \
   `bash ./install_conda.sh` 

### Manual

1. Requirements/Dependencies \
   OS tested (Linux 3.10.0, Ubuntu 14.04, Ubuntu 16.04) \
   See `scripts/longread-UMI-pipeline_version_dump.txt`
2. Clone from github in terminal \
   `git clone https://github.com/SorenKarst/longread-UMI-pipeline.git`
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
1. Test the initialization command in terminal:  
  `longread_umi -h` or `/path/to/longread_umi.sh -h`
2. Test the nanopore_pipeline in terminal:  
  `longread_umi nanopore_pipeline -h` or `/path/to/longread_umi.sh nanopore_pipeline -h`
3. Test longread_umi nanopore_pipeline and qc_pipeline on test data:  
   Go to /path/to/longread-UMI-pipeline/test_data and open a terminal in the directory.
4. Run nanopore pipeline:  
   `longread_umi nanopore_pipeline \`  
     `-d test_reads.fq \`  
     `-o test \`  
     `-v 30 \`  
     `-w rrna_operon \`  
     `-t 1 \`  
     `-c 3 \`  
     `-p 1 \`  
     `-q r941_min_high_g330`  
   or  
   `longread_umi nanopore_pipeline \`  
     `-d test_reads.fq \`  
     `-v 30 \`  
     `-o test \`  
     `-s 90 \`  
     `-e 90 \`  
     `-m 3500 \`  
     `-M 6000 \`  
     `-f CAAGCAGAAGACGGCATACGAGAT \`  
     `-F AGRGTTYGATYMTGGCTCAG \`  
     `-r AATGATACGGCGACCACCGAGATC \`  
     `-R CGACATCGAGGTGCCAAAC \`  
     `-c 3 \`  
     `-p 1 \`  
     `-q r941_min_high_g330 \`  
     `-t 1`
5. Run qc pipeline:  
   `longread_umi qc_pipeline \`  
     `-d test_reads.fq \`  
     `-c consensus_raconx3_medakax1.fa \`  
     `-r zymo_curated \`  
     `-t 1`  

### Run pipeline on Zymo mock Nanopore data and perform qc
1. Create a working directory and open a terminal
2. Download the Zymo mock Nanopore R9.4.1 fastq data and decompress:  
   `wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR333/003/ERR3336963/ERR3336963_1.fastq.gz`  
   `gunzip -c ERR3336963_1.fastq.gz > reads.fq`
3. Download the SILVA database and decompress:  
   `wget https://www.arb-silva.de/fileadmin/silva_databases/release_132/Exports/SILVA_132_SSURef_Nr99_tax_silva.fasta.gz`  
   `gunzip -c SILVA_132_SSURef_Nr99_tax_silva.fasta.gz | sed '/^>/! s/U/T/g' >  silva_db.fasta` 
4. Open a terminal in the directory and run:  
   `longread_umi nanopore_pipeline \`  
     `-d reads.fq \`  
     `-o . \`  
     `-v 30  \`  
     `-w rrna_operon \`  
     `-t <Number-of-threads> \`  
     `-q r941_min_high_330`  
5. Open a terminal in the directory and run:  
   `longread_umi qc_pipeline \`  
     `-d "umi_binning/trim/reads_tf.fq;reads.fq" \`  
     `-c "consensus_racon_medaka_medaka.fa;variants.fa" \`  
     `-r "zymo_curated;zymo_vendor;variants.fa" \`  
     `-s silva_db.fasta \`  
     `-t <Number-of-threads>`  
