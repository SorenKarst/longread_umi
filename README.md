# longread_umi 

A collection of scripts for processing longread UMI data.

**Table of contents**
- [Installation](#installation)
- [Quick start](#quick-start)
- [Data and examples](#data-and-examples)
- [Usage](#usage)

**Citation**  
SM Karst, RM Ziels, RH Kirkegaard, EA Sørensen, D. McDonald, Q Zhu, R Knight, & M Albertsen. (2020). Enabling high-accuracy long-read amplicon sequences using unique molecular identifiers with Nanopore or PacBio sequencing. [bioRxiv, 6459039](https://www.biorxiv.org/content/10.1101/645903v3).

## Installation

### Conda

1. Requirements/Dependencies  
   OS tested (Linux 3.10.0, Ubuntu 14.04, Ubuntu 16.04)  
    `usearch` >=10
2. Download installer script from terminal  
   ```
   wget https://raw.githubusercontent.com/SorenKarst/longread_umi/master/scripts/install_conda.sh
   ```
3. Run installation script from terminal and follow instructions (< 10 min on desktop)  
   `bash ./install_conda.sh` 
4. Initiate conda and refresh terminal before using pipeline.  
   `conda init; source ~/.bashrc`  
5. Activate and deactivate conda environment
   ```
   conda activate longread_umi
   ...
   conda deactivate
   ```

### Manual

1. Requirements/Dependencies  
   OS tested (Linux 3.10.0, Ubuntu 14.04, Ubuntu 16.04)  
   See `scripts/longread-UMI-pipeline_version_dump.txt`
2. Clone from github in terminal  
   `git clone https://github.com/SorenKarst/longread_umi.git`
3. Make bash scripts executable  
   `find ./longread_umi -name "*.sh" -exec chmod +x {} \;`
4. Install dependencies  
   See `./longread_umi/scripts/install_dependencies.sh` for inspiration.
5. Change paths to dependencies  
   Modify `./longread_umi/scripts/dependencies.sh` in a texteditor.
6. Customize porechop adaptors.py to be able to detect custom primers  
   Replace current `adapters.py` with `./longread_umi/scripts/adapters.py`

## Quick start

### Test data
1. Test the initialization command in terminal  
    `longread_umi -h` or `/path/to/longread_umi.sh -h`
2. Test the nanopore_pipeline in terminal  
    `longread_umi nanopore_pipeline -h` or `/path/to/longread_umi.sh nanopore_pipeline -h`
3. Test longread_umi nanopore_pipeline and qc_pipeline on test data:  
   Go to /path/to/longread-UMI-pipeline/test_data and open a terminal in the directory.
4. Run nanopore pipeline (< 10 minutes on desktop)
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
   
   Expected output
   - `consensus_raconx3_medakax1.fa` containing 8 UMI consensus sequences
   - `variants.fa` containing 2 variant consensus sequences

7. Run qc pipeline (< 5 minutes on desktop)
   ```
   longread_umi qc_pipeline \
     -d test_reads.fq \
     -c test/consensus_raconx3_medakax1.fa \
     -r zymo_curated \
     -t 1 \
     -u test \
     -o test/qc
   ```
   Expected output
   
   - ...


### Zymomock rRNA operon data
1. Download the Zymomock rRNA operon Nanopore R9.4.1 fastq data and decompress
   ```
   wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR333/003/ERR3336963/ERR3336963_1.fastq.gz 
   gunzip -c ERR3336963_1.fastq.gz > reads.fq
   ```
2. Run nanopore pipeline (~ 500 CPU hours)
   ```
   longread_umi nanopore_pipeline \
     -d reads.fq \
     -o umi_out \
     -v 30 \
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
     -t <Number-of-threads>
   ```
5. Run qc pipeline
   ```
   longread_umi qc_pipeline \
     -d "umi_out/umi_binning/trim/reads_tf.fq;reads.fq" \
     -c "umi_out/consensus_raconx3_medakax1.fa;umi_out/variants.fa" \
     -r "zymo_curated" \
     -u umi_out \
     -o umi_out/qc \
     -t <Number-of-threads> 
   ```

## Data and examples

- [Example data](docs/DATA.md)
  
- [ONT R10 Zymomock rRNA - generate UMI consensus sequences and validate data](https://htmlpreview.github.io/?https://github.com/SorenKarst/longread_umi/blob/master/docs/ONT_R10_ZYMO_rRNA.html)
  
  
  

## Usage

```

-- longread_umi: pipelines and tools for longread UMI processing.

usage: longread_umi [-h] [ name ...]

where:
    -h   Show this help text.
    name Name of tool or pipeline.
    ...  Commands for tool or pipeline.

Pipelines:

   nanopore_pipeline   Generate UMI consensus sequences from Nanopore data
   pacbio_pipeline     Generates UMI consensus sequences from PacBio CCS data.
   qc_pipeline         UMI consensus data statistics and compare to references

Tools:

   consensus_racon          Generate UMI consensus sequence with racon
   demultiplex              Dual barcode demultiplexing
   demultiplex_3end         3'-end dual barcode demultiplexing
   nanopore_settings_test   Test impact of polishing rounds on UMI consensus.
   polish_medaka            Nanopore UMI consensus polishing with Medaka
   primer_position          Locate adapter and primer positions in read data
   trim_amplicon            Trimming sequences based on primers
   umi_binning              Longread UMI detection and read binning.
   variants                 Phase and call variants from UMI consensus sequences.

For help with a specific tool or pipeline:
longread_umi <name> -h
```
  
  

```

-- longread_umi consensus_racon: Generate UMI consensus sequence with racon
   Raw read centroid found with usearch and used as seed for
   (r) x times racon polishing.

usage: consensus_racon [-h] (-d dir -o dir -r value -t value -n file) 

where:
    -h  Show this help text.
    -d  Directory containing UMI read bins in the format
        'umi*bins.fastq'. Recursive search.
    -o  Output directory.
    -r  Number of racon polishing rounds.
    -t  Number of threads to use.
    -n  Process n number of bins. If not defined all bins
        are processed.
```
  
  

```

-- longread_umi demultiplex: Dual barcode demultiplexing

   Script for demultiplexing UMI consensus sequences based on 
   custom barcodes. The script demultiplexes raw read data
   and assigns the consensus sequences to a sample by majority vote
   of the raw read assignments. Post processing demultiplxing optimizes 
   consensus yield. The script expects dual barcodes in a barcode file.
   If the same barcode is used in both ends simply repeat barcode.

usage: demultiplex [-h] (-c file -r file -u file -o dir -b file)
(-p string -n range -t value) 

where:
    -h  Show this help text.
    -c  UMI consensus sequences that need demultiplexing.
    -r  Raw read sequences that were used to generate
        the consensus sequences.
    -u  List of raw read names and UMI bin assignments.
    -o  Output directory.
    -b  File containing barcodes. 
        [Default = longread_umi/scripts/barcodes.tsv].
    -p  Barcode name prefix [Default = 'barcode'].
    -n  Barcode numbers used. [Default  = '1-120'].
    -t  Number of threads used.
```
  
  

```

-- longread_umi demultiplex_3end: 3'-end dual barcode demultiplexing

    Script for demultiplexing UMI consensus sequences based on 
    custom barcodes. The script demultiplexes raw read data
    and assigns the consensus sequences to a sample by majority vote
    of the raw read assignments. Post processing demultiplxing optimizes 
    consensus yield. The script expects dual barcodes in a barcode file.
    If the same barcode is used in both ends simply repeat barcode. This
    version of the script only looks for barcodes in the 3' end. This demultiplexing
    is for data types which are truncated in the 5' end. Dual barcoding is
    still used.
	
usage: demultiplex_3end [-h] (-c file -r file -u file -o dir -b file -p string)
(-n range -m value -t value) 

where:
    -h  Show this help text.
    -c  UMI consensus sequences that need demultiplexing.
    -r  Raw read sequences that were used to generate
        the consensus sequences.
    -u  List of raw read names and UMI bin assignments.
    -o  Output directory.
    -b  File containing barcodes.
        Default is longread_umi/scripts/barcodes.tsv
    -p  Barcode name prefix. [Default = 'barcode'].
    -n  Barcode numbers used. [Default = '1-120'].
    -m  Minimum number of barcodes found to demultiplex
        sequences. Default 2.
    -t  Number of threads used.
```
  
  

```

-- longread_umi nanopore_pipeline: Generate UMI consensus sequences from Nanopore data
   
usage: nanopore_pipeline [-h] [-w string] (-d file -v value -o dir -s value) 
(-e value -m value -M value -f string -F string -r string -R string )
( -c value -p value -n value -u dir -t value -T value ) 

where:
    -h  Show this help text.
    -d  Single file containing raw Nanopore data in fastq format.
    -v  Minimum read coverage for using UMI consensus sequences for 
        variant calling.
    -o  Output directory.
    -s  Check start of read up to s bp for UMIs.
    -e  Check end of read up to f bp for UMIs.
    -m  Minimum read length.
    -M  Maximum read length.
    -f  Forward adaptor sequence. 
    -F  Forward primer sequence.
    -r  Reverse adaptor sequence.
    -R  Reverse primer sequence.
    -c  Number of iterative rounds of consensus calling with Racon.
    -p  Number of iterative rounds of consensus calling with Medaka.
    -q  Medaka model used for polishing. r941_min_high, r10_min_high etc.
    -w  Use predefined workflow with settings for s, e, m, M, f, F, r, R.
        rrna_operon [70, 80, 3500, 6000, CAAGCAGAAGACGGCATACGAGAT,
        AGRGTTYGATYMTGGCTCAG, AATGATACGGCGACCACCGAGATC, CGACATCGAGGTGCCAAAC]
        Overwrites other input.
    -n  Process n number of bins. If not defined all bins are processed.
        Pratical for testing large datasets.
    -u  Directory with UMI binned reads.
    -t  Number of threads to use.
    -T  Number of medaka jobs to start. Threads pr. job is threads/jobs.
        [Default = 1].
```
  
  

```

-- longread_umi nanopore_settings_test: Test impact of polishing rounds on UMI consensus.

usage: nanopore_settings_test [-h] (-d file -n value -c value -o dir -s value -e value) 
(-m value -M value -f string -F string -r string -R string -t value -T value) 
(-x value -y value -q string -p -u dir ) 

where:
    -h  Show this help text.
    -d  Single file containing raw Nanopore data in fastq format.
    -o  Output directory.
    -s  Check start of read up to s bp for UMIs.
    -e  Check end of read up to f bp for UMIs.
    -m  Minimum read length.
    -M  Maximum read length.
    -f  Forward adaptor sequence. 
    -F  Forward primer sequence.
    -r  Reverse adaptor sequence.
    -R  Reverse primer sequence.
    -n  Process n number of bins. If not defined all bins are processed.
        Pratical for testing large datasets.
    -w  Use predefined workflow with settings for s, e, m, M, f, F, r, R.
        rrna_operon [70, 80, 3500, 6000, CAAGCAGAAGACGGCATACGAGAT,
        AGRGTTYGATYMTGGCTCAG, AATGATACGGCGACCACCGAGATC, CGACATCGAGGTGCCAAAC]
    -t  Number of threads to use.
    -T  Number of medaka jobs to start. Threads pr. job is threads/jobs.
        [Default = 1].
    -x  Test Racon consensus rounds from 1 to <value>.
    -y  Test Medaka polishing rounds from 1 to <value>.
    -q  Medaka model used for polishing. r941_min_high, r10_min_high etc.
    -p  Flag to disable Nanopore trimming and filtering.
    -u  Directory with UMI binned reads.

Test run:
longread_umi nanopore_settings_test 
  -d test_reads.fq 
  -o settings_test 
  -w rrna_operon 
  -t 100 
  -T 20 
  -x 4 
  -y 3 
  -n 1000
```
  
  

```

-- longread_umi pacbio_pipeline: Generates UMI consensus sequences from PacBio CCS data.
   
usage: pacbio_pipeline [-h] (-d file -v value -o dir -s value -e value) 
(-m value -M value -f string -F string -r string -R string -c value -w string)
(-n value -u dir -t value) 

where:
    -h  Show this help text.
    -d  Single file containing PacBio CCS read data in fastq format.
    -v  Minimum read coverage for using UMI consensus sequences for 
        variant calling.
    -o  Output directory.
    -s  Check start of read up to s bp for UMIs.
    -e  Check end of read up to f bp for UMIs.
    -m  Minimum read length.
    -M  Maximum read length.
    -f  Forward adaptor sequence. 
    -F  Forward primer sequence.
    -r  Reverse adaptor sequence.
    -R  Reverse primer sequence.
    -c  Number of iterative rounds of consensus calling with Racon.
    -w  Use predefined workflow with settings for s, e, m, M, f, F, r, R, c.
        rrna_operon [70, 80, 3500, 6000, CAAGCAGAAGACGGCATACGAGAT,
        AGRGTTYGATYMTGGCTCAG, AATGATACGGCGACCACCGAGATC, CGACATCGAGGTGCCAAAC, 2]
    -n  Process n number of bins. If not defined all bins are processed.
        Pratical for testing large datasets.
    -u  Directory with UMI binned reads.
    -t  Number of threads to use.
```
  
  

```

-- longread_umi polish_medaka: Nanopore UMI consensus polishing with Medaka
   
usage: polish_medaka [-h] [-l value T value] 
(-c file -m string -d dir -o dir -t value -n file -T value)

where:
    -h  Show this help text.
    -c  File containing consensus sequences.
    -m  Medaka model.
    -l  Expected minimum chunk size. [Default = 6000]
    -d  Directory containing UMI read bins in the format
        'umi*bins.fastq'. Recursive search.
    -o  Output directory.
    -t  Number of threads to use.
    -n  Process n number of bins. If not defined all bins
        are processed.
    -t  Number of Medaka jobs to run. [Default = 1].
```
  
  

```

-- longread_umi primer_position: Locate adapter and primer positions in read data
    Script for checking position of adapters and gene specific primers flanking
    UMI sequences in read terminals. Relevant if using custom UMI adapters/primers,
    sample barcoding or if basecalling/processing truncates reads.
   
usage: primer_position [-h] [-e value -n value ] (-d value -o dir -t value)
(-f string -F string -r string -R string ) 

where:
    -h  Show this help text.
    -d  Raw fastq reads.
    -o  Output directory
    -t  Number of threads to use.
    -f  Forward adaptor sequence. 
    -F  Forward primer sequence.
    -r  Reverse adaptor sequence.
    -R  Reverse primer sequence.
    -e  Length of terminal end to search for primers. [Default = 500]
    -n  Subset reads before search. [Default = 100000]
```
  
  

```

-- longread_umi qc_pipeline: UMI consensus data statistics and compare to references
   Calculates data statistics (UMI bins size, UMI cluster size, yield, read length etc.).
   Mapping of read data and UMI consensus sequences to reference sequences to allow for 
   error profiling. Detects chimeras using uchime2_ref. Detects contamination by
   comparing mapping results to known references and the SILVA database - only works
   if reference database contains all expected sequences. Alternatively, use variants.fa
   as reference database.
   
usage: qc_pipeline [-h] (-d files -c files -r files -s file -u dir -o dir -t value) 

where:
    -h  Show this help text.
    -d  List of read files seperated by ';'
        i.e. 'reads.fq;trim/reads_tf.fq'
        First read file used for read classification. 'reads_tf.fq' recommended.
    -c  List of consensus files seperated by ';'
        i.e. 'consensus_medaka_medaka.fa;racon/consensus_racon.fa'
        First consensus file used for comparison to alternative refs
        and for chimera checking. Subsequent consensus sequences only mapped to
        first reference.
    -r  List of reference files seperated by ';'.
        First reference is used for all mappings. Subsequent references
        only used for mapping first consensus file.
        'zymo_curated' refers to:
        longread_umi/scripts/zymo-ref-uniq_2019-10-28.fa
        'zymo_vendor' refers to:
        longread_umi/scripts/zymo-ref-uniq_vendor.fa
    -s  SILVA reference database in fasta format used for detecting contamination.
    -u  UMI consensus output folder.
    -o  Output folder. Default 'qc'.
    -t  Number of threads to use.

Example of SILVA database download:
wget https://www.arb-silva.de/fileadmin/silva_databases/
release_132/Exports/SILVA_132_SSURef_Nr99_tax_silva.fasta.gz
gunzip SILVA_132_SSURef_Nr99_tax_silva.fasta.gz
```
  
  

```

-- longread_umi trim_amplicon: Trimming sequences based on primers
   
usage: trim_amplicon [-h] (-d dir(s) -p string(s) -o dir )
(-F string -R string -m value -M value -t value -l dir) 

where:
    -h  Show this help text.
    -d  Directory to look for sequence files using pattern (-p).
        Mutliple directories can be seperated by ';'.
    -p  File pattern(s) to look for. Multiple patterns
        can be separared by ';' and flanked by '"..."'.
    -o  Output directory.
    -F  Forward primer sequence.
    -R  Reverse primer sequence.
    -m  Minimum read length.
    -M  Maximum read length.
    -t  Number of threads to use.
    -l  Log directory
```
  
  

```

-- longread_umi umi_binning: Longread UMI detection and read binning.
   Tool requires UMIs in both ends of the read flanked by defined
   adaptor regions.

usage: umi_binning [-h] (-d file -o dir -m value -M value )
(-s value -e value -f string -F string -r string -R string -p -t value) 

where:
    -h  Show this help text.
    -d  Reads in fastq format.
    -o  Output directory.
    -m  Minimum read length.
    -M  Maximum read length.
    -s  Check start of read up to s bp for UMIs.
    -e  Check end of read up to f bp for UMIs.
    -f  Forward adaptor sequence. 
    -F  Forward primer sequence.
    -r  Reverse adaptor sequence.
    -R  Reverse primer sequence.
    -p  Flag to disable Nanopore trimming and filtering. Use with PacBio reads.
    -t  Number of threads to use.
```
  
  

```

-- longread_umi variants: Phase and call variants from UMI consensus sequences.
   This is a naive variant caller, which phases UMI consensus sequences
   based on SNPs and calls a variant with >=3x coverage. Reads are initially 
   grouped by read clustering at 99.5% identity and a centroid sequence is picked.
   The centroid sequence is used as a mapping reference for all reads in the cluster
   to detect SNPs for phasing and variant calling. Before read clustering homopolymers
   are masked and then reintroduced before variant calling.
   
usage: variants [-h -b] (-c file -o dir -t value ) 

where:
    -h  Show this help text.
    -c  UMI consensus file.
    -o  Output directory.
    -t  Number of threads to use. [Default = 1]
    -b  Debug flag. Keep temp files. [Default = NO]
```
  
  


## License
[GNU General Public License, version 3](LICENSE)
