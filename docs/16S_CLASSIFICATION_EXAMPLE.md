## Exploration of different 16S amplicon data types for microbial community profiling

Søren Karst, 2020-02-18

 ### Introduction

The aim of this hands-on session is to compare three different types of amplicon data (Illumina, Nanopore and Nanopore UMI) for microbial community profiling using the 16S marker gene. We will use data from a simple mock community of 8 bacterial species to investigate resolution and taxonomic classification and discuss the impact on the observed microbial profiles. Afterwards we will look at the negative impact of missing closely related reference sequences in the databases used for taxonomic classification, and how generating our own reference sequences mitigate some of these problems (and create new ones). 



#### Getting started

* Open a terminal and create a working directory

  ```
  mkdir amplicon_umi
  cd amplicon_umi
  ```

* Create conda environment containing the required software for the data analysis. Remember the environment has to be activated before you can perform the analysis.

  ```
  conda create -n amplicon_umi --channel bioconda cutadapt vsearch seqtk gawk 
  
  conda activate amplicon_umi
  ...
  conda deactivate
  ```

  

* Generate example data

  All the data actually originates from Nanopore UMI full rRNA operon data set. Raw 16S data set is simply the raw Nanopore data before UMI processing where the 16S part has been extracted, the 16S UMI data is the result of UMI processing of the raw data with extraction of the 16S, and lastly the Illumina data is generated in silico by extracting the V4 region from the 16S UMI data.

  ```
  # Download data
  wget "ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR381/ERR3813594/smkj412_zymo_D6306_rrna_umi_ont_min_r10_g344hac.fq.gz"
  wget "ftp://ftp.sra.ebi.ac.uk/vol1/ERZ128/ERZ1284843/smkj412_zymo_D6306_rrna_umi_ont_r10_g344hac_con.fasta.gz"
  
  # Decompress data
  gunzip \
    -c smkj412_zymo_D6306_rrna_umi_ont_min_r10_g344hac.fq.gz \
    > ont_rrna_raw.fastq
  gunzip \
    -c smkj412_zymo_D6306_rrna_umi_ont_r10_g344hac_con.fasta.gz \
    > ont_rrna_umi.fasta
  rm -rf *.gz
    
  # Generate 16S data ------------------------------------------------------------
  
  # Extract 16S from raw
  cutadapt \
    --discard-untrimmed \
    -j 100 \
    -m 1200 \
    -M 2000 \
    -e 0.2 \
    -g AGRGTTYGATYMTGGCTCAG...TGYACWCACCGCCCGTC \
    -g GACGGGCGGTGWGTRCA...GACGGGCGGTGWGTRCA \
    ont_rrna_raw.fastq > ont_16s_raw_all.fastq  
  gawk '
    NR%4==1 && NR<=100000{
      print(">" substr($1,2) ";sample=ont_16s_raw;")
      getline; print
    }
  ' ont_16s_raw_all.fastq \
  > ont_16s_raw.fasta
  
  # Extract 16S from UMI  
  cutadapt \
    --discard-untrimmed \
    -j 100 \
    -m 1200 \
    -M 2000 \
    -e 0.2 \
    -a TGYACWCACCGCCCGTC \
    ont_rrna_umi.fasta > ont_16s_umi_all.fasta
  gawk '
    NR%2==1 && NR<=50000{
      print($1 ";sample=ont_16s_umi;")
      getline; print
    }
  ' ont_16s_umi_all.fasta \
  > ont_16s_umi.fasta
  
  # Extract V4 from UMI  
  cutadapt \
    --discard-untrimmed \
    -j 100 \
    -m 200 \
    -M 300 \
    -e 0.2 \
    -g CAGCMGCCGCGGTAA...GGATTAGATACCCBNGTA \
    ont_rrna_umi.fasta > ont_v4_umi_all.fasta
  gawk '
    NR%2==1 && NR<=50000{
      print($1 ";sample=ilm_v4;")
      getline; print
    }
  ' ont_v4_umi_all.fasta \
  > ilm_v4.fasta
  ```

* Fetch a minimal version of the SILVA database containing 16S reference sequences and taxonomy. This is a minimal version of the full database.

  ```
  wget https://www.drive5.com/sintax/ltp_16s_v123.fa.gz
  gunzip -c ltp_16s_v123.fa.gz > silva.fasta
  ```

  


#### About the example amplicon data

* The microbial community

  It is the ZymoBIOMICS Microbial Community DNA Standard (D6306, lot no. ZRC190811), which consists of 8 bacterial species and 2 fungal species. The amplicon primers used to generate the sequencing libraries only target bacteria and hence we should only see these in the data: *Bacillus subtilis*, *Enterococcus faecalis*, *Escherichia coli*, *Lactobacillus fermentum*, *Listeria monocytogenes*, *Pseudomonas aeruginosa*, *Salmonella enterica* and *Staphylococcus aureu*s. In total there are 43 unique rRNA operons between the 8 bacteria, and depending on the resolution of the different data types these will be aggregated into fewer sequences e.g. for 16S there should be 27 unique sequences.

* The sequence data

  Illumina short-read 16S V4 amplicons (~250 bp, < 0.1% error rate): ilm_v4.fasta

  Nanopore raw full-length 16S amplicons (~1350 bp, ~ 10% error rate): ont_16s_raw.fasta

  Nanopore UMI full-length 16S amplicons (~ 1350 bp, ~ < 0.01% error rate): ont_16s_umi.fasta

  All data types have been subset to 25,000 reads.



#### Generating count tables and performing classification with a good reference database

*Illumina short-read 16S V4 amplicons*

1. Dereplicate the data. This step reports all unique sequences and their count.

```
vsearch \
   --derep_fulllength ilm_v4.fasta \
   --minuniquesize 2 \
   --sizeout \
   --fasta_width 0 \
   --output ilm_v4_derep.fasta
```

2. Denoise dereplicated data. This step removes erroneous sequences that are derived from more abundant true sequences. The output is `amplicon sequence variants` or ASVs. The ASV sequences represent all sequence types believe to be found in our dataset.

```
vsearch \
  --cluster_unoise ilm_v4_derep.fasta \
  --sizein \
  --sizeout \
  --relabel asv \
  --centroids ilm_v4_asv.fasta
```

3. Remove chimeras.

```
vsearch \
   --uchime3_denovo ilm_v4_asv.fasta \
   --sizein \
   --sizeout \
   --fasta_width 0 \
   --nonchimeras ilm_v4_asvf.fasta
```

4. Taxonomic classification of ASVs. 

```
vsearch \
  --sintax ilm_v4_asvf.fasta \
  --db silva.fasta \
  --sintax_cutoff 0.80 \
  --tabbedout ilm_v4_asv.tax
```

5. `gawk` magic to add taxonomy to ASV headers.

```
gawk '
  # Load taxonomic classifications
  FILENAME ~ /_asv.tax/{
    TAX[$1]=$4
  }
  # Load ASV sequences and add taxonomy to headers
  FILENAME ~ /_asvf.fasta/{
    if ($0 ~ /^>/){
      print $0 ";tax=" TAX[substr($0, 2)] ";"
    }
    if ($0 !~ /^>/){
      print $0
    }
  } ' \
ilm_v4_asv.tax \
ilm_v4_asvf.fasta \
> ilm_v4_asv_tax.fasta
```

6. Generate count table with taxonomic classifications

```
vsearch \
  -usearch_global ilm_v4.fasta \
  --db ilm_v4_asv_tax.fasta \
  --id 0.97 \
  --otutabout ilm_v4_countable.tsv
```

7. Open `ilm_v4_countable.tsv` with a spreadsheet editor and take a look at the community profile. The profile is listed as three columns here. 1 = `#OTU id` is the ID of the ASV sequences found in the data. 2 =  `ilm_v4` is the name of the sample and this column contains the counts of the respective ASVs in this particular sample. If we had analyzed more samples, there would be more columns like this. 3 = `taxonomy` the taxonomic assignment of the ASV.  

* Did we find the bacteria we expect to find? 
* How does the number of ASVs compare to what we know about the number of bacterial operons? What does that say about resolution?
* Did we find any chimeras? If not try to create some chimeras yourself and add to the dataset and rerun the analysis.



*Nanopore raw full-length 16S amplicons*

The error rate in raw Nanopore data is too high to use conventional marker gene pipelines to generate ASVs or even OTUs. So we have to rely on our database reference sequences for making the count table for this type of data. There are different strategies for doing this, but here we will simply classify each sequence and create a count table by counting how many times we see a specific classification. This would be similar to using Oxford Nanopore WIMP workflow.

1. Classify raw Nanopore 16S sequences.

```
vsearch \
  --sintax ont_16s_raw.fasta \
  --db silva.fasta \
  --sintax_cutoff 0.80 \
  --tabbedout ont_16s_raw.tax
```

2. `gawk` magic to create a count table from the classified raw reads

```
gawk '
  {TAX[$4";"]++}
  END{
    print "ont_16s_raw\ttaxonomy"
    for (i in TAX){
      print i"\t"TAX[i]
    }
  }
' ont_16s_raw.tax \
> ont_16s_raw_countable.tsv
```

3. Open `ont_16s_raw_countable.tsv` with a spreadsheet editor and take a look at the community profile. 

   * Did we find the bacteria we expect to find? Did we find more?

   * What are the risks of defining and counting different bacteria by taxonomic classification alone?

   * What type of questions is this type of results useful for?

   * If you got time try to optimize the vsearch settings to get a better result.
   
     

*Nanopore UMI full-length 16S amplicons*

Amplicons generated by using UMIs have a very low error rate, even better than the base Illumina error rate. Hence, we can use the standard marker gene workflows for generating our count tables. 

1. Dereplicate the data.

```
vsearch \
   --derep_fulllength ont_16s_umi.fasta \
   --minuniquesize 2 \
   --sizeout \
   --fasta_width 0 \
   --output ont_16s_umi_derep.fasta
```

2. Denoise dereplicated data. 

```
vsearch \
  --cluster_unoise ont_16s_umi_derep.fasta \
  --sizein \
  --sizeout \
  --relabel asv \
  --centroids ont_16s_umi_asv.fasta
```

3. Remove chimeras.

```
vsearch \
   --uchime3_denovo ont_16s_umi_asv.fasta \
   --sizein \
   --sizeout \
   --fasta_width 0 \
   --nonchimeras ont_16s_umi_asvf.fasta
```

4. Taxonomic classification of ASVs. 

```
vsearch \
  --sintax ont_16s_umi_asvf.fasta \
  --db silva.fasta \
  --sintax_cutoff 0.80 \
  --tabbedout ont_16s_umi_asv.tax
```

5. `gawk` magic to add taxonomy to ASV headers.

```
gawk '
  # Load taxonomic classifications
  FILENAME ~ /_asv.tax/{
    TAX[$1]=$4
  }
  # Load ASV sequences and add taxonomy to headers
  FILENAME ~ /_asvf.fasta/{
    if ($0 ~ /^>/){
      print $0 ";tax=" TAX[substr($0, 2)] ";"
    }
    if ($0 !~ /^>/){
      print $0
    }
  } ' \
ont_16s_umi_asv.tax \
ont_16s_umi_asvf.fasta \
> ont_16s_umi_tax.fasta
```

6. Generate count table with taxonomic classifications

```
vsearch \
  -usearch_global ont_16s_umi.fasta \
  --db ont_16s_umi_tax.fasta \
  --id 0.97 \
  --otutabout ont_16s_umi_countable.tsv
```

7. Open `ont_16s_umi_countable.tsv` with a spreadsheet editor and take a look at the community profile. 
   * Did we find the bacteria we expect to find?
   * Did the classification become better? Compared to ONT raw 16S and ILM V4?
   * Is the resolution too high now?
   * What impact does the full 16S have on the observed relative abundance?
   * What type of questions is this type of results useful for?

#### Classification with poor reference databases

A database can be poor in many ways, but here we define poor as missing closely related reference sequences at genus level. We can simulate such an example our selves.

1. Sabotage the database by removing all reference sequences from the same genus as our mock community species are from.

```
seqtk \
  seq -l0 silva.fasta |\
  gawk \
    -v REMOVE_SEQ="$REMOVE_SEQ" '
    /^>/ && !/Bacillus|Enterococcus|Escherichia|Lactobacillus|Listeria|Pseudomonas|Salmonella|Staphylococcus/{
      print
      getline; print
    }
  ' > silva_poor.fasta
```

2. Perform classification for the three different datasets

```
# Illumina V4
vsearch \
  --sintax ilm_v4_asvf.fasta \
  --db silva_poor.fasta \
  --sintax_cutoff 0.80 \
  --tabbedout ilm_v4_asv_poor_tax.tsv
  
# Nanopore raw 16S
vsearch \
  --sintax ont_16s_raw.fasta \
  --db silva_poor.fasta \
  --sintax_cutoff 0.80 \
  --tabbedout ont_16s_raw_poor_tax.tsv
  
gawk '
  {TAX[$4";"]++}
  END{
    print "ont_16s_raw\ttaxonomy"
    for (i in TAX){
      print i"\t"TAX[i]
    }
  }
' ont_16s_raw_poor_tax.tsv \
> ont_16s_raw_poor_tax_countable.tsv

# Nanopore UMI 16S
vsearch \
  --sintax ont_16s_umi_asvf.fasta \
  --db silva_poor.fasta \
  --sintax_cutoff 0.80 \
  --tabbedout ont_16s_umi_asv_poor_tax.tsv
```

Take a look at the output files. We know we won't get genus level classifications as we removed all the correct genus level references. However, we should be able to get the correct family classifications as the closest family relatives should still be in the databases.

* Are the family classifications correct?

* What is the most disturbing issue with the classification?

* What kind of real world scenarios would we see similar effect?

* If we know/suspect that the databases are incomplete, what can we use the results for?

  

#### Classification with home made reference sequences

If you are working with poorly explored ecosystems there is a big risk that the targets of your analysis are underrepresented or completely missing from the databases. For microbial profiling with 16S amplicon analysis a good example is usually soil or sediment samples.

What can we do about that? Make references and taxonomy for our ecosystem of interest ourselves. For 16S taxonomy you have to use full length 16S sequences to be able to make decent taxonomy. This is where Nanopore UMI sequencing really has its advantage, as it can sequence the full length and with a high accuracy. The next issue is to assign taxonomy and this is no simple matter, and has historically been performed manually by experts with years of domain knowledge. However, frameworks are being made for creating de novo taxonomy automatically (e.g. AutoTax). This is not perfect but it creates a crude taxonomy that can be used until better taxonomy is defined. 

Let us imagine we are domain experts and we have generated long 16S sequences for our ecosystems in the form of `ont_16s_umi_asvf.fasta`. Based on either an automatic pipeline or by making a phylogenetic tree and manually defining groups, we have been able to assign taxonomic levels to the different ASVs. 

1. Make a copy of `ont_16s_umi_asvf.fasta`.

   ```
   cp ont_16s_umi_asvf.fasta custom_db.fasta
   ```

   

2. Open it in a text editor and add taxonomic levels in the format seen below. Use the `ont_16s_umi_tax.fasta` to guide you to make the custom levels, but provide your own name for each of the different levels. If you are pressed for time got to step 3.

   ```
   >asv1;size=4878;tax=f:f_ecoli,g:g_ecoli,s:s_ecoli1;
   GATGAACGCTGGCGGCGTGCCTAATACATGCAAGTCGAGCGAACGGACGAGAAGCTTGCTTCTCTGATGTTAGCGGCGGAC.....
   .....
   
   ```

3. Or download an example of manual curation here.

   ```
   wget \
     -O custom_db.fasta \
     "https://www.dropbox.com/s/3lsvpq36kqdba9h/custom_db.fasta?dl=1"
   ```

4. Now let us try to make classifications with our new database.

   ```
   # Illumina V4
   vsearch \
     --sintax ilm_v4_asvf.fasta \
     --db custom_db.fasta \
     --sintax_cutoff 0.80 \
     --tabbedout ilm_v4_asv_custom_tax.tsv
     
   # Nanopore raw 16S
   vsearch \
     --sintax ont_16s_raw.fasta \
     --db custom_db.fasta \
     --sintax_cutoff 0.80 \
     --tabbedout ont_16s_raw_custom_tax.tsv
     
   gawk '
     {TAX[$4";"]++}
     END{
       print "ont_16s_raw\ttaxonomy"
       for (i in TAX){
         print i"\t"TAX[i]
       }
     }
   ' ont_16s_raw_custom_tax.tsv \
   > ont_16s_raw_custom_tax_countable.tsv
   
   # Nanopore UMI 16S
   vsearch \
     --sintax ont_16s_umi_asvf.fasta \
     --db custom_db.fasta \
     --sintax_cutoff 0.80 \
     --tabbedout ont_16s_umi_asv_custom_tax.tsv
   ```



Take a look at the output files. The ecosystem specific database seems to perform well for our data. Is the solution just to make ecosystem specific databases?

* What is the problem of populating your own reference database with sequences? What happens if you are missing references?

* What is the problem of generating you own taxonomy or groupings? What is the benefits of doing it any way?

* What happens if you add your new references to the `silva_poor.fasta` database and perform the classification again? Merge the databases by doing:

  ```bash
  cat silva_poor.fasta custom_db.fasta > silva_custom.fasta
  ```



