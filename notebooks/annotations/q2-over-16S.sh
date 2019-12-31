source activate qiime2-2019.10

if [[ -z $sequence_dir ]];
then
    echo "Please specify the directory containing the full length"
    echo "16S sequences in the \$sequence_dir variable"
    exit 1
fi

find . -name "${sequence_dir}/sequences.fa" -exec cat {} \; > q2_seqs_for_import.fa
qiime tools import --input-path q2_seqs_for_import.fa --output-path q2_seqs.qza --type SampleData[Sequences]
qiime vsearch dereplicate-sequences \
    --i-sequences q2_seqs.qza \
    --o-dereplicated-table q2_seqs_derep_table.qza \
    --o-dereplicated-sequences q2_seqs_derep.qza

qiime vsearch cluster-features-de-novo \
  --i-table q2_seqs_derep_table.qza \
  --i-sequences q2_seqs_derep.qza \
  --p-perc-identity 0.99 \
  --o-clustered-table q2_seqs-table-dn-99.qza \
  --o-clustered-sequences q2_seqs-seqs-dn-99.qza \
  --p-threads 4

qiime feature-classifier classify-consensus-vsearch \
    --i-query q2_seqs-seqs-dn-99.qza \
    --i-reference-reads ~/ResearchWork/greengenes_release/gg_13_8_otus/99_repset.qza \
    --i-reference-taxonomy ~/ResearchWork/greengenes_release/gg_13_8_otus/99_taxonomy.qza \
    --p-threads 4 \
    --p-strand plus \
    --p-query-cov 0.9 \
    --p-perc-identity 0.9 \
    --o-classification q2_seqs-seqs-dn-99-classified.qza
