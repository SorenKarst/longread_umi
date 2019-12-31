if [[ -z $cpus ]];
then
    echo "Please set the number of CPUs in \$cpus"
    exit 1
fi

if [[ -z $id ]];
then
    echo "Please set the sample ID of fasta to process in \$id"
    exit 1
fi

if [[ -z $db ]];
then
    echo "Please set the database path in \$db"
    exit 1
fi

blastn -num_threads $cpus -query $id.fa -db $db -evalue 1e-5 -max_target_seqs 1000 -outfmt 6 -out $id.out
