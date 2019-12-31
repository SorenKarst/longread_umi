if [[ -z $rnammer ]];
then
    echo "Please set the path to RNAmmer in \$rnammer"
    exit 1
fi

if [[ -z $id ]];
then
    echo "Please set the sample ID of fasta to process in \$id"
    exit 1
fi

$rnammer -S bac -m lsu,ssu,tsu -gff $id.gff $id.fa
