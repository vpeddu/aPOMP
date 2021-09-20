
# arg 1 is nodes.dmp
grep 'genus' $1 | cut -f1 >> genus_taxids.txt

for i in `cat genus_taxids.txt`
do
    esearch -db genome -query "$i [Organism]"|elink -target nuccore|efilter -query "genbank"|efetch -format fasta >> $i.fasta
    pigz $i.fasta
done
