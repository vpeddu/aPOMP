#create filestructure
mkdir eggnog_db krakenuniq_db ribosome_trna star_host kraken2_db minimap2_host plasmid_db sourmash taxdump accession2taxid/

# nt download
`wget': seq -w 000 157 | \
	parallel -j 32 \
	wget https://ftp.ncbi.nlm.nih.gov/blast/db/nt.{}.tar.gz

# build nt genus lookup table
python3 accession_to_genus.py \
	nucl_gb.accession2taxid \
	/private/home/vpeddu/.etetoolkit/taxa.sqlite

# build nt index
python3 genus_to_fasta.py \
	nucl2gb_lookup.pkl \
	nt/
mkdir final_index/
ls genus_organized/ |\
	cut -f1 -d . |\
	sort |\
	uniq |\
	parallel -j 32 \
	"cat genus_organized/{}.nt.* > final_index/{}.genus.fasta.gz"

mv nucl_gb.accession2taxid accession2taxid/

# for gtnadb go to https://gtrnadb.ucsc.edu/search.html, enter a blank query, and download all
pigz gtrnadb-*

# Download ribosome 
wget  ftp://ftp.ncbi.nlm.nih.gov/refseq/TargetedLoci/Archaea/*.fna.gz
wget  ftp://ftp.ncbi.nlm.nih.gov/refseq/TargetedLoci/Bacteria/*.fna.gz
wget  ftp://ftp.ncbi.nlm.nih.gov/refseq/TargetedLoci/Fungi/*.fna.gz
cat *rRNA.fna.gz > ribosome.fa.gz

cat gtrnadb-* ribosome.fa.gz > ribosome_trna.combined.fa.gz
seqkit rmdup -s ribosome_trna.combined.fa.gz | pigz > all_trna.fa.gz

# download kraken2 pluspfp16 
aws s3 cp --no-sign-request s3://genome-idx/kraken/k2_pluspfp_16gb_20240605.tar.gz kraken2_db
tar -xvzf kraken2_db/k2_pluspfp_16gb_20240605.tar.gz

# download minikraken krakenuniq database (this is very old)
wget https://ccb.jhu.edu/software/kraken/dl/minikraken_20171019_4GB.tgz
tar -xvzf minikraken_20171019_4GB.tgz
mv minikraken_20171013_4GB*/* krakenuniq_db/
rm -rf minikraken_20171013_4GB* 

# build plsdb
# download from https://ccb-microbe.cs.uni-saarland.de/plsdb/plasmids/download/
bzip2 -d plsdb.fna.bz2

# downlaod and extract taxdump
wget https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz
tar -xvzf taxdump.tar.gz
