# aPOMP: a Pore Metagenomic Pipeline

## Introduction
APOMP is a portable metagenomics pipeline designed for use with Oxford Nanopore long read sequencing data. 
## Database download 
The full aPOMP database is 427GB and can be downloaded at <insert tarball link>. Smaller indexes are avaiable at <build smaller indexes>
* HG38 human host depletion (Minimap2, STAR)
* Full BLAST NT (downloaded 09/2021)
* Kraken PFP (downloaded 10/2021) 
* NCBI taxonomy (downloaded 09/20021)
* Eggnog DB (downloaded 12/2021)

## Nanopore workflow
run command:  
```
  nextflow run vpeddu/ev-meta \
    --NANOPORE \
	--INPUT_FOLDER <Input Folder> \
	--OUTPUT new_index_test \
	--INDEX <Index Location>  \
	--NUCL_TYPE <RNA or DNA> \
	-with-docker ubuntu:18.04 
    -latest \
    -with-report \
```
## Illumina workflow 
run command:  
```
  nextflow run vpeddu/ev-meta \
    --ILLUMINA \
	--INPUT_FOLDER <Input Folder> \
	--OUTPUT new_index_test \
	--INDEX <Index Location>  \
	--NUCL_TYPE <RNA or DNA> \
	-with-docker ubuntu:18.04 
    -latest \
    -with-report \
```