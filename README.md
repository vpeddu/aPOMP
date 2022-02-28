# aPOMP: a Pore Metagenomic Pipeline

## Introduction
APOMP is a portable metagenomics pipeline designed for use with Oxford Nanopore long read sequencing data. 

## Workflow 
aPOMP options 

# *required flags
flag| description
-----|-----
*--NANOPORE/--ILLUMINA|sequencing platform used
*--INPUT\_FOLDER|folder containing input FASTQs (illumina must be paired)  
*--INDEX|path to aPOMP index folder
*--OUTPUT|output folder
--IDENTIFY\_RESISTANCE\_PLASMIDS|assemble and identify resistant plasmids 
--CLEAN\_RIBOSOME\_TRNA|filter out ribosomal and tRNA sequences before classification 
--EGGNOG|identify orthologous groups in unclassified reads using  the Eggnog-mapper
--MINIMAPSPLICE|run Minimap2 with preset -ax splice (default -ax map-ont)
--NANOFILT\_QUALITY [int]|minimum quality threshold for nanofilt (default 10) 
--KRAKEN2\_THRESHOLD [int]|discard Kraken2 results containing less than this many reads at or below the genus level (default 10)
--LOW_COMPLEXITY_FILTER_NANOPORE|run low complexity filtering on nanopore samples
--METAFLYE|run metaflye before all classifications. Not for use on small memory machines
--NANOFILT_MINLENGTH [int]| Have Nanofilt filter out any reads smaller than this number 
--NANOFILT_MAXLENGTH [int]| Have Nanofilt filter out any reads larger than this number 

--help|display help message
## Database download 
The full aPOMP database is ~250GB and can be downloaded at <insert tarball link>. Smaller indexes are avaiable at <build smaller indexes>
* HG38 human host depletion (Minimap2, STAR)
* Full BLAST NT (downloaded 09/2021)
* Kraken PFP (downloaded 10/2021) 
* NCBI taxonomy (downloaded 09/20021)
* Eggnog DB (downloaded 12/2021)

## Nanopore workflow
quick run command:  
```
  nextflow run vpeddu/ev-meta \		
		 --NANOPORE \
		 --IDENTIFY_RESISTANCE_PLASMIDS \
   -profile standard \
		 --EGGNOG \
		 --INPUT_FOLDER <input_fastq_folder> \
		 --OUTPUT Zymo-GridION-EVEN-BB-SN_out \
		 --INDEX <index_path> \
		 --CLEAN_RIBOSOME_TRNA \
		 -with-docker 'ubuntu:18.04' \
		 -with-tower \
		 -with-report \
		 -latest \
		 -resume
```
## Illumina workflow (need to update)
quick run command:  
```
  nextflow run vpeddu/ev-meta \
    --ILLUMINA \
	--INPUT_FOLDER <Input folder> \
	--OUTPUT <Output folder> \
	--INDEX <Index location>  \
	--NUCL_TYPE <RNA or DNA> \
	-with-docker ubuntu:18.04 
    -latest \
    -with-report 
```

