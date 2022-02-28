# aPOMP: a Pore Metagenomic Pipeline

## Introduction
aPOMP is a portable metagenomics pipeline designed for use with Oxford Nanopore long read sequencing data. 

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
--KRAKEN2\_THRESHOLD [int]|discard Kraken2 results containing less than this many reads at or below the genus level (default 10)
--LOW_COMPLEXITY_FILTER_NANOPORE|run low complexity filtering on nanopore samples
--METAFLYE|run metaflye before all classifications. Not for use on small memory machines
--NANOFILT\_QUALITY [int]|minimum quality threshold for nanofilt (default 10) 
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
## Workflow 
### Read filtering 
1. Low complexity filtering (bbDuk.sh) if `--LOW_COMPLEXITY_FILTER_NANOPORE` specified
2. Read quality filtering (NanoFilt). 
..*Quality threshold adjustable with `--NANOFILT_QUALITY [int]` (default: 10) 
..*Min readlength adjustable with `--NANOFILT_MINLENGTH [int]` (default: 200)
..*Max readlength adjustable with `--NANOFILT_MAXLENGTH [int]` (default: 5000)
3. tRNA filtering (Minimap2) if `--CLEAN_RIBOSOME_TRNA` specified 
..* Reference tRNA database downloaded from http://gtrnadb.ucsc.edu/cgi-bin/GtRNAdb2-search.cgi
4. Host filtering (Minimap2, HG38 default). To specify a different host, replace the fasta in the index folder (/path_to_index/minimap2_host/new_host.fa) 
5. Plasmid extraction (Minimap2) done with alignment against plsDB v.2021_06_23_v2
..* if --IDENTIFY_RESISTANCE_PLASMIDS specified, plasmid reads are first assembled (`Flye --plasmid`), and then run against NCBI AMRfinder 

![alt text](https://github.com/vpeddu/ev-meta/blob/main/img/read_filtering.png)

### Alignment to NT 
![alt text](https://github.com/vpeddu/ev-meta/blob/main/img/alignment.png)
### Classification 
![alt text](https://github.com/vpeddu/ev-meta/blob/main/img/classification.png)

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

