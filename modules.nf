process Trimming_FastP { 
//conda "${baseDir}/env/env.yml"
publishDir "${params.OUTPUT}/fastp_PE/${base}", mode: 'symlink', overwrite: true
container "bromberglab/fastp"
beforeScript 'chmod o+rw .'
cpus 6
input: 
    tuple val(base), file(r1), file(r2)
output: 
    tuple val(base), file("${base}.trimmed.R1.fastq.gz"), file("${base}.trimmed.R2.fastq.gz")
    tuple val(base), file("${base}.trimmed.R1.fastq.gz")
    file "*"


script:
"""
#!/bin/bash
#logging
echo "ls of directory" 
ls -lah 
echo "running fastp on ${base}"
fastp -w ${task.cpus} \
    -i ${r1} \
    -I ${r2} \
    -o ${base}.trimmed.R1.fastq.gz \
    -O ${base}.trimmed.R2.fastq.gz
"""
}

process Host_depletion { 
publishDir "${params.OUTPUT}/Star_PE/${base}", mode: 'symlink', overwrite: true
container "quay.io/biocontainers/star:2.7.9a--h9ee0642_0"
beforeScript 'chmod o+rw .'
cpus 8
input: 
    tuple val(base), file(r1), file(r2)
    file starindex
output: 
    file "${base}.star*"
    file "${base}.starAligned.out.bam"
    tuple val("${base}"), file("${base}.starUnmapped.out.mate1.fastq.gz"), file("${base}.starUnmapped.out.mate2.fastq.gz")
script:
"""
#!/bin/bash
#logging
echo "ls of directory" 
ls -lah 
STAR   \
    --runThreadN ${task.cpus}  \
    --genomeDir ${starindex}   \
    --readFilesIn ${r1} ${r2} \
    --readFilesCommand zcat      \
    --outSAMtype BAM Unsorted \
    --outReadsUnmapped Fastx \
    --outFileNamePrefix ${base}.star  

mv ${base}.starUnmapped.out.mate1 ${base}.starUnmapped.out.mate1.fastq
mv ${base}.starUnmapped.out.mate2 ${base}.starUnmapped.out.mate2.fastq

gzip ${base}.starUnmapped.out.mate1.fastq
gzip ${base}.starUnmapped.out.mate2.fastq
"""
}

process Kraken_prefilter { 
publishDir "${params.OUTPUT}/Interleave_FASTQ/${base}", mode: 'symlink', overwrite: true
container "staphb/kraken2"
beforeScript 'chmod o+rw .'
cpus 8
input: 
    tuple val(base), file(r1), file(r2)
    file kraken2_db
output: 
    tuple val("${base}"), file("${base}.kraken2.report")
script:
"""
#!/bin/bash
#logging
echo "ls of directory" 
ls -lah 

kraken2 --db ${kraken2_db} \
    --threads 24 \
    --classified-out ${base}.kraken2.classified \
    --output ${base}.kraken2.output \
    --report ${base}.kraken2.report \
    --gzip-compressed \
    --unclassified-out ${base}.kraken2.unclassified \
    ${r1} ${r2}

"""
}

process Extract_db { 
//publishDir "${params.OUTPUT}//${base}", mode: 'symlink', overwrite: true
//container "quay.io/biocontainers/star:2.7.9a--h9ee0642_0"
conda 'Metalign'
beforeScript 'chmod o+rw .'
cpus 8
input: 
    tuple val(base), file(report)
    file fastadb
    file extract_script
output: 
    tuple val("${base}"), file("${base}.species.fasta")


script:
"""
#!/bin/bash
#logging
echo "ls of directory" 
ls -lah 

python3 ${extract_script} ${report}

mv species.fasta ${base}.species.fasta

"""
}

process Minimap2 { 
//conda "${baseDir}/env/env.yml"
publishDir "${params.OUTPUT}/Minimap2/${base}", mode: 'symlink'
container "biocontainers/minimap2:v2.15dfsg-1-deb_cv1"
beforeScript 'chmod o+rw .'
cpus 8
input: 
    tuple val(base), file(r1), file(r2)
    tuple val(base), file(species_fasta)
output: 
    file "${base}.minimap2.aligned.sam"

script:
"""
#!/bin/bash

#logging
echo "ls of directory" 
ls -lah 

echo "running Minimap2 on ${base}"
#TODO: FILL IN MINIMAP2 COMMAND 
minimap2 \
    -ax sr \
    -t ${task.cpus} \
    -2 \
    ${species_fasta} \
    ${r1} ${r2} > \
    ${base}.minimap2.sam
"""
}

process Metalign_profiling { 
publishDir "${params.OUTPUT}/Profiling/${base}", mode: 'symlink', overwrite: true
container "biocontainers/samtools"
conda 'Metalign'
beforeScript 'chmod o+rw .'
cpus 8
input: 
    tuple val(base), file(sam)
    file metalign_db
output: 
    //file "${base}.metalign.tempdir"

script:
"""
#!/bin/bash
#logging
echo "ls of directory" 
ls -lah 

map_and_profile.py ${interleaved_fastq} ${metalign_db} \
    --length_normalize \
    --strain_level \
    --threads ${task.cpus} \
    --output ${base}.abundances.tsv \
    --verbose

"""
}