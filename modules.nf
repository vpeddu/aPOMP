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

process Interleave_FASTQ { 
publishDir "${params.OUTPUT}/Star_PE/${base}", mode: 'symlink', overwrite: true
container "staphb/bbtools"
beforeScript 'chmod o+rw .'
cpus 8
input: 
    tuple val(base), file(r1), file(r2)
output: 
    tuple val("${base}"), file("${base}.unmapped.interleaved.fastq.gz")
script:
"""
#!/bin/bash
#logging
echo "ls of directory" 
ls -lah 

reformat.sh in1=${r1} in2=${r2} out1=${base}.unmapped.interleaved.fastq.gz
"""
}

process Metalign_db_selection { 
publishDir "${params.OUTPUT}/Star_PE/${base}", mode: 'symlink', overwrite: true
//container "quay.io/biocontainers/star:2.7.9a--h9ee0642_0"
conda 'Metalign'
beforeScript 'chmod o+rw .'
cpus 8
input: 
    tuple val(base), file(interleaved_fastq)
    file metalign_db
output: 
    file "${base}.metalign.tempdir"



script:
"""
#!/bin/bash
#logging
echo "ls of directory" 
ls -lah 

select_db.py ${interleaved_fastq} ${metalign_db} \
    --keep_temp_files \
    --strain_level
    --temp_dir ${base}.metalign.tempdir

"""
}

process Minimap2 { 
//conda "${baseDir}/env/env.yml"
publishDir "${params.OUTPUT}/Minimap2/${base}", mode: 'symlink'
container "biocontainers/minimap2:v2.15dfsg-1-deb_cv1"
beforeScript 'chmod o+rw .'
cpus 4
input: 
    tuple val(base), file(r1)
    file metalign_db
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
    -ax splice \
    -t ${task.cpus} \
    -L \
    -Y \
    -2 \
    ${metalign_db}/cmashed_db.fna \
    ${r1} > \
    ${base}.minimap2.sam
"""
}
