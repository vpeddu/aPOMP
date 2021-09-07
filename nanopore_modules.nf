process NanoFilt { 
//conda "${baseDir}/env/env.yml"
publishDir "${params.OUTPUT}/Nanofilt/${base}", mode: 'symlink', overwrite: true
container " quay.io/biocontainers/nanofilt:2.8.0"
beforeScript 'chmod o+rw .'
cpus 6
input: 
    tuple val(base), file(r1)
output: 
    tuple val(base), file("${base}.filtered.fastq.gz")
    file "*"


script:
"""
#!/bin/bash
#logging
echo "ls of directory" 
ls -lah 
echo "running Nanofilt on ${base}"

gunzip ${r1} | \
    NanoFilt -q 9 \
        --maxlength 5000 \
        --length 200 | \
        gzip > ${base}.filtered.fastq.gz

"""
}

process NanoPlot { 
//conda "${baseDir}/env/env.yml"
publishDir "${params.OUTPUT}/NanoPlot/${base}", mode: 'symlink', overwrite: true
container "quay.io/biocontainers/nanoplot:1.38.1--pyhdfd78af_0"
beforeScript 'chmod o+rw .'
cpus 2
input: 
    tuple val(base), file(r1) 
output: 
    "*"

script:
"""
#!/bin/bash
#logging
echo "ls of directory" 
ls -lah 

NanoPlot -t ${task.cpus} \
    -p ${base} \
    --fastq ${r1} \
    --title ${base} \
    --summary 
"""
}

process Host_depletion_nanopore { 
publishDir "${params.OUTPUT}/Host_filtered/${base}", mode: 'symlink', overwrite: true
container "staphb/minimap2"
beforeScript 'chmod o+rw .'
cpus 8
input: 
    tuple val(base), file(r1), file(r2)
    file minimap2_host_index
output: 
    tuple val("${base}"), file("${base}.starUnmapped.out.mate1.fastq.gz")
script:
"""
#!/bin/bash
#logging
echo "ls of directory" 
ls -lah 

minimap2 \
    -ax map-ont \
    -t ${task.cpus} \
    -2 \
    ${minimap2_host_index}
    ${r1} > \
    ${base}.host_filtered.sam
"""
}

process Host_depletion_extraction_nanopore { 
publishDir "${params.OUTPUT}/Host_filtered/${base}", mode: 'symlink', overwrite: true
container "staphb/samtools"
beforeScript 'chmod o+rw .'
cpus 8
input: 
    tuple val(base), file(sam)
output: 
    tuple val("${base}"), file("${base}.host_filtered.fastq.gz")
script:
"""
#!/bin/bash
#logging
echo "ls of directory" 
ls -lah 

samtools fastq -n -f 4 ${sam} | gzip > ${base}.host_filtered.fastq.gz
"""
}

process Kraken_prefilter { 
publishDir "${params.OUTPUT}/Interleave_FASTQ/${base}", mode: 'symlink', overwrite: true
container "staphb/kraken2"
beforeScript 'chmod o+rw .'
cpus 8
input: 
    tuple val(base), file(r1)
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
    --threads ${task.cpus} \
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
container 'quay.io/vpeddu/evmeta'
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

python3 ${extract_script} ${report} ${fastadb}

mv species.fasta ${base}.species.fasta

"""
}

process Minimap2 { 
//conda "${baseDir}/env/env.yml"
publishDir "${params.OUTPUT}/Minimap2/${base}", mode: 'symlink'
container "staphb/minimap2"
beforeScript 'chmod o+rw .'
cpus 8
input: 
    tuple val(base), file(r1), file(r2)
    tuple val(base), file(species_fasta)
output: 
    tuple val("${base}"), file("${base}.minimap2.sam")

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

process Sam_conversion { 
publishDir "${params.OUTPUT}/Profiling/${base}", mode: 'symlink', overwrite: true
container "staphb/samtools"
beforeScript 'chmod o+rw .'
cpus 8
input: 
    tuple val(base), file(sam)
output: 
    tuple val("${base}"), file("${base}.sorted.bam"), file("${base}.sorted.bam.bai")

script:
"""
#!/bin/bash
#logging
echo "ls of directory" 
ls -lah 

samtools view -Sb -@  ${task.cpus} -F 4 ${sam} > ${base}.bam
samtools sort -@ ${task.cpus} ${base}.bam > ${base}.sorted.bam
samtools index ${base}.sorted.bam
"""
}

process Classify { 
publishDir "${params.OUTPUT}/Profiling/${base}", mode: 'symlink', overwrite: true
container 'quay.io/vpeddu/evmeta'
beforeScript 'chmod o+rw .'
cpus 8
input: 
    tuple val(base), file(bam), file(bamindex)
    file taxdump
    file classify_script
output: 
    tuple val("${base}"), file("${base}.prekraken.tsv")

script:
"""
#!/bin/bash
#logging
echo "ls of directory" 
ls -lah 
#mv taxonomy/taxdump.tar.gz .
#tar -xvzf taxdump.tar.gz
cp viral/*.dmp .
python3 ${classify_script} ${bam} ${base}
"""
}

process Write_report { 
publishDir "${params.OUTPUT}/", mode: 'symlink', overwrite: true
container "evolbioinfo/krakenuniq:v0.5.8"
beforeScript 'chmod o+rw .'
cpus 8
input: 
    tuple val(base), file(prekraken)
    file krakenuniqdb
output: 
    file "${base}.final.report.tsv"

script:
"""
#!/bin/bash
#logging
echo "ls of directory" 
ls -lah 

krakenuniq-report --db ${krakenuniqdb} \
 --taxon-counts \
  ${prekraken} > ${base}.final.report.tsv
"""
}