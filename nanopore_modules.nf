process NanoFilt { 
//conda "${baseDir}/env/env.yml"
publishDir "${params.OUTPUT}/Nanofilt/${base}", mode: 'symlink', overwrite: true
container " quay.io/biocontainers/nanofilt:2.8.0--py_0"
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

gunzip -c ${r1} | NanoFilt -q 9 \
        --maxlength 5000 \
        --length 200 | gzip > ${base}.filtered.fastq.gz

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
    file "*"

script:
"""
#!/bin/bash
#logging
echo "ls of directory" 
ls -lah 

NanoPlot -t ${task.cpus} \
    -p ${base} \
    --fastq ${r1} \
    --title ${base} 
"""
}

process Host_depletion_nanopore { 
publishDir "${params.OUTPUT}/Host_filtered/${base}", mode: 'symlink', overwrite: true
container "staphb/minimap2"
beforeScript 'chmod o+rw .'
cpus 8
input: 
    tuple val(base), file(r1)
    file minimap2_host_index
output: 
    tuple val("${base}"), file("${base}.host_filtered.sam")
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
    ${minimap2_host_index} \
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

//TODO flag for nano-corr or nano-hq
process MetaFlye { 
publishDir "${params.OUTPUT}/MetaFlye/${base}", mode: 'symlink', overwrite: true
container "quay.io/biocontainers/flye:2.9--py27h6a42192_0"
beforeScript 'chmod o+rw .'
errorStrategy 'ignore'
cpus 16
input: 
    tuple val(base), file(unassigned_bam), file(unassigned_fastq)
output: 
    tuple val("${base}"), file("${base}.flye.fasta"), file("${unassigned_fastq}")
script:
"""
#!/bin/bash
#logging
echo "ls of directory" 
ls -lah 

flye --nano-corr ${unassigned_fastq} \
    --out-dir ${base}.flye \
    -t ${task.cpus} \
    --meta 

mv ${base}.flye/assembly.fasta ${base}.flye.fasta

"""
}


process Kraken_prefilter_nanopore { 
publishDir "${params.OUTPUT}/Kraken_prefilter/${base}", mode: 'symlink', overwrite: true
container "staphb/kraken2"
beforeScript 'chmod o+rw .'
cpus 8
input: 
    tuple val(base), file(flye_assembly)
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
    ${flye_assembly} 

"""
}


process Minimap2_nanopore { 
//conda "${baseDir}/env/env.yml"
publishDir "${params.OUTPUT}/Minimap2/${base}", mode: 'symlink'
container "quay.io/vpeddu/evmeta:latest"
beforeScript 'chmod o+rw .'
cpus 24
input: 
    tuple val(base), file(r1), file(species_fasta)
output: 
    tuple val("${base}"), file("${base}.sorted.filtered.bam"), file("${base}.sorted.filtered.bam.bai")
    tuple val("${base}"), file("${base}.unclassified.bam"), file ("${base}.unclassified.fastq.gz")

script:
"""
#!/bin/bash

#logging
echo "ls of directory" 
ls -lah 

echo "running Minimap2 on ${base}"
#TODO: FILL IN MINIMAP2 COMMAND 
minimap2 \
    -ax map-ont \
    -t "\$((${task.cpus}-4))" \
    -2 \
    --split-prefix \
    -K16G \
    ${species_fasta} \
    ${r1} | samtools view -Sb -@ 4 - > ${base}.bam

samtools view -Sb -F 4 ${base}.bam > ${base}.filtered.bam
samtools sort ${base}.filtered.bam -o ${base}.sorted.filtered.bam 
samtools index ${base}.sorted.filtered.bam
# output unclassified reads
samtools view -Sb -@  ${task.cpus} -f 4 ${base}.bam > ${base}.unclassified.bam

# cleanup intermediate file
rm ${base}.bam

samtools fastq -@ ${task.cpus} ${base}.unclassified.bam | gzip > ${base}.unclassified.fastq.gz

"""
}

// TODO: UPDATE INDEX SO WE CAN USE NEWEST VERSION OF DIAMOND
process Diamond_translated_alignment_unclassified { 
publishDir "${params.OUTPUT}/Diamond_unclassified_translated/${base}", mode: 'symlink', overwrite: true
container "quay.io/biocontainers/diamond:2.0.13--hdcc8f71_0"
beforeScript 'chmod o+rw .'
cpus 20
input: 
    tuple val(base), file(unassigned_bam), file(unassigned_fastq)
    //tuple val(base), file(unclassified_bam), file(unclassified_fastq)
    file diamond_protein_db
output: 
    tuple val("${base}"), file("*.diamond.out*")
script:
"""
#!/bin/bash
#logging
echo "ls of directory" 
ls -lah 

    # diamond out formats
	#0 = BLAST pairwise
	#5 = BLAST XML
	#6 = BLAST tabular
	#100 = DIAMOND alignment archive (DAA)
	#101 = SAM
if [[ -s ${unassigned_fastq} ]] 
    then
        echo "HERE"
        #https://currentprotocols.onlinelibrary.wiley.com/doi/full/10.1002/cpz1.59
        diamond blastx \
            --query ${unassigned_fastq} \
            --db ${diamond_protein_db} \
            --out ${base}.diamond.out \
            --outfmt 101 \
            --threads ${task.cpus} \
            --compress 1 \
            --unal 1 \
            --un ${base}.diamond.unaligned \
            --top 10 \
            -F 15 \
            --range-culling
    else
        echo "THERE"
        touch ${base}.diamond.out.blankinput

fi
"""
}

process Mmseq2_translated_alignment_unclassified { 
publishDir "${params.OUTPUT}/Mmsesq2_unclassified_translated/${base}", mode: 'symlink', overwrite: true
container "quay.io/biocontainers/mmseqs2:13.45111--h95f258a_1"
beforeScript 'chmod o+rw .'
cpus 8
input: 
    tuple val(base), file(unassigned_bam), file(unassigned_fastq)
    //tuple val(base), file(unclassified_bam), file(unclassified_fastq)
    file diamond_protein_db
output: 
    tuple val("${base}"), file("${base}.read_clu_rep.fasta"), file("${base}.mmseq.report.tsv")
script:
"""
#!/bin/bash
#logging
echo "ls of directory" 
ls -lah 

    # diamond out formats
	#0 = BLAST pairwise
	#5 = BLAST XML
	#6 = BLAST tabular
	#100 = DIAMOND alignment archive (DAA)
	#101 = SAM
if [[ -s ${unassigned_fastq} ]] 
    then
        echo "HERE"
        #https://currentprotocols.onlinelibrary.wiley.com/doi/full/10.1002/cpz1.59
        #diamond blastx \
        #    --query ${unassigned_fastq} \
        #    --db ${diamond_protein_db} \
        #    --out ${base}.diamond.out \
        #    --outfmt 101 \
        #    --threads ${task.cpus} \
        #   --compress 1 \
        #    --unal 1 \
        #    --un ${base}.diamond.unaligned \
        #    --top 10 \
        #    -F 15 \
        #    --range-culling

        # create read database
        mmseqs createdb ${unassigned_fastq} reads

        #cluster with linclust 
        mmseqs linclust reads reads_clu tmp 
        mmseqs createsubdb reads_clu reads reads_clu_rep 
        
        # extract clustered fasta
        mmseqs convert2fasta reads_clu_rep ${base}.read_clu_rep.fasta

        mmseqs taxonomy reads_clu_rep ${diamond_protein_db}/swissprot lca_result tmp -s 2 --threads ${task.cpus}
        mmseqs taxonomyreport ${diamond_protein_db}/swissprot  lca_result ${base}.mmseq.report.tsv 

    else
        echo "THERE"
        touch ${base}.diamond.out.blankinput

fi
"""
}

process Eggnog_mapper { 
publishDir "${params.OUTPUT}/Eggnog/${base}", mode: 'symlink', overwrite: true
container "quay.io/biocontainers/eggnog-mapper:2.1.6--pyhdfd78af_0"
beforeScript 'chmod o+rw .'
cpus 8
input: 
    tuple val(base), file(unassigned_bam), file(unassigned_fastq)
    //tuple val(base), file(unclassified_bam), file(unclassified_fastq)
    file eggnog_db
output: 
    tuple val("${base}"), file("*")
script:
"""
#!/bin/bash
#logging
echo "ls of directory" 
ls -lah 
if [[ -s ${unassigned_fastq} ]] 
    then
        echo "HERE"
        # convert unclassified fastq to fasta for eggnog
        gunzip -f ${unassigned_fastq} 
        sed -n '1~4s/^@/>/p;2~4p' ${base}.unclassified.fastq > ${base}.unclassified.fasta
        emapper.py \
            --dmnd_frameshift 15 \
            #-m mmseqs \
            --itype metagenome \
            -i ${base}.unclassified.fasta \
            --cpu ${task.cpus} \
            --data_dir ${eggnog_db} \
            -o test 
    else
        echo "THERE"
        touch ${base}.diamond.out.blankinput

fi
"""
}

process Extract_true_novel { 
//conda "${baseDir}/env/env.yml"
publishDir "${params.OUTPUT}/novel_reads/${base}", mode: 'symlink'
container "quay.io/vpeddu/evmeta:latest"
beforeScript 'chmod o+rw .'
cpus 24
input: 
    tuple val(base), file(unassigned_fastq), file(metaflye_contigs)
output: 
    tuple val("${base}"), file("${base}.unassembled.unclassified.fastq.gz")


script:
"""
#!/bin/bash

#logging
echo "ls of directory" 
ls -lah 

echo "remapping ${base} to contigs to find unassembled reads"
#TODO: FILL IN MINIMAP2 COMMAND 
minimap2 \
    -ax map-ont \
    -t "\$((${task.cpus}-4))" \
    -2 \
    --split-prefix \
    -K16G \
    ${unassigned_fastq} \
    ${metaflye_contigs}| samtools view -Sb -f 4 -@ 4 - > ${base}.unassembled.unclassified.bam

samtools fastq -@ ${task.cpus} ${base}.unassembled.unclassified.bam | gzip > ${base}.unassembled.unclassified.fastq.gz

"""
}
