//TODO: rebuild database but with taxids or full names appended 

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
process Low_complexity_filtering { 
//conda "${baseDir}/env/env.yml"
publishDir "${params.OUTPUT}/fastp_PE/${base}", mode: 'symlink', overwrite: true
container "quay.io/biocontainers/bbmap:38.76--h516909a_0"
beforeScript 'chmod o+rw .'
cpus 6
input: 
    tuple val(base), file(r1), file(r2)
output: 
    tuple val(base), file("${base}.lcf_filtered.R1.fastq.gz"), file("${base}.lcf_filtered.R2.fastq.gz")

script:
"""
#!/bin/bash
#logging
echo "ls of directory" 
ls -lah 

bbduk.sh \
    in1=${r1} in2=${r2} \
    out1=${base}.lcf_filtered.R1.fastq.gz out2=${base}.lcf_filtered.R2.fastq.gz \
    entropy=0.7 \
    entropywindow=50 \
    entropyk=4 
"""
}

//TODO: delete intermediate bam for space savings
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
    tuple val("${base}"), file("${base}.species.fasta.gz")


script:
"""
#!/bin/bash
#logging
echo "ls of directory" 
ls -lah 

# python3 ${extract_script} ${report} ${fastadb}

#grep -P "\tG\t" ${report} | cut -f5 | parallel {}.genus.fasta.gz /scratch/vpeddu/genus_level_download/test_index/

for i in `grep -P "\tG\t" ${report} | cut -f5`
do
echo adding \$i
if [[ -f ${fastadb}/\$i.genus.fasta.gz ]]; then
    cat ${fastadb}/\$i.genus.fasta.gz >> species.fasta.gz
fi
done


mv species.fasta.gz ${base}.species.fasta.gz

"""
}

//TODO: create containre with Minimap2 and samtools so we can get rid of intermediate sam for space savings
process Minimap2 { 
//conda "${baseDir}/env/env.yml"
publishDir "${params.OUTPUT}/Minimap2/${base}", mode: 'symlink'
container "staphb/minimap2"
beforeScript 'chmod o+rw .'
cpus 8
input: 
    tuple val(base), file(r1), file(r2), file(species_fasta)
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
    -K 16G \
    --split-prefix \
    -2 \
    ${species_fasta} \
    ${r1} ${r2} > \
    ${base}.minimap2.sam
"""
}

process Sam_conversion { 
publishDir "${params.OUTPUT}/sam_conversion/${base}", mode: 'symlink', overwrite: true
container "staphb/samtools"
beforeScript 'chmod o+rw .'
cpus 8
input: 
    tuple val(base), file(sam)
output: 
    tuple val("${base}"), file("${base}.sorted.bam"), file("${base}.sorted.bam.bai")
    file "${base}.unclassfied.bam"

script:
"""
#!/bin/bash
#logging
echo "ls of directory" 
ls -lah 

samtools view -Sb -@  ${task.cpus} -F 4 ${sam} > ${base}.bam
samtools sort -@ ${task.cpus} ${base}.bam > ${base}.sorted.bam
samtools index ${base}.sorted.bam

samtools view -Sb -@  ${task.cpus} -f 4 ${sam} > ${base}.unclassfied.bam


"""
}

//TODO: ADD ACCESSION DNE OUTPUT BACK IN 
process Classify { 
publishDir "${params.OUTPUT}/Classification/${base}", mode: 'symlink', overwrite: true
container 'quay.io/vpeddu/evmeta'
beforeScript 'chmod o+rw .'
errorStrategy 'ignore'
cpus 8
input: 
    tuple val(base), file(bam), file(bamindex), file(unclassified_bam), file(unclassified_fastq)
    file taxdump
    file classify_script
    file accessiontotaxid
output: 
    tuple val("${base}"), file("${base}.prekraken.tsv")
    //file "${base}.accession_DNE.txt"

script:
"""
#!/bin/bash
#logging
echo "ls of directory" 
ls -lah 
#mv taxonomy/taxdump.tar.gz .
#tar -xvzf taxdump.tar.gz
cp taxdump/*.dmp .
python3 ${classify_script} ${bam} ${base} ${accessiontotaxid}/nucl_gb.accession2taxid

samtools sort ${unclassified_bam} -o ${base}.unclassified.sorted.bam
samtools index ${base}.unclassified.sorted.bam
echo -e "0\\t `samtools view -c ${base}.unclassified.sorted.bam`"  >> ${base}.prekraken.tsv


"""
}

process Classify_orthologs { 
publishDir "${params.OUTPUT}/Classify_orthologs/${base}", mode: 'symlink', overwrite: true
container 'quay.io/vpeddu/evmeta'
beforeScript 'chmod o+rw .'
errorStrategy 'ignore'
cpus 8
input: 
    tuple val(base), file(bam), file(bamindex), file(unclassified_bam), file(unclassified_fastq)
    file taxdump
    file classify_script
    file accessiontotaxid
output: 
    tuple val("${base}"), file("${base}.orthologs.prekraken.tsv")
    //file "${base}.accession_DNE.txt"

script:
"""
#!/bin/bash
#logging
echo "ls of directory" 
ls -lah 
#mv taxonomy/taxdump.tar.gz .
#tar -xvzf taxdump.tar.gz
cp taxdump/*.dmp .
python3 ${classify_script} ${base}.emapper.annotations ${base} ${accessiontotaxid}/nucl_gb.accession2taxid

"""
}
process Write_report { 
publishDir "${params.OUTPUT}/", mode: 'symlink', overwrite: true
container "evolbioinfo/krakenuniq:v0.5.8"
beforeScript 'chmod o+rw .'
errorStrategy 'ignore'
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