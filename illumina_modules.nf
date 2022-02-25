params.FASTA_SPLIT_CHUNKS = 10
params.KRAKEN2_THRESHOLD = 10
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

//BIG BUG: if kraken can't classify past the phylum level we aren't 
// pulling any genuses so those reads are left unclassified
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
    file("${base}__*")


script:
"""
#!/bin/bash
#logging
echo "ls of directory" 
ls -lah 

# python3 ${extract_script} ${report} ${fastadb}

#grep -P "\tG\t" ${report} | cut -f5 | parallel {}.genus.fasta.gz /scratch/vpeddu/genus_level_download/test_index/
# could filter by kraken report column 2 for all above some parameter (if > 25)
for i in `grep -P "\tG\t" ${report} | awk '\$2>${params.KRAKEN2_THRESHOLD}' | cut -f5`
do
echo adding \$i
if [[ -f ${fastadb}/\$i.genus.fasta.gz ]]; then
    ##cat ${fastadb}/\$i.genus.fasta.gz >> species.fasta.gz
    cp ${fastadb}/\$i.genus.fasta.gz ${base}__\$i.genus.fasta.gz
fi
done

# TODO need to optimize this 
##mv species.fasta.gz ${base}.species.fasta.gz

##gunzip ${base}.species.fasta.gz

##pyfasta split -n ${params.FASTA_SPLIT_CHUNKS} ${base}.species.fasta

##pigz ${base}.species.fasta 


##for f in *.fasta; do mv "\$f" "${base}__-\$f"; done

##for i in *.fasta; do pigz \$i; done
"""
}

process Minimap2_illumina { 
//conda "${baseDir}/env/env.yml"
publishDir "${params.OUTPUT}/Minimap2/${base}", mode: 'symlink'
container "quay.io/vpeddu/evmeta:latest"
beforeScript 'chmod o+rw .'
cpus 8
input: 
    tuple val(base), file(r1), file(r2), file(species_fasta)
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
minimap2 \
    -ax sr \
    -t ${task.cpus} \
    -K 16G \
    --split-prefix \
    -2 \
    ${species_fasta} \
    ${r1} ${r2} | samtools view -Sb -@ 4 - > ${base}.bam

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
//TODO: CHANGE TO LCA* 
//TODO: do something about the stupid cp taxdump/*.dmp step
process Classify { 
publishDir "${params.OUTPUT}/Classification/${base}", mode: 'symlink', overwrite: true
container 'vpeddu/nanopore_metagenomics:latest'
beforeScript 'chmod o+rw .'
errorStrategy 'ignore'
cpus 8
input: 
    tuple val(base), file(bam), file(bamindex), file(unclassified_fastq), file(plasmid_fastq), val(plasmid_count)
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
python3 ${classify_script} ${bam} ${base} 

# counting unassigned reads to add back into final report
#samtools sort -@ ${task.cpus}-o ${base}.unclassified.sorted.bam
#samtools index ${base}.unclassified.sorted.bam
#echo -e "0\\t `samtools view -c ${base}.unclassified.sorted.bam`"  >> ${base}.prekraken.tsv

 echo \$(zcat ${unclassified_fastq} | wc -l)/4 | bc >> ${base}.prekraken.tsv
 linecount=\$(zcat ${unclassified_fastq} | wc -l)
 fastqlinecount=\$(echo \$linecount/4|bc)
 echo -e "0\\t\$fastqlinecount" >> ${base}.prekraken.tsv

echo -e "36549\\t${plasmid_count}" >> ${base}.prekraken.tsv

echo \$fastqlinecount unclassified reads 
"""
}

process Classify_orthologs { 
publishDir "${params.OUTPUT}/Classify_orthologs/${base}", mode: 'symlink', overwrite: true
container 'quay.io/vpeddu/evmeta'
beforeScript 'chmod o+rw .'
errorStrategy 'ignore'
cpus 24
input: 
    // need to fix input cardinality coming from eggnog
    // need to calculate number of still unassigned reads and output those too 
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
python3 ${classify_script} ${base}.emapper.annotations ${base} 

"""
}

process Write_report { 
publishDir "${params.OUTPUT}/", mode: 'symlink', overwrite: true
container "evolbioinfo/krakenuniq:v0.5.8"
beforeScript 'chmod o+rw .'
errorStrategy 'ignore'
cpus 1
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

process Write_report_orthologs { 
publishDir "${params.OUTPUT}/ortholog_reports/", mode: 'symlink', overwrite: true
container "evolbioinfo/krakenuniq:v0.5.8"
beforeScript 'chmod o+rw .'
errorStrategy 'ignore'
cpus 1
input: 
    tuple val(base), file(prekraken)
    file krakenuniqdb
output: 
    file "${base}.orthologs.final.report.tsv"

script:
"""
#!/bin/bash
#logging
echo "ls of directory" 
ls -lah 

krakenuniq-report --db ${krakenuniqdb} \
--taxon-counts \
${prekraken} > ${base}.orthologs.final.report.tsv
"""
}