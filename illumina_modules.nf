// Ignore genuses with less than this many reads in prefiltering
params.KRAKEN2_THRESHOLD = 10

//trim short reads with fastp
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

// run low complexity filtering on illumina samples
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


process Host_depletion_illumina { 
publishDir "${params.OUTPUT}/Host_filtered/${base}", mode: 'symlink', overwrite: true
container "vpeddu/nanopore_metagenomics:latest"
beforeScript 'chmod o+rw .'
cpus 8
input: 
    tuple val(base), file(r1), file(r2)
    file minimap2_host_index
    file ribosome_trna
    file minimap2_plasmid_db
output: 
    tuple val("${base}"), file("${base}.host_filtered_R1.fastq.gz") file("${base}.host_filtered_R2.fastq.gz")
    file "${base}.host_mapped.bam"
    file "${base}.trna.mapped.bam"
    tuple val("${base}"), file("${base}.plasmid.fastq.gz"), file("${base}.plasmid_read_ids.txt"), file("${base}.plasmid_extraction.bam")

script:
// if CLEAN_RIBOSOME_TRNA FLAG
if ( "${params.CLEAN_RIBOSOME_TRNA}" == true) {
    """
    #!/bin/bash
    #logging
    echo "ls of directory" 
    ls -lah 

    #cat ${minimap2_host_index} ${ribosome_trna} > host.fa

    minimap2 \
        -ax sr \
        -t "\$((${task.cpus}-2))" \
        -2 \
        ${ribosome_trna} \
        ${r1} ${r2}| samtools view -Sb -@ 2 - > ${base}.trna.bam

        samtools fastq -@ 4 -n -f 4 ${base}.trna.bam | pigz > ${base}.trna_filtered.fastq.gz
        samtools fastq -@ 4 -n -F 4 ${base}.trna.bam > ${base}.trna.mapped.bam

    minimap2 \
        -ax sr \
        -t "\$((${task.cpus}-2))" \
        -2 \
        ${minimap2_host_index} \
        ${base}.trna_filtered.fastq.gz | samtools view -Sb -@ 2 - > ${base}.host_mapped.bam
        samtools fastq -@ 4 -n -f 4 -1 ${base}.host_filtered_R1.fastq  -2 ${base}.host_filtered_R2.fastq ${base}.host_mapped.bam 
        pigz ${base}.host_filtered_R1.fastq 
        pigz ${base}.host_filtered_R2.fastq 

    minimap2 \
        -ax sr \
        -t ${task.cpus} \
        --sam-hit-only \
        ${minimap2_plasmid_db} \
        ${base}.host_filtered_R1.fastq ${base}.host_filtered_R2.fastq  | samtools view -F 4 -Sb - > ${base}.plasmid_extraction.bam

    samtools view ${base}.plasmid_extraction.bam | cut -f1 | sort | uniq > ${base}.plasmid_read_ids.txt

    /usr/local/miniconda/bin/seqkit grep -f ${base}.plasmid_read_ids.txt ${base}.host_filtered.fastq.gz | pigz > ${base}.plasmid.fastq.gz 

"""
    }
// If leaving tRNA in file
else {
    """
    #!/bin/bash
    #logging
    echo "ls of directory" 
    ls -lah 

    #cat ${minimap2_host_index} ${ribosome_trna} > host.fa


    minimap2 \
        -ax map-ont \
        -t "\$((${task.cpus}-2))" \
        -2 \
        ${minimap2_host_index} \
        ${r1} | samtools view -Sb -@ 2 - > ${base}.host_mapped.bam
        samtools fastq -@ 4 -n -f 4 ${base}.host_mapped.bam | pigz > ${base}.host_filtered.fastq.gz
    
    minimap2 \
        -ax map-ont \
        -t ${task.cpus} \
        --sam-hit-only \
        ${minimap2_plasmid_db} \
        ${base}.host_filtered.fastq.gz | samtools view -F 4 -Sb - > ${base}.plasmid_extraction.bam


    samtools view ${base}.plasmid_extraction.bam | cut -f1 | sort | uniq > ${base}.plasmid_read_ids.txt

    plasmid_count=`cat ${base}.plasmid_read_ids.txt | wc -l`
    echo "\$plasmid_count sequences mapped to plasmid" 

    /usr/local/miniconda/bin/seqkit grep -f ${base}.plasmid_read_ids.txt ${base}.host_filtered.fastq.gz | pigz > ${base}.plasmid.fastq.gz 
    /usr/local/miniconda/bin/seqkit grep -v -f ${base}.plasmid_read_ids.txt ${base}.host_filtered.fastq.gz | pigz > ${base}.host_filtered.plasmid_removed.fastq.gz 

    """
    }  
}

// run kraken2 prefiltering
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
    tuple val("${base}"), file("${base}.sorted.filtered.*.bam"), file("${base}.sorted.filtered.*.bam.bai")
    tuple val("${base}"), file ("${base}.*.unclassified_reads.txt")
script:
"""
#!/bin/bash

    #logging
    echo "ls of directory" 
    ls -lah 

    species_basename=`basename ${species_fasta} | cut -f1 -d .`

    # if this is the first attempt at running an alignment against this reference for this sample proceed

    if [ "${task.attempt}" -eq "1" ]
    then
        echo "running Minimap2 on ${base}"
        # run minimap2 and pipe to bam output 
        minimap2 \
            -ax sr \
            -t "\$((${task.cpus}-4))" \
            -2 \
            -K 25M \
            --split-prefix ${base}.split \
            ${species_fasta} \
            ${r1} ${r2} | samtools view -Sb -@ 4 - > ${base}.bam

        # extract mapped reads and sort 
        samtools view -Sb -F 4 ${base}.bam > ${base}.filtered.bam
        samtools sort -@ ${task.cpus} ${base}.filtered.bam -o ${base}.sorted.filtered.bam 

        # output unclassified reads
        samtools view -Sb -@  ${task.cpus} -f 4 ${base}.bam > ${base}.unclassified.bam

        # cleanup intermediate file
        # TODO uncomment later
        rm ${base}.bam

        # gather the read IDs of unassigned reads to extract from host filtered fastq downstream
        samtools view ${base}.unclassified.bam | cut -f1 > ${base}.\$species_basename.\$RANDOM.unclassified_reads.txt
        
        # adding random identifier to species bams to avoid filename collisions while merging later
        mv ${base}.sorted.filtered.bam ${base}.sorted.filtered.\$species_basename.\$RANDOM.bam

        #index merged bam 
        samtools index ${base}.sorted.filtered.*.bam

        # stats for reads mapped and unmapped
        readsmapped=`samtools view -c ${base}.filtered.bam`
        readsunmapped=`samtools view -c  ${base}.unclassified.bam`
        echo "reads in filtered bam"
        echo \$readsmapped

        echo "reads in unclassified bam"
        echo \$readsunmapped

        # removing unclassified bam to save space
        rm ${base}.unclassified.bam
        
        # for some reason if Minimap2 fails because it ran out of memory it doesn't exit the process
        # To check for failed Minimap2, mapped and unmapped reads will both be 0, in which case the process crashes and reattempts
        if [ "\$readsmapped" -eq "0" -a "\$readsunmapped" -eq "0" ]
        then
            echo "minimap2 ran out of memory but failed to crash for ${base} retrying with fasta split"
            exit 1
        fi
    
    # if process reattempts because it ran out of memory 
    else

        # split the fasta into chunks smaller chunks depending on how many times the process has been attempted 

        echo "running Minimap2 RNA on ${base} attempt ${task.attempt}"
        #echo "fasta being split \$split_num times"
        
        # faSplit has some weird splitting activity but it works
        /usr/local/miniconda/bin/faSplit sequence ${species_fasta} ${task.attempt} genus_split
        
        #NEED TO FIX: check within the loop for blank output. Minimap2 running out of memory might not crash the loop
        # something like if bam empty, exit 1
        for f in `ls genus_split*`
        do
            minimap2 \
                -ax sr \
                -t "\$((${task.cpus}-4))" \
                -2 \
                --split-prefix ${base}.split \
                \$f \
                ${r1} ${r2} | samtools view -Sb -@ 4 - > ${base}.\$f.bam
            samtools sort -@ ${task.cpus} ${base}.\$f.bam -o ${base}.sorted.temp.bam
            bamcount=`samtools view -c ${base}.sorted.temp.bam`
            # check if an individual split caused an out of memory error and exit the process 
            if [ "\$bamcount" -eq "0" ]
                then
                echo "minimap2 ran out of memory but failed to crash for ${base} retrying with fasta split"
                exit 1
            fi
            mv ${base}.sorted.temp.bam ${base}.sorted.\$RANDOM.bam
        done

        # merge the fasta split alignments 
        samtools merge ${base}.merged.bam ${base}.sorted.*.bam

        # extract mapped reads and sort the bam 
        samtools view -Sb -F 4 ${base}.merged.bam > ${base}.filtered.bam
        samtools sort -@ ${task.cpus} ${base}.filtered.bam -o ${base}.sorted.filtered.bam 

        # output unclassified reads
        samtools view -Sb -@  ${task.cpus} -f 4 ${base}.merged.bam > ${base}.unclassified.bam

        # cleanup intermediate file to save space
        rm ${base}.merged.bam

        ##samtools fastq -@ 4 ${base}.unclassified.bam | pigz > ${base}.unclassified.fastq.gz
        samtools view ${base}.unclassified.bam | cut -f1 > ${base}.\$species_basename.\$RANDOM.unclassified_reads.txt
        
        mv ${base}.sorted.filtered.bam ${base}.sorted.filtered.\$species_basename.\$RANDOM.bam
        samtools index ${base}.sorted.filtered.*.bam

        readsmapped=`samtools view -c ${base}.filtered.bam`
        readsunmapped=`samtools view -c  ${base}.unclassified.bam`
        echo "reads in filtered bam"
        echo \$readsmapped

        echo "reads in unclassified bam"
        echo \$readsunmapped
        
        # removing unclassified bam to save space
        rm ${base}.unclassified.bam
        
        # for some reason if Minimap2 fails because it ran out of memory it doesn't exit the process
        # To check for failed Minimap2, mapped and unmapped reads will both be 0, in which case the process crashes and reattempts
        if [ "\$readsmapped" -eq "0" -a "\$readsunmapped" -eq "0" ]
        then
            echo "minimap2 ran out of memory but failed to crash for ${base} retrying with fasta split"
            exit 1
        fi
    fi
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

//TODO: ADD ACCESSION DNE OUTPUT BACK IN- this accounts for taxa that don't exist in our taxdump file but exist in the nt database or vice versa

// run LCA algorithm
process Classify { 
publishDir "${params.OUTPUT}/Classification/${base}", mode: 'symlink', overwrite: true
container 'vpeddu/nanopore_metagenomics:latest'
beforeScript 'chmod o+rw .'
errorStrategy 'ignore'
cpus 8
input: 
    tuple val(base), file(bam), file(bamindex), file(unclassified_fastq) //,  file(plasmid_fastq), val(plasmid_count)
    file taxdump
    file classify_script
    file accessiontotaxid

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

# for whatever reason if we don't copy the taxdump file, the original gets modified which breaks every other classification process
cp taxdump/*.dmp .

# run LCA script
python3 ${classify_script} ${bam} ${base} 

# counting unassigned reads to add back into final report
#echo \$(zcat ${unclassified_fastq} | wc -l)/4 | bc >> ${base}.prekraken.tsv
linecount=\$(zcat ${unclassified_fastq} | wc -l)
fastqlinecount=\$(awk -v lc=\$linecount 'BEGIN {  print (lc/4) }')
echo -e "0\\t\$fastqlinecount" >> ${base}.prekraken.tsv

# add plasmid count back into results

echo \$fastqlinecount \$linecount unclassified reads 

"""
}

// run LCA for eggnog-mapper output 
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

// write pavian style report
process Write_report { 
publishDir "${params.OUTPUT}/", mode: 'copy', overwrite: true
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

// write pavian style report for orthologs
process Write_report_orthologs { 
publishDir "${params.OUTPUT}/ortholog_reports/", mode: 'copy', overwrite: true
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