
// set defaults
params.MINIMAPSPLICE = false 
params.NANOFILT_MAXLENGTH = 20000
params.NANOFILT_MINLENGTH = 100
params.MINIMAP2_RETRIES = 10 
params.NANOFILT_QUALITY = 15

process NanoFilt { 
publishDir "${params.OUTPUT}/Nanofilt/${base}", mode: 'symlink', overwrite: true
// need to change this to the nanopore metagenomics container
//TODO: change container to metagenomics container

//process will exit 1 if filtered file is empty
//can happen with runs that were basecalled with lower accuracy 
//execution of the other samples will continue
errorStrategy 'ignore'

container " quay.io/biocontainers/nanofilt:2.8.0--py_0"
beforeScript 'chmod o+rw .'
cpus 6
input: 
    tuple val(base), file(r1)
output: 
    tuple val(base), file("${base}.filtered.fastq.gz")
    file "*"


script:
if ( params.REALTIME ){
    """
    #!/bin/bash
    #logging
    echo "ls of directory" 
    ls -lah 
    echo "running Nanofilt on ${base}"


    cat *.fastq >  tmp.merged.fastq

    # nanofilt doesn't have gzip support so we have to pipe in from gunzip
    cat tmp.merged.fastq | NanoFilt -q ${params.NANOFILT_QUALITY} \
            --maxlength ${params.NANOFILT_MAXLENGTH} \
            --length ${params.NANOFILT_MINLENGTH} | gzip > ${base}.filtered.fastq.gz

    if [[ \$(gunzip -c ${base}.filtered.fastq.gz | head -c1 | wc -c) == "0" ]] 
        then
            echo "${base}.filtered.fastq.gz is empty"
            exit 1
        fi

    """
}
else {
    """
    #!/bin/bash
    #logging
    echo "ls of directory" 
    ls -lah 
    echo "running Nanofilt on ${base}"

    # nanofilt doesn't have gzip support so we have to pipe in from gunzip
    gunzip -c ${r1} | NanoFilt -q ${params.NANOFILT_QUALITY} \
            --maxlength ${params.NANOFILT_MAXLENGTH} \
            --length ${params.NANOFILT_MINLENGTH} | gzip > ${base}.filtered.fastq.gz

    if [[ \$(gunzip -c ${base}.filtered.fastq.gz | head -c1 | wc -c) == "0" ]] 
        then
            echo "${base}.filtered.fastq.gz is empty"
            exit 1
        fi

    """
}
}


process NanoFilt_RT { 

publishDir "${params.OUTPUT}/Nanofilt/${base}", mode: 'symlink', overwrite: true
// need to change this to the nanopore metagenomics container
//TODO: change container to metagenomics container

//process will exit 1 if filtered file is empty
//can happen with runs that were basecalled with lower accuracy 
//execution of the other samples will continue
errorStrategy 'ignore'

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


cat *.fastq >  tmp.merged.fastq

# nanofilt doesn't have gzip support so we have to pipe in from gunzip
cat tmp.merged.fastq | NanoFilt -q ${params.NANOFILT_QUALITY} \
        --maxlength ${params.NANOFILT_MAXLENGTH} \
        --length ${params.NANOFILT_MINLENGTH} | gzip > ${base}.filtered.fastq.gz

if [[ \$(gunzip -c ${base}.filtered.fastq.gz | head -c1 | wc -c) == "0" ]] 
    then
        echo "${base}.filtered.fastq.gz is empty"
        exit 1
    fi

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

# run nanoplot 
NanoPlot -t ${task.cpus} \
    -p ${base} \
    --fastq ${r1} \
    --title ${base} 
"""
}

process Host_depletion_nanopore { 
publishDir "${params.OUTPUT}/Host_filtered/${base}", mode: 'symlink', overwrite: true
container "vpeddu/nanopore_metagenomics:latest"
beforeScript 'chmod o+rw .'
cpus 8
input: 
    tuple val(base), file(r1)
    file minimap2_host_index
    file ribosome_trna
    file minimap2_plasmid_db
output: 
    tuple val("${base}"), file("${base}.host_filtered.fastq.gz")
    file "${base}.host_mapped.bam"
    file "${base}.trna.mapped.bam"
    tuple val("${base}"), file("${base}.plasmid.fastq.gz"), file("${base}.plasmid_read_ids.txt"), file("${base}.plasmid_extraction.bam")

script:
// if CLEAN_RIBOSOME_TRNA FLAG
if ( "${params.LEAVE_TRNA_IN}" == true) {
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
// if filtering out tRNA and other stuff (default)
else 
    {
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
        ${ribosome_trna} \
        ${r1} | samtools view -Sb -@ 2 - > ${base}.trna.bam

        samtools fastq -@ 4 -n -f 4 ${base}.trna.bam | pigz > ${base}.trna_filtered.fastq.gz
        samtools fastq -@ 4 -n -F 4 ${base}.trna.bam > ${base}.trna.mapped.bam

    minimap2 \
        -ax map-ont \
        -t "\$((${task.cpus}-2))" \
        -2 \
        ${minimap2_host_index} \
        ${base}.trna_filtered.fastq.gz | samtools view -Sb -@ 2 - > ${base}.host_mapped.bam
        samtools fastq -@ 4 -n -f 4 ${base}.host_mapped.bam | pigz > ${base}.host_filtered.fastq.gz

    minimap2 \
        -ax asm5 \
        -t ${task.cpus} \
        --sam-hit-only \
        ${minimap2_plasmid_db} \
        ${base}.host_filtered.fastq.gz | samtools view -F 4 -Sb - > ${base}.plasmid_extraction.bam

    samtools view ${base}.plasmid_extraction.bam | cut -f1 | sort | uniq > ${base}.plasmid_read_ids.txt

    /usr/local/miniconda/bin/seqkit grep -f ${base}.plasmid_read_ids.txt ${base}.host_filtered.fastq.gz | pigz > ${base}.plasmid.fastq.gz 

"""
    }
}

// idenitfy resistant plasmids
process Identify_resistant_plasmids { 
publishDir "${params.OUTPUT}/plasmid_identification/${base}", mode: 'symlink', overwrite: true
container "vpeddu/nanopore_metagenomics:latest"
beforeScript 'chmod a+rw .'
cpus 8

// Process will fail if no plasmids are assembled. 
// Errorstrategy set to ignore so it doesn't cause the pipeline to exit
errorStrategy 'ignore'
input: 
    tuple val(base), file(plasmid_fastq), val(plasmid_read_ids), file(plasmid_bam)
    file amrdb
output: 
    tuple val("${base}"), file("${base}.amrfinder.out.txt"), file("${base}.plasmid.flye/assembly.fasta")
script:
"""
#!/bin/bash
#logging
echo "ls of directory" 
ls -lah 

# assemble plasmids with flye
# meta and plasmid flags are used here to find plasmids from a metagenomics sample 
# need error handling for if nothing is assembled
/Flye/bin/flye --plasmids \
    --meta \
    -t ${task.cpus} \
    -o ${base}.plasmid.flye \
    --nano-hq ${plasmid_fastq}

if [[ -f "${base}.plasmid.flye/assembly.fasta" ]]; then 
        # run amrfinder on flye assembly
        /amrfinder/amrfinder \
            -n ${base}.plasmid.flye/assembly.fasta \
            --threads ${task.cpus} \
            -d ${amrdb}/2021-12-21.1/ \
            -o ${base}.amrfinder.out.txt

    else 
        echo "flye assembly failed. Running amrfinder on plasmid fastq converted to fasta"
        seqtk seq -a ${plasmid_fastq} > ${base}.plasmid.flye/assembly.fasta
        /amrfinder/amrfinder \
                -n ${base}.plasmid.flye/assembly.fasta \
                --threads ${task.cpus} \
                -d ${amrdb}/2021-12-21.1/ \
                -o ${base}.amrfinder.out.txt
fi
"""
}

process MetaFlye { 
publishDir "${params.OUTPUT}/MetaFlye/${base}", mode: 'symlink', overwrite: true
container "quay.io/biocontainers/flye:2.9--py27h6a42192_0"
beforeScript 'chmod o+rw .'
errorStrategy 'ignore'
cpus 16
input: 
    tuple val(base), file(unassigned_fastq)
output: 
    tuple val("${base}"), file("${base}.flye.fasta.gz")
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

if [[ -f ${base}.flye/assembly.fasta ]]
then
    echo "flye assembled reads"
    mv ${base}.flye/assembly.fasta ${base}.flye.fasta
else
    echo "flye did not assemble reads" 
    mv ${unassigned_fastq} ${base}.flye.fasta
fi

mv ${base}.flye/assembly.fasta ${base}.flye.fasta

gzip ${base}.flye.fasta
"""
}

// should update with user definable flags
process Low_complexity_filtering_nanopore { 
publishDir "${params.OUTPUT}/low_comnplexity_filter_nanopore/${base}", mode: 'symlink', overwrite: true
container "quay.io/biocontainers/bbmap:38.76--h516909a_0"
beforeScript 'chmod o+rw .'
cpus 6
input: 
    tuple val(base), file(r1)
output: 
    tuple val(base), file("${base}.lcf_filtered.fastq.gz")

script:
"""
#!/bin/bash
#logging
echo "ls of directory" 
ls -lah 

# run bbduk for low complexity filtering
bbduk.sh \
    in1=${r1} \
    out1=${base}.lcf_filtered.fastq.gz \
    --ignorebadquality \
    entropy=0.7 \
    qin=33 \
    entropywindow=50 \
    entropyk=4 
"""
}

// run kraken prefiltering 
process Kraken_prefilter_nanopore { 
publishDir "${params.OUTPUT}/Kraken_prefilter/${base}", mode: 'symlink', overwrite: true
//#TODO need to fix container
container "staphb/kraken2"
beforeScript 'chmod o+rw .'
cpus 8
input: 
    tuple val(base), file(fastq)
    file kraken2_db
output: 
    tuple val("${base}"), file("${base}.kraken2.report")
script:
"""
#!/bin/bash
#logging
echo "ls of directory" 
ls -lah 

# run kraken2
kraken2 --db ${kraken2_db} \
    --threads ${task.cpus} \
    --classified-out ${base}.kraken2.classified \
    --output ${base}.kraken2.output \
    --report ${base}.kraken2.report \
    --gzip-compressed \
    --unclassified-out ${base}.kraken2.unclassified \
    ${fastq} 

"""
}

process Sourmash_prefilter_nanopore { 
publishDir "${params.OUTPUT}/Sourmash_prefilter/${base}", mode: 'symlink', overwrite: true
//#TODO need to fix container
container "vpeddu/nanopore_metagenomics:latest"
beforeScript 'chmod o+rw .'
cpus 8
input: 
    tuple val(base), file(fastq)
    file sourmash_db
    file taxdump
    file taxonomy_parse_script
output: 
    tuple val("${base}"), file("${base}.sourmash_to_genus.txt")
script:
"""
#!/bin/bash
#logging
echo "ls of directory" 
ls -lah 

/usr/local/bin/sourmash sketch dna -p scaled=1000,k=31 ${fastq} --name-from-first

/usr/local/bin/sourmash lca summarize --db ${sourmash_db} --query ${fastq}.sig -o ${base}.sourmash_lca_summ.csv --threshold 2 

cat ${base}.sourmash_lca_summ.csv | cut -f1,7 -d , | sed  '/^\$/d' > ${base}.sourmash_lca_summ.genus.csv

python3.7 ${taxonomy_parse_script} ${base}.sourmash_lca_summ.csv ${base}

"""
}



// run minimap2 
process Minimap2_nanopore { 
//conda "${baseDir}/env/env.yml"
publishDir "${params.OUTPUT}/Minimap2/${base}", mode: 'symlink'
container "vpeddu/nanopore_metagenomics"
beforeScript 'chmod o+rw .'
cpus 28
errorStrategy 'retry'
maxRetries params.MINIMAP2_RETRIES
input: 
    tuple val(base), file(species_fasta), file(r1)
output: 
    tuple val("${base}"), file("${base}.sorted.filtered.*.bam"), file("${base}.sorted.filtered.*.bam.bai")
    tuple val("${base}"), file ("${base}.*.unclassified_reads.txt")

script:
    // Default is False
    if ( params.MINIMAPSPLICE ) {
    """
    #!/bin/bash

    #logging
    echo "ls of directory" 
    ls -lah 

    echo "running Minimap2 DNA on ${base}"
    #TODO: FILL IN MINIMAP2 COMMAND 
    minimap2 \
        -ax splice \
        -t "\$((${task.cpus}-2))" \
        -2 \
        --split-prefix ${base}.split \
        ${species_fasta} \
        ${r1} | samtools view -Sb -@ 2 - > ${base}.bam

    samtools view -Sb -F 4 -q 40 ${base}.bam > ${base}.filtered.bam
    samtools sort -@ ${task.cpus} ${base}.filtered.bam -o ${base}.sorted.filtered.bam 
    #samtools index ${base}.sorted.filtered.bam
    # output unclassified reads
    samtools view -Sb -@  ${task.cpus} -f 4 ${base}.bam > ${base}.unclassified.bam

    # cleanup intermediate file
   # rm ${base}.bam

    species_basename=`basename ${species_fasta} | cut -f1 -d .`

    ##samtools fastq -@ 4 ${base}.unclassified.bam | gzip > ${base}.unclassified.fastq.gz
    samtools view -Sb ${base}.unclassified.bam | cut -f1 > ${base}.\$species_basename.unclassified_reads.txt
    #echo "reads in filtered bam"
    #samtools view -c ${base}.filtered.bam

    #echo "reads in unclassified bam"
    #samtools view -c  ${base}.unclassified.bam
    """
        }

    else {
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
            -ax asm20 \
            -t "\$((${task.cpus}-4))" \
            -2 \
            -K 25M \
            --split-prefix ${base}.split \
            ${species_fasta} \
            ${r1} | samtools view -Sb -@ 4 - > ${base}.bam

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
                -ax asm20 \
                -t "\$((${task.cpus}-4))" \
                -2 \
                --split-prefix ${base}.split \
                \$f \
                ${r1} | samtools view -Sb -@ 4 - > ${base}.\$f.bam
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
}

process Extract_fungi { 
//conda "${baseDir}/env/env.yml"
container "vpeddu/nanopore_metagenomics"
beforeScript 'chmod o+rw .'
cpus 4
errorStrategy 'retry'
maxRetries params.MINIMAP2_RETRIES
input: 
    file fungi_list
    file fastadb
output: 
    file "all_fungi.fasta.gz"

script:
    """
    #!/bin/bash

for i in `cat ${fungi_list}`
    do
    echo adding \$i
    if [[ -f ${fastadb}/\$i.genus.fasta.gz ]]; then
        ##cat ${fastadb}/\$i.genus.fasta.gz >> species.fasta.gz
        cp ${fastadb}/\$i.genus.fasta.gz \$i.fungi.genus.fasta.gz
    fi
done

    cat *.fungi.genus.fasta.gz > all_fungi.fasta.gz

    """
}

process Align_fungi { 
//conda "${baseDir}/env/env.yml"
publishDir "${params.OUTPUT}/Minimap2/${base}", mode: 'symlink'
container "vpeddu/nanopore_metagenomics"
beforeScript 'chmod o+rw .'
cpus 28
errorStrategy 'retry'
maxRetries params.MINIMAP2_RETRIES
input: 
    tuple val(base),file(r1)
    file fungi_db
output: 
    tuple val("${base}"), file("${base}.sorted.filtered.*.bam"), file("${base}.sorted.filtered.*.bam.bai")
    tuple val("${base}"), file ("${base}.*.unclassified_reads.txt")

script:
    """
    #!/bin/bash

    #logging
    echo "ls of directory" 
    ls -lah 

    # if this is the first attempt at running an alignment against this reference for this sample proceed

    species_basename='fungi'

    if [ "${task.attempt}" -eq "1" ]
    then
        echo "running Minimap2 on ${base}"
        # run minimap2 and pipe to bam output 
        minimap2 \
            -ax asm20 \
            -t "\$((${task.cpus}-4))" \
            -2 \
            -K 25M \
            --split-prefix ${base}.split \
            ${fungi_db} \
            ${r1} | samtools view -Sb -@ 4 - > ${base}.bam

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
        /usr/local/miniconda/bin/faSplit sequence ${fungi_db} ${task.attempt} genus_split
        
        #NEED TO FIX: check within the loop for blank output. Minimap2 running out of memory might not crash the loop
        # something like if bam empty, exit 1
        for f in `ls genus_split*`
        do
            minimap2 \
                -ax asm20 \
                -t "\$((${task.cpus}-4))" \
                -2 \
                --split-prefix ${base}.split \
                \$f \
                ${r1} | samtools view -Sb -@ 4 - > ${base}.\$f.bam
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

        readsmapped=`samtools view -c ${base}.filtered.bam`ssh
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


// Collect minimap2 alignments from each sample and merge into one large bam
process Collect_alignment_results{ 
publishDir "${params.OUTPUT}/Minimap2/${base}", mode: 'symlink'
container "vpeddu/nanopore_metagenomics"
beforeScript 'chmod o+rw .'
cpus 16
input: 
    tuple val(base), file(filtered_bam), file(bam_index), file(plasmid_fastq), file(plasmid_read_ids), file(plasmid_bam)
output: 
    tuple val("${base}"), file("${base}.merged.sorted.bam"), file("${base}.merged.sorted.bam.bai")

script:
    """
    #!/bin/bash

    # merging plasmid bam in here so it goes into LCA algorithm
    samtools sort -@ ${task.cpus} ${plasmid_bam} -o ${base}.sorted.filtered.plasmid.bam
    samtools index ${base}.sorted.filtered.plasmid.bam


    #samtools merge ${base}.merged.filtered.bam *.sorted.filtered.*.bam
    #samtools sort -@ ${task.cpus} ${base}.merged.filtered.bam -o ${base}.merged.sorted.bam
    #samtools index ${base}.merged.sorted.bam 

    # speeding up samtools merge using gnu parallel  
    find . -name '*.sorted.filtered.*.bam' |
        parallel -j${task.cpus} --tmpdir . -N4 -m --files samtools merge -u - |
        parallel --xargs samtools merge -@2 ${base}.merged.filtered.bam {}";" rm {}

    samtools sort -@ ${task.cpus} ${base}.merged.filtered.bam -o ${base}.merged.sorted.bam
    samtools index ${base}.merged.sorted.bam 


    """
}

// Collect unaligned reads for each sample and extract the unique reads from the host filtered fastq
process Collect_unassigned_results{ 
//conda "${baseDir}/env/env.yml"
publishDir "${params.OUTPUT}/Minimap2/${base}", mode: 'symlink'
container "vpeddu/nanopore_metagenomics"
beforeScript 'chmod o+rw .'
cpus 4
input: 
    tuple val(base), file(unclassified_fastq), file(depleted_fastq), file(fungi_unclassified)
    file filter_unassigned_reads
    //tuple val(base), file(plasmid_fastq), file(plasmid_read_ids)

    
output: 
    tuple val("${base}"), file ("${base}.merged.unclassified.fastq.gz") //, file("${base}.plasmid_unclassified_intersection.fastq.gz") , env(plasmid_count)

script:
    """
    #!/bin/bash

    #cat *.unclassified_reads.txt | sort | uniq > unique_unclassified_read_ids.txt
    python3 ${filter_unassigned_reads}
    /usr/local/miniconda/bin/seqtk subseq ${depleted_fastq} true_unassigned_reads.txt | gzip > ${base}.merged.unclassified.fastq.gz


    """
}

// TODO: UPDATE INDEX SO WE CAN USE NEWEST VERSION OF DIAMOND
//legacy and not being used anymore. Could add back in later as an alternative to eggnog but it uses so much memory 
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

//Cluster unclassified reads using mmseqs cluster
// could use linclust 
process Cluster_unclassified_reads { 
publishDir "${params.OUTPUT}/Mmsesq2_unclassified_translated/${base}", mode: 'symlink', overwrite: true
container "quay.io/biocontainers/mmseqs2:13.45111--h95f258a_1"
beforeScript 'chmod o+rw .'
cpus 16
input: 
    tuple val(base), file(unassigned_fastq)
    //tuple val(base), file(unclassified_bam), file(unclassified_fastq)
output: 
    tuple val("${base}"), file("${base}.mmseq.clustered.fasta")
script:
"""
#!/bin/bash
#logging
echo "ls of directory" 
ls -lah 

#create mmseq db
mmseqs createdb ${unassigned_fastq} ${base}.mmseq.DB

# cluster with mmseq cluster2
mmseqs cluster --threads ${task.cpus} ${base}.mmseq.DB ${base}.mmseq.DB_clu tmp
#extract representative sequences and convert back to fasta
mmseqs createsubdb ${base}.mmseq.DB_clu ${base}.mmseq.DB ${base}.mmseq.clu_rep
mmseqs convert2fasta ${base}.mmseq.clu_rep ${base}.mmseq.clustered.fasta

"""
}

// Map orthologous groups using eggnog mapper 
process Eggnog_mapper { 
publishDir "${params.OUTPUT}/Eggnog/${base}", mode: 'symlink', overwrite: true
container "quay.io/biocontainers/eggnog-mapper:2.1.6--pyhdfd78af_0"
beforeScript 'chmod o+rw .'
cpus 24
input: 
    tuple val(base), file(unassigned_fastq)
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
        # convert unclassified fastq to fasta for eggnog
        gunzip -f ${unassigned_fastq} 
        #sed -n '1~4s/^@/>/p;2~4p' ${base}.unclassified.fastq > ${base}.unclassified.fasta
        emapper.py \
            --dmnd_frameshift 15 \
            --itype metagenome \
            -i ${base}.flye.fasta \
            --cpu ${task.cpus} \
            --data_dir ${eggnog_db} \
            --report_orthologs \
            --decorate_gff yes \
            -o ${base} 
    else
        touch ${base}.diamond.out.blankinput

fi
# removed 
#-m mmseqs 
#--excel 
"""
}
// legacy don't need this anymore
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
    ${unassigned_fastq} \
    ${metaflye_contigs}| samtools view -Sb -f 4 -@ 4 - > ${base}.unassembled.unclassified.bam

samtools fastq -@ 4 ${base}.unassembled.unclassified.bam | gzip > ${base}.unassembled.unclassified.fastq.gz

"""
}

process Merge_classifiy_RT { 
publishDir "${params.OUTPUT}/merge_classify_tmp/${base}", mode: 'symlink'
input: 
    file prekraken
output: 
    file outprekraken

container "vpeddu/nanopore_metagenomics"
beforeScript 'chmod o+rw .'
script: 
"""
#!/bin/bash

cat *.prekraken > combined_prekraken.txt

awk '{arr[\$1]+=\$2} END {for (i in arr) {print i,arr[i]}}' combined_prekraken.txt

"""
}


// write pavian style report
process Accumulate_reports { 
publishDir "${params.OUTPUT}/Accumulate/", mode: 'copy', overwrite: true
container "vpeddu/nanopore_metagenomics:latest"
beforeScript 'chmod o+rw .'
errorStrategy 'ignore'
cpus 1
input: 
    file prekraken
    //file krakenuniqdb
    //file mergescript
output: 
    file "accumulated.*.prekraken.tsv"
    //file krakenuniqdb
    //file "*.rt.report.tsv"
    //file mergescript

script:
"""
#!/bin/bash
#logging
echo "ls of directory" 
ls -lah 
prev=\$((${task.index}-1))
echo \$prev

echo "iteration ${task.index}"
cat ${prekraken} > accumulated.${task.index}.prekraken.tsv

"""
}

process Write_report_RT { 
publishDir "${params.OUTPUT}/RT_out/", mode: 'copy', overwrite: true
//container "evolbioinfo/krakenuniq:v0.5.8"
container "vpeddu/nanopore_metagenomics:latest"
beforeScript 'chmod o+rw .'
errorStrategy 'ignore'
cpus 1
input: 
    file prekraken
    file krakenuniqdb
    file mergescript
    //file krakenuniqdb
    //file mergescript
output: 
    file "*.rt.report.tsv"
    //file "*.rt.report.tsv"

script:
"""
#!/bin/bash
#logging
echo "ls of directory" 
ls -lah 

timestamp=\$( date +%T )
echo \$timestamp

python3 ${mergescript} ${prekraken} \$timestamp

/usr/local/miniconda/bin/krakenuniq-report --db ${krakenuniqdb} \
--taxon-counts \$timestamp.merged.prekraken.tsv > \$timestamp.rt.report.tsv
#mv \$timestamp.merged.prekraken.tsv \$timestamp.rt.report.tsv
"""
}
process Combine_fq {
//publishDir "${params.OUTPUT}/", mode: 'copy', overwrite: true
container "vpeddu/nanopore_metagenomics:latest"
beforeScript 'chmod o+rw .'
errorStrategy 'ignore'
cpus 1
input:
//    tuple val(base), val(fq)
      path 'fq'
output:
    tuple env(basename), file("*.combined_fq.fastq")

script:
"""
#!/bin/bash
#logging
echo "ls of directory"
ls -lah

basename="\$RANDOM"

echo "new filename is \$basename.combined_fq.fastq"
cat fq* > \$basename.combined_fq.fastq
"""
}