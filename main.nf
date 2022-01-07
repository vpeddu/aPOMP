#!/usr/bin/env nextflow
nextflow.enable.dsl=2

def helpMessage() {
// log.info "foo" prints to STDOUT and nextflow.log
    log.info"""
    TE_Flow: Pipeline that combines standard RNA-seq quantification approaches with TE-aware
    steps to generate a comprehensive accounting of transcriptomic activity.
    Usage:
    nextflow run rreggiar/TE_Flow \\
        --INPUT_FOLDER test_fastq_paired_small/ 
        --ADAPTER_CHOICE intra 
        --PAIRED_END True 
        --OUTPUT testout 
        --SALMON_INDEX test/test_salmon_index 
        -with-docker Ubuntu:18.04
    """.stripIndent()
}

// show help message
params.help = false
// The params scope allows you to define parameters that will be accessible in the pipeline script
if (params.help){
    helpMessage()
    exit 0
}

//Nanopore mode on by default 
//params.NANOPORE = true
//Minimap2 -ax splice off by default 
params.MINIMAPSPLICE = false
//Eggnog off by default 
params.EGGNOG = false
//Flye assembly off be default
params.METAFLYE = false

include { Trimming_FastP } from './illumina_modules.nf'
include { Low_complexity_filtering } from './illumina_modules.nf'
include { Host_depletion } from './illumina_modules.nf'
include { Kraken_prefilter } from './illumina_modules.nf'
include { Extract_db } from './illumina_modules.nf'
include { Minimap2_illumina } from './illumina_modules.nf'
include { Sam_conversion } from './illumina_modules.nf'
include { Classify } from './illumina_modules.nf'
include { Write_report } from './illumina_modules.nf'

include { NanoFilt } from './nanopore_modules.nf'
include { NanoPlot } from './nanopore_modules.nf'
include { Host_depletion_nanopore } from './nanopore_modules.nf'
include { Host_depletion_extraction_nanopore } from './nanopore_modules.nf'
include { Minimap2_nanopore } from './nanopore_modules.nf'
include { MetaFlye } from './nanopore_modules.nf'
include { Kraken_prefilter_nanopore } from './nanopore_modules.nf'
include { Diamond_translated_alignment_unclassified } from './nanopore_modules.nf'
include { Cluster_unclassified_reads } from './nanopore_modules.nf'
include { Eggnog_mapper } from './nanopore_modules.nf'
include { Extract_true_novel } from './nanopore_modules.nf'
include { Classify_orthologs } from './illumina_modules.nf'
include { Write_report_orthologs } from './illumina_modules.nf'


Star_index_Ch = Channel
            .fromPath("${params.INDEX}/star_host/")

Kraken2_db = Channel
            .fromPath("${params.INDEX}/kraken2_db/")

NT_db = Channel
            .fromPath("${params.INDEX}/nt_fasta_index/")

Taxdump = Channel
            .fromPath("${params.INDEX}/taxdump/")

Krakenuniq_db = Channel
            .fromPath("${params.INDEX}/krakenuniq_db/")

Eggnog_db = Channel
            .fromPath("${params.INDEX}/eggnog_db/")

// Accession_to_taxid = Channel
//                     .fromPath("${params.INDEX}/accession2taxid/")

// Minimap2_host_index = Channel
                    // .fromPath("${params.INDEX}/minimap2_host/minimap2_hg38.mmi")

workflow{
    if ( params.NANOPORE ){
        input_read_Ch = Channel
            .fromPath("${params.INPUT_FOLDER}**.fastq.gz")
            .map { it -> [it.name.replace(".fastq.gz", ""), file(it)]}
        NanoFilt(
            input_read_Ch
        )
        NanoPlot (
            NanoFilt.out[0]
        )
        Host_depletion_nanopore( 
            NanoFilt.out[0],
            //Minimap2_host_index
            file("${params.INDEX}/minimap2_host/hg38.fa"),
            file("${params.INDEX}/ribosome_trna/all_trna.fa")
        )
        if( params.METAFLYE ) {
            MetaFlye(
                Host_depletion_nanopore.out[0]
            )
            Kraken_prefilter_nanopore(
                MetaFlye.out,
                Kraken2_db.collect()
            )
        }
        else{
            Kraken_prefilter_nanopore(
                Host_depletion_nanopore.out[0],
                Kraken2_db.collect()
            )
        }
            Extract_db(
                Kraken_prefilter_nanopore.out,
                NT_db.collect(),
                file("${baseDir}/bin/extract_seqs.py")
                )
            Minimap2_nanopore( 
                Host_depletion_nanopore.out[0].groupTuple(size:1).join(
                    Extract_db.out)
                )

            if (params.EGGNOG){
                Cluster_unclassified_reads(
                    Minimap2_nanopore.out[1],
                )
                MetaFlye(
                    Cluster_unclassified_reads.out
                )
                // need to fix this so it doesn't collide with pre-classification metaflye
                Eggnog_mapper(
                    MetaFlye.out,
                    Eggnog_db.collect()
                )
                Classify_orthologs(
                    Eggnog_mapper.out, 
                    Taxdump.collect(),
                    file("${baseDir}/bin/orthologs_to_pavian.py"),
                    file("${params.INDEX}/accession2taxid/")
                )
                Write_report_orthologs(
                    Classify_orthologs.out[0],
                    Krakenuniq_db.collect()
                )
            }
            // Diamond_translated_alignment_unclassified(
            //     Minimap2_nanopore.out[1],
            //     Diamond_protein_db
            // )
            // MetaFlye( 
            //     Minimap2_nanopore.out[1]
            // )
            // Extract_true_novel(
            //     MetaFlye.out
            // )
            Classify ( 
                // works but can clean up groupTuple later
                Minimap2_nanopore.out[0].groupTuple(size:1).join(
                Minimap2_nanopore.out[1]), 
                Taxdump.collect(),
                file("${baseDir}/bin/classify_reads.py"),
                file("${params.INDEX}/accession2taxid/")
                )
            Write_report(
                Classify.out[0],
                Krakenuniq_db.collect()
            )
            //}
    }
    else if (params.ILLUMINA) {
        // Workflow assumes reads are paired 
        // Only runs classification with NT. No protein level classification
        // input channel structure: {val(base), file(R1), file(R2)}
        input_read_Ch = Channel
            .fromFilePairs("${params.INPUT_FOLDER}**_R{1,2}*.fastq.gz")
            .map { it -> [it[0], it[1][0], it[1][1]] }
        if (params.SKIP_TRIMMING){
            Low_complexity_filtering(
                input_read_Ch
            )
        }
        else{
            Trimming_FastP(
                input_read_Ch
            )
            Low_complexity_filtering(
                Trimming_FastP.out[0],
                )
        }
        Host_depletion(
            Low_complexity_filtering.out[0],
            Star_index_Ch.collect()
            )
        Kraken_prefilter(
            Host_depletion.out[2],
            Kraken2_db.collect()
            )
        Extract_db(
            Kraken_prefilter.out,
            NT_db.collect(),
            file("${baseDir}/bin/extract_seqs.py")
            )
        Minimap2_illumina( 
            Host_depletion.out[2].groupTuple(size:1).join(
                Extract_db.out) 
            )
        Classify ( 
            // works but can clean up groupTuple later
            Minimap2_illumina.out[0].groupTuple(size:1).join(
            Minimap2_illumina.out[1]), 
            Taxdump.collect(),
            file("${baseDir}/bin/classify_reads.py"),
            file("${params.INDEX}/accession2taxid/")
            )
        Write_report(
            Classify.out[0],
            Krakenuniq_db.collect()
        )
    }
}