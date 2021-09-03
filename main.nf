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
// accessible via nextflow run $project --help
params.help = false
// The params scope allows you to define parameters that will be accessible in the pipeline script
// my guess is params.help is a global boolean that depends on --help at CLI
if (params.help){
    // if help == T
    helpMessage()
    // run helpMessage()
    exit 0
    // clean exit
}

include { Trimming_FastP } from './modules.nf'
include { Host_depletion } from './modules.nf'
include { Interleave_FASTQ } from './modules.nf'
include { Metalign_db_selection } from './modules.nf'
include { Minimap2 } from './modules.nf'
include { Metalign_profiling } from './modules.nf'



Star_index_Ch = Channel
            .fromPath(params.HOST_STAR_INDEX)

Kraken2_db = Channel
            .fromPath(params.KRAKEN2_DB)

NT_db = Channel
            .fromPath(params.NT_DB)

Taxdump = Channel
            .fromPath(params.TAXDUMP)

Krakenuniq_db = Channel
            .fromPath(params.KRAKENUNIQUE_DB)

input_read_Ch = Channel
    .fromFilePairs("${params.INPUT_FOLDER}**_R{1,2}*.fastq.gz")
    .map { it -> [it[0], it[1][0], it[1][1]]}
workflow{
    // defined at CLI    
    Trimming_FastP(
        input_read_Ch
        )
    Host_depletion(
        // access output of preceeding process
        Trimming_FastP.out[0],
        // collects all items emitted by a channel to a list, return
        Star_index_Ch.collect()
        )
    Kraken_prefilter(
        Host_depletion.out[2],
        KRAKEN2_DB.collect()
        )
    Extract_db(
        Kraken_prefilter.out,
        NT_db.collect(),
        file("${basedir}/bin/extract_seqs.py")
        )
    Minimap2( 
        Host_depletion.out[2],
        Extract_db.out
        )
    Sam_conversion (
        Minimap2.out
        )
    Classify ( 
        Sam_conversion.out, 
        Taxdump.collect(),
        file("${basedir}/bin/classify_reads.py")
        )
    Write_report(
        Classify.out,
        Krakenuniq_db.collect()
    )
    }