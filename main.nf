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
params.NANOPORE = false

include { Trimming_FastP } from './illumina_modules.nf'
include { Low_complexity_filtering } from './illumina_modules.nf'
include { Host_depletion } from './illumina_modules.nf'
include { Kraken_prefilter } from './illumina_modules.nf'
include { Extract_db } from './illumina_modules.nf'
include { Minimap2 } from './illumina_modules.nf'
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



workflow{
    if ( params.NANOPORE){
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
            params.MINIMAP2_HOST_INDEX
        )
        Host_depletion_extraction_nanopore( 
            Host_depletion_nanopore.out,
        )
        // if metaflye specified
        if (params.METAFLYE){
            MetaFlye(
                Host_depletion_extraction_nanopore.out
            )
            Kraken_prefilter_nanopore(
                MetaFlye.out,
                Kraken2_db.collect()
            )
            Extract_db(
                Kraken_prefilter_nanopore.out,
                NT_db.collect(),
                file("${baseDir}/bin/extract_seqs.py")
                )
            Minimap2_nanopore( 
                Host_depletion_extraction_nanopore.out.groupTuple(size:1).join(
                    Extract_db.out)
                )
            Sam_conversion (
                Minimap2_nanopore.out
                )
            Classify ( 
                Sam_conversion.out[0], 
                Taxdump.collect(),
                file("${baseDir}/bin/classify_reads.py"),
                file("${params.ACCESSIONTOTAXID}")
                )
            Write_report(
                Classify.out[0],
                Krakenuniq_db.collect()
            )
            }
        // metaflye not specified
        else {
            Kraken_prefilter_nanopore(
                Host_depletion_extraction_nanopore.out,
                Kraken2_db.collect()
            )
                Kraken_prefilter_nanopore.out[1]
                    .splitCsv()
                    .combine(Kraken_prefilter_nanopore.out[0].splitCsv()).toList()
                    .join(
                    Host_depletion_extraction_nanopore.out).view()
            Minimap2_nanopore( 
                Kraken_prefilter_nanopore.out[1]
                    .splitCsv()
                    .combine(Kraken_prefilter_nanopore.out[0].splitCsv())
                    .join(
                    Host_depletion_extraction_nanopore.out),
                NT_db.collect()
                )
            Sam_conversion (
                Minimap2_nanopore.out.groupTuple()
                )
            Classify ( 
                Sam_conversion.out[0], 
                Taxdump.collect(),
                file("${baseDir}/bin/classify_reads.py"),
                file("${params.ACCESSIONTOTAXID}")
                )
            Write_report(
                Classify.out[0],
                Krakenuniq_db.collect()
            )
            }
    }
    else {
        input_read_Ch = Channel
            .fromFilePairs("${params.INPUT_FOLDER}**_R{1,2}*.fastq.gz")
            .map { it -> [it[0], it[1][0], it[1][1]] }
        // defined at CLI    
        Trimming_FastP(
            input_read_Ch
            )
        Low_complexity_filtering(
            Trimming_FastP.out[0],
        )
        Host_depletion(
            // access output of preceeding process
            Low_complexity_filtering.out[0],
            // collects all items emitted by a channel to a list, return
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
        Minimap2( 
            Host_depletion.out[2].groupTuple(size:1).join(
                Extract_db.out) 
            )
        Sam_conversion (
            Minimap2.out
            )
        Classify ( 
            Sam_conversion.out[0], 
            Taxdump.collect(),
            file("${baseDir}/bin/classify_reads.py"),
            file("${params.ACCESSIONTOTAXID}")
            )
        Write_report(
            Classify.out[0],
            Krakenuniq_db.collect()
        )
    }
    }