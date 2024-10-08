#!/usr/bin/env nextflow
nextflow.enable.dsl=2
nextflow.preview.recursion=true
def helpMessage() {
    log.info"""
aPOMP: A portable, open-source metagenomics pipeline.

Usage:
    --NANOPORE/--ILLUMINA sequencing platform used
    --INPUT_FOLDER folder containing input FASTQs (illumina must be paired)  
    --INDEX path to aPOMP index folder
    --OUTPUT output folder
    --IDENTIFY_RESISTANCE_PLASMIDS assemble and identify resistant plasmids 
    --CLEAN_RIBOSOME_TRNA filter out ribosomal and tRNA sequences before classification 
    --EGGNOG identify orthologous groups in unclassified reads using  the Eggnog-mapper
    --MINIMAPSPLICE run Minimap2 with preset -ax splice (default -ax map-ont)
    --NANOFILT_QUALITY [int] minimum quality threshold for nanofilt (default 10) 
    --KRAKEN2_THRESHOLD [int] discard Kraken2 results containing less than this many reads at or below the genus level (default 10)
    --LOW_COMPLEXITY_FILTER_NANOPORE run low complexity filtering on nanopore samples
quick run command: 
  nextflow run vpeddu/ev-meta \
		 --NANOPORE \
		 --IDENTIFY_RESISTANCE_PLASMIDS \
		 -profile standard \
		 --EGGNOG \
		 --INPUT_FOLDER <input_fastq_folder> \
		 --OUTPUT Zymo-GridION-EVEN-BB-SN_out \
		 --INDEX <index_path> \
		 --CLEAN_RIBOSOME_TRNA \
		 -with-docker 'ubuntu:18.04' \
		 -with-tower \
		 -with-report \
		 -latest \
		 -resume
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
//Eggnog off by default 
params.EGGNOG = false
//Flye assembly off be default
params.METAFLYE = false
params.NANOPORE = true
params.ALIGN_ALL_FUNGI = false
params.LEAVE_TRNA_IN = false
params.REALTIME = false
params.PREFILTER_THRESHOLD = 10
params.LOW_COMPLEXITY_FILTER_NANOPORE = false
params.IDENTIFY_RESISTANCE_PLASMIDS = false

// Import modules from modules files
include { Trimming_FastP } from './illumina_modules.nf'
include { Low_complexity_filtering } from './illumina_modules.nf'
include { Host_depletion_illumina } from './illumina_modules.nf'
include { Kraken_prefilter } from './illumina_modules.nf'
include { Extract_db } from './illumina_modules.nf'
include { Minimap2_illumina } from './illumina_modules.nf'
include { Sam_conversion } from './illumina_modules.nf'
include { Classify } from './illumina_modules.nf'
include { Classify_RT } from './illumina_modules.nf'
include { Write_report } from './illumina_modules.nf'
include { Write_report_RT } from './nanopore_modules.nf'
include { Accumulate_reports } from './nanopore_modules.nf'
include { Sourmash_prefilter_nanopore } from './nanopore_modules.nf'

include { Combine_fq } from './nanopore_modules.nf'
include { Chopper } from './nanopore_modules.nf'
include { Chopper_RT } from './nanopore_modules.nf'
include { NanoPlot } from './nanopore_modules.nf'
include { Low_complexity_filtering_nanopore } from './nanopore_modules.nf'
include { Host_depletion_nanopore } from './nanopore_modules.nf'
include { Minimap2_nanopore } from './nanopore_modules.nf'
include { Identify_resistant_plasmids } from './nanopore_modules.nf'
include { Collect_alignment_results } from './nanopore_modules.nf'
include { Collect_unassigned_results } from './nanopore_modules.nf'
include { Collect_alignment_results_RT } from './nanopore_modules.nf'
include { Collect_unassigned_results_RT } from './nanopore_modules.nf'
include { Collect_unassigned_results_illumina } from './illumina_modules.nf'
include { MetaFlye } from './nanopore_modules.nf'
include { Extract_fungi } from './nanopore_modules.nf'
include { Align_fungi } from './nanopore_modules.nf'
include { Kraken_prefilter_nanopore } from './nanopore_modules.nf'
include { Diamond_translated_alignment_unclassified } from './nanopore_modules.nf'
include { Cluster_unclassified_reads } from './nanopore_modules.nf'
include { Eggnog_mapper } from './nanopore_modules.nf'
include { Extract_true_novel } from './nanopore_modules.nf'
include { Classify_orthologs } from './illumina_modules.nf'
include { Write_report_orthologs } from './illumina_modules.nf'


// Define input channels 
Star_index_Ch = Channel
            .fromPath("${params.INDEX}/star_host/bowtie2_index/")

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

Amrfinder_db = Channel
            .fromPath("${params.INDEX}/plasmid_db/amrfinder/")



// Workflow logic
workflow{
    // Nanopore workflow (default)
    if ( params.NANOPORE ){
        // Create input tuple for each fastq in input folder
        //tuple is val(base),file(fastq)
        //basename is anything before ".fastq.gz"

        if ( params.REALTIME ) {
            params.PREFILTER_THRESHOLD = 1
            params.CHOPPER_QUALITY = 10
            input_read_Ch = Channel.watchPath("${params.FAST5_FOLDER}*.fastq")
            .map { it -> file(it) }
//            .map { it -> [it.name.replace(".fastq", ""), file(it)]}
            .buffer( size: 4, remainder: true)
            // input_read_Ch.view()
            Combine_fq(
                input_read_Ch
            )
            Chopper(
                Combine_fq.out[0]
                )
        } else { input_read_Ch = Channel
            .fromPath("${params.INPUT_FOLDER}**.fastq.gz")
            .map { it -> [it.name.replace(".fastq.gz", ""), file(it)]}
        //run nanofilt
            Chopper(
            input_read_Ch
                )
        }

        //run nanoplot on nanofilt output
        NanoPlot (
            Chopper.out[0]
        )
        // if --LOW_COMPLEXITY_FILTER_NANOPORE, run low complexity filtering on nanofilt  output 
        if( params.LOW_COMPLEXITY_FILTER_NANOPORE){
            Low_complexity_filtering_nanopore(
                Chopper.out[0]
            )
            // run host depletion after low complexity filtering
            Host_depletion_nanopore(
                Low_complexity_filtering_nanopore.out[0],
                file("${params.INDEX}/minimap2_host/hg38.fa"),
                file("${params.INDEX}/ribosome_trna/all_trna.fa"),
                file("${params.INDEX}/plasmid_db/plsdb.mmi")
            )
        }
        // if --LOW_COMPLEXITY_FILTER_NANOPORE not specified, run host depletion on nanofilt output 
        else {
            // tRNA depletion is evaluated within the process
            Host_depletion_nanopore( 
                Chopper.out[0],
                //Minimap2_host_index
                file("${params.INDEX}/minimap2_host/hg38.fa"),
                file("${params.INDEX}/ribosome_trna/all_trna.fa"),
                file("${params.INDEX}/plasmid_db/plsdb.mmi")
        )
        }
        // Identify resistant plasmids if --IDENTIFY_RESISTANCE_PLASMIDS specified
        if ( params.IDENTIFY_RESISTANCE_PLASMIDS ){ 
            Identify_resistant_plasmids(
                Host_depletion_nanopore.out[3],
                Amrfinder_db.collect()
            )
        }
        // Assemble reads before classification
        // very resource intensive 
        if( params.METAFLYE ) {
            MetaFlye(
                Host_depletion_nanopore.out[0]
            )
            Kraken_prefilter_nanopore(
                Host_depletion_nanopore.out[0],
                Kraken2_db.collect()
            )
        }
        // Run kraken2 prefiltering 
        else{
            if ( params.KRAKEN_PREFILTER ) { 
            Kraken_prefilter_nanopore(
                Host_depletion_nanopore.out[0],
                Kraken2_db.collect()
            )
            prefilterCh =
                Kraken_prefilter_nanopore.out[0]
                .filter{ it -> it[1].countLines() > 0 }

            prefilterfailCh= Kraken_prefilter_nanopore.out[0]
                .filter{ it ->  it[1].countLines() == 0}
                .map{ it -> it[0]}

                prefilterfailCh.view(f -> "\u001B[33m" + "$f failed prefiltering- skipping" + "\u001B[0m")

            } else { 
            Sourmash_prefilter_nanopore( 
                Host_depletion_nanopore.out[0],
                file("${params.INDEX}/sourmash/sourmash.nt.k31.lca.json"),
                file("${params.INDEX}/taxdump/taxa.sqlite"),
                file("${baseDir}/bin/sourmash_to_taxonomy.py"),
                params.PREFILTER_THRESHOLD
            )
            prefilterCh =
                Sourmash_prefilter_nanopore.out[0]
                .filter{ it -> it[1].countLines() > 0 }

            prefilterfailCh= Sourmash_prefilter_nanopore.out[0]
                .filter{ it ->  it[1].countLines() == 0}
                .map{ it -> it[0]}

                prefilterfailCh.view(f -> "\u001B[33m" + "$f failed prefiltering- skipping" + "\u001B[0m")

            }
        }
            // extract databases from NT index 
            if ( params.KRAKEN_PREFILTER ) { 
            Extract_db(
                prefilterCh,
                NT_db.collect(),
                file("${baseDir}/bin/extract_seqs.py"),
                file("${baseDir}/bin/fungi_genera_list.txt"),
                params.PREFILTER_THRESHOLD
                ) 
            } else { 
            Extract_db(
                prefilterCh,
                NT_db.collect(),
                file("${baseDir}/bin/extract_seqs.py"),
                file("${baseDir}/bin/fungi_genera_list.txt"),
                params.PREFILTER_THRESHOLD
                ) 
            }
            // run Minimap2 on each individual genus 
            // Extract_db.out.flatten().map{
            //     it -> [it.name.split("__")[0], it]}.map{
            //         it -> [it[0],it[1],it[1].name.split("--")[1].split('.fasta.gz')[0]]}.
            //         combine(Host_depletion_nanopore.out[0], by:0).view()
            // Extract_db.out.flatten().map{
            //         it -> [it.name.split("__")[0], it]}.view()
                    
                    //.flatMap{
                    //    it -> [it[0], it[1], it.size()]}.view()

            Minimap2_nanopore( 
            Extract_db.out.flatten().map{
                it -> [it.name.split("__")[0], it]}.map{
                    it -> [it[0],it[1],it[1].name.split("--")[1].split('.fasta.gz')[0]]}.
                    combine(Host_depletion_nanopore.out[0], by:0)

                // Extract_db.out.flatten().map{
                //     it -> [it.name.split("__")[0], it]}.combine(Host_depletion_nanopore.out[0], by:0)
                )
            // Collection hold for each sample's Minimap2 aligned results
            // BAMs are merged at this step for each sample 

            // if ( params.ALIGN_ALL_FUNGI ) { 
            //     Extract_fungi(
            //         file("${baseDir}/bin/fungi_genera_list.txt"),
            //         NT_db.collect()
            //     )
            //     // Host_depletion_nanopore.out[0].view()
            //     // prefilterCh.view()
            //     Minimap2_nanopore( 
            //         Host_depletion_nanopore.out[0].join(Extract_fungi.out[0]).map{
            //             it -> [it[0], it[2], it[1]]}
            //     )
            //     // fungiCh = Channel( 
            //     //     Minimap2_nanopore.out[0].mix(Align_fungi.out[0])
            //     // ) 
            //     /*
            //     Collect_alignment_results(
            //     Minimap2_nanopore.out[0].mix(Align_fungi.out[0])groupTuple().join(
            //     Host_depletion_nanopore.out[3]
            //     ) 
            // )*/
            // // Collection hold for each sample's Minimap2 unaligned results
            // // Unique read IDs found to be unassignable are extracted from the host filtered fastq here for downstream classification
            //  /* Collect_unassigned_results(
            //     Minimap2_nanopore.out[1].mix(Align_fungi.out[1]).groupTuple().join(
            //     Host_depletion_nanopore.out[0]
            //     ),
            //     file("${baseDir}/bin/filter_unassigned_reads.py")
            // ) */

            // Collect_alignment_results(
            //     Minimap2_nanopore.out[0].map{ key, bam, bai -> tuple( groupKey(key, key.size()), bam, bai ) }.join(
            //     Host_depletion_nanopore.out[3]

            //     )
            // )
            // // Collection hold for each sample's Minimap2 unaligned results
            // // Unique read IDs found to be unassignable are extracted from the host filtered fastq here for downstream classification
            // Collect_unassigned_results(
            //     Minimap2_nanopore.out[1].map{ key, unclassified -> tuple( groupKey(key, key.size()), unclassified ) }.join(
            //     Host_depletion_nanopore.out[0]),
            //     file("${baseDir}/bin/filter_unassigned_reads.py")
            // )


            //}
            //else{

        if ( params.REALTIME ) {

            Minimap2_nanopore.out[0].map{ key, bam, bai, aln_num -> tuple( groupKey(key, aln_num.toInteger()), bam, bai, aln_num ) }.groupTuple().join(
                Host_depletion_nanopore.out[3]
                ).view()

            Collect_alignment_results_RT(
                Minimap2_nanopore.out[0].map{ key, bam, bai, aln_num -> tuple( groupKey(key, aln_num.toInteger()), bam, bai, aln_num ) }.groupTuple().join(
                Host_depletion_nanopore.out[3]
                )
            )
                // Collect_alignment_results_RT(
                // Minimap2_nanopore.out[0].groupTuple().join(
                // Host_depletion_nanopore.out[3]

                // //Minimap2_nanopore.out[0].groupTuple().join(
                // //Host_depletion_nanopore.out[3]
                // )
            //)
            // Collection hold for each sample's Minimap2 unaligned results
            // Unique read IDs found to be unassignable are extracted from the host filtered fastq here for downstream classification
            Collect_unassigned_results_RT(
                Minimap2_nanopore.out[1].map{ key, unclassified_file, aln_num -> tuple( groupKey(key, aln_num.toInteger()), unclassified_file, aln_num ) }.groupTuple().join(
                Host_depletion_nanopore.out[0]),
                // ).join( Align_fungi.out[1].map{ key, unclassified -> tuple( groupKey(key, key.size()), unclassified ) }
                // ),
                file("${baseDir}/bin/filter_unassigned_reads.py")
            )//}
            Classify_RT ( 
                Collect_alignment_results_RT.out.join(
                Collect_unassigned_results_RT.out), 
                Taxdump.collect(),
                file("${baseDir}/bin/classify_reads.py")
                        )
        } else {

            Collect_alignment_results(
                Minimap2_nanopore.out[0].map{ key, bam, bai, aln_num -> tuple( groupKey(key, aln_num.toInteger()), bam, bai, aln_num ) }.groupTuple().join(
                Host_depletion_nanopore.out[3]

                //Minimap2_nanopore.out[0].groupTuple().join(
                //Host_depletion_nanopore.out[3]
                )
            )
            // Collection hold for each sample's Minimap2 unaligned results
            // Unique read IDs found to be unassignable are extracted from the host filtered fastq here for downstream classification
            Collect_unassigned_results(
                Minimap2_nanopore.out[1].map{ key, unclassified_file, aln_num -> tuple( groupKey(key, aln_num.toInteger()), unclassified_file, aln_num ) }.groupTuple().join(
                Host_depletion_nanopore.out[0]),
                // ).join( Align_fungi.out[1].map{ key, unclassified -> tuple( groupKey(key, key.size()), unclassified ) }
                // ),
                file("${baseDir}/bin/filter_unassigned_reads.py")
            )//}
            // run LCA for each sample 
            Classify ( 
                Collect_alignment_results.out.join(
                Collect_unassigned_results.out).groupTuple(), 
                Taxdump.collect(),
                file("${baseDir}/bin/classify_reads.py")
                    )
        }
            // if --EGGNOG run clustering, metaflye, and the eggnog OG search
            if (params.EGGNOG){
                Cluster_unclassified_reads(
                    Collect_unassigned_results.out,
                )
                MetaFlye(
                    Cluster_unclassified_reads.out
                )
                Eggnog_mapper(
                    MetaFlye.out,
                    Eggnog_db.collect()
                )
                // Eggnog provides taxIDs with each OG. Here the taxids are extracted and run through LCA 
                Classify_orthologs(
                    Eggnog_mapper.out, 
                    Taxdump.collect(),
                    file("${baseDir}/bin/orthologs_to_pavian.py")
                )
                // write OG classification results to pavian file
                Write_report_orthologs(
                    Classify_orthologs.out[0],
                    Krakenuniq_db.collect()
                )
            }
            // write pavian report for each sample 
            if ( params.REALTIME ) { 
            Accumulate_reports.scan ( 
                Classify_RT.out[1] 
                ) 
            Write_report_RT( 
                Accumulate_reports.out,
                file("${params.INDEX}/taxdump/taxa.sqlite"),
                file("${baseDir}/bin/merge_realtime_prekraken.py"),
                file("${baseDir}/bin/write_report.py"),
                )
                //Krakenuniq_db.collect() //,
                //file("${baseDir}/bin/merge_rt_reports.py")
                
                // Write_report_RT( 
                //     Channel.watchPath("${params.OUTPUT}/Classification/**.prekraken.tsv"),
                //     Krakenuniq_db.collect()
                // )
            } else { 
            Write_report(
                Classify.out[0],
                file("${params.INDEX}/taxdump/taxa.sqlite"),
                file("${baseDir}/bin/write_report.py")
            )
            }
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
        Host_depletion_illumina(
            Low_complexity_filtering.out[0],
            Star_index_Ch.collect(),
            file("${params.INDEX}/ribosome_trna/all_trna.fa"),
            file("${params.INDEX}/plasmid_db/plsdb.mmi")
            )
        Kraken_prefilter(
            Host_depletion_illumina.out[0],
            Kraken2_db.collect()
            )
        Extract_db(
            Kraken_prefilter.out,
            NT_db.collect(),
            file("${baseDir}/bin/extract_seqs.py"),
            file("${baseDir}/bin/fungi_genera_list.txt")
            )
        Minimap2_illumina( 
            Extract_db.out.flatten().map{
                it -> [it.name.split("__")[0], it]}.combine(Host_depletion_illumina.out[0], by:0)
            )
        Collect_alignment_results(
            Minimap2_illumina.out[0].groupTuple().join(
            Host_depletion_illumina.out[3]
            )
        )
            // Collection hold for each sample's Minimap2 unaligned results
            // Unique read IDs found to be unassignable are extracted from the host filtered fastq here for downstream classification
        Collect_unassigned_results_illumina(
            Minimap2_illumina.out[1].groupTuple().join(
            Host_depletion_illumina.out[0]
            ),
            file("${baseDir}/bin/filter_unassigned_reads.py")
            )
            // run LCA for each sample 
        Classify ( 
            Collect_alignment_results.out.join(
            Collect_unassigned_results_illumina.out).groupTuple(), 
            Taxdump.collect(),
            file("${baseDir}/bin/classify_reads.py")
                )
        Write_report(
            Classify.out[0],
            file("${params.INDEX}/taxdump/taxa.sqlite"),
            file("${baseDir}/bin/write_report.py")
        )
    }
}