//
// Check input samplesheet and get read channels
//

include { BWA_MEM          } from '../../modules/local/bwa_mem'
include { IVAR_CONSENSUS   } from '../../modules/local/ivar_consensus'
include { IRMA             } from '../../modules/local/irma'
include { BAM_STATS        } from '../../modules/local/bam_stats'
include { MAPPED_FASTQ     } from '../../modules/local/get_mapped_fastq'
include { NEXTCLADE_CONFIG } from '../../modules/local/nexclade/config/main'
include { NEXTCLADE_RUN    } from '../../modules/local/nexclade/run/main'


workflow ASSEMBLE {
    take:
    ref_list // channel: [meta, ref_id, ref_path, reads]

    main:

    ch_versions = Channel.empty()
    ch_bam = Channel.empty()

    /* 
    =============================================================================================================================
        ASSEMBLY OPTION 1: IRMA
    =============================================================================================================================
    */
    
    if (params.consensus_assembler == "irma"){
        // MODULE: Run IRMA
        IRMA (
            ref_list.groupTuple(by:0).map{ meta, ref_ids, ref_paths, reads -> [ meta, ref_paths, reads.get(0) ]},
            file("${baseDir}/assets/IRMA_MODULE/", checkIfExists: true)
        )

        IRMA
            .out
            .bam
            .flatMap{ meta, files -> files.collect{ file -> [ meta, file ] }}
            .map{ meta, bam -> [meta, file(bam).getSimpleName(), bam] }
            .set{ ch_bam }
        IRMA
            .out
            .consensus
            .flatMap{ meta, files -> files.collect( file -> [ meta, file ] )}
            .map{ meta, consensus -> [meta, file(consensus).getSimpleName(), consensus] }
            .set{ ch_consensus }
    }

    
    /* 
    =============================================================================================================================
        ASSEMBLY OPTION 2: IVAR
    =============================================================================================================================
    */

    if (params.consensus_assembler == "ivar"){
        // MODULE: Run BWA MEM
        BWA_MEM (
            ref_list
        )
        ch_versions = ch_versions.mix(BWA_MEM.out.versions.first())
        BWA_MEM.out.bam.set{ ch_bam }

        // MODULE: Run Ivar
        IVAR_CONSENSUS (
            ch_bam
        )
        ch_versions = ch_versions.mix(IVAR_CONSENSUS.out.versions.first())
        IVAR_CONSENSUS.out.consensus.set{ch_consensus}
    }

    /* 
    =============================================================================================================================
        GET MAPPING STATS & SPLIT FASTQS
    =============================================================================================================================
    */

    // MODULE: Summarize Samtool stats
    BAM_STATS (
        ch_bam
    )

    // MODULE: Convert mapped reads to FASTQ
    MAPPED_FASTQ (
        BAM_STATS.out.read_list.join(ref_list, by: [0,1])
    )
    

    /* 
    =============================================================================================================================
        ASSEMBLY QC METRICS
    =============================================================================================================================
    */
    // MODULE: Configure Nextclade QC file
    NEXTCLADE_CONFIG (
        ref_list.map{ meta, ref_id, ref_path, reads -> [ ref_id, ref_path ] }.unique(),
        file("${baseDir}/assets/nextclade-template.json", checkIfExists: true)
    )
    // MODULE: Run Nextclade
    NEXTCLADE_RUN (
        ch_consensus
            .combine(NEXTCLADE_CONFIG.out.config.map{ ref_id, ref, config -> [ ref, ref_id, config ] }, by: 1)
            .map{ ref_id, meta, consensus, ref, config -> [ meta, ref_id, consensus, ref, config ] }
    )

    emit:
    samtoolstats2tbl = BAM_STATS.out.stats    // channel: [ val(meta), val(ref), path(stats)) ]
    nextclade        = NEXTCLADE_RUN.out.tsv // channel: [ val(meta), val(ref), path(stats)) ]
}
