//
// Check input samplesheet and get read channels
//

include { MINIMAP2_ALIGN } from '../../modules/local/minimap2_align'
include { PAF_TO_SAM     } from '../../modules/local/paf_to_sam'
include { IVAR_CONSENSUS } from '../../modules/local/ivar_consensus'
include { IRMA           } from '../../modules/local/irma'
include { CONDENSE       } from '../../modules/local/condense'
include { BAM_STATS      } from '../../modules/local/bam_stats'
include { MAPPED_FASTQ   } from '../../modules/local/get_mapped_fastq'
include { NEXTCLADE      } from '../../modules/local/nextclade'
include { EXPORT_ALN     } from '../../modules/local/export_aln'


workflow ASSEMBLE {
    take:
    ch_ref_list // channel: [meta, ref_id, ref_path, reads]

    main:

    ch_versions = Channel.empty()
    ch_bam      = Channel.empty()

    /* 
    =============================================================================================================================
        ASSEMBLY OPTION 1: IRMA
    =============================================================================================================================
    */
    
    if (params.cons_assembler == "irma"){
        // MODULE: Run IRMA
        IRMA (
            ch_ref_list,
            file("${baseDir}/assets/IRMA_MODULE/", checkIfExists: true)
        )
        ch_versions = ch_versions.mix(IRMA.out.versions)

        // Create consensus and bam channels
        IRMA.out.bam.set{ ch_bam }
        IRMA.out.consensus.set{ ch_consensus }
    }
    
    /* 
    =============================================================================================================================
        ASSEMBLY OPTION 2: IVAR (with MINIMAP2)
    =============================================================================================================================
    */

    if (params.cons_assembler == "ivar"){
        // MODULE: Run MINIMAP2 ALIGN
        MINIMAP2_ALIGN (
            ch_ref_list
        )
        ch_versions = ch_versions.mix(MINIMAP2_ALIGN.out.versions)

        // MODULE: Convert PAF to SAM/BAM (if MINIMAP2 outputs PAF format)
        // Note: If MINIMAP2_ALIGN is modified to output BAM directly with -a flag,
        // this conversion step can be removed
        PAF_TO_SAM (
            MINIMAP2_ALIGN.out.paf
        )
        ch_versions = ch_versions.mix(PAF_TO_SAM.out.versions)
        PAF_TO_SAM.out.bam.set{ ch_bam }

        // MODULE: Run Ivar
        IVAR_CONSENSUS (
            ch_bam
        )
        ch_versions = ch_versions.mix(IVAR_CONSENSUS.out.versions)
        IVAR_CONSENSUS.out.consensus.set{ch_consensus}
    }

    /* 
    =============================================================================================================================
        GET MAPPING STATS & EXPORT MAPPED READS TO FASTQ
    =============================================================================================================================
    */
    // Module: Get BAM stats
    BAM_STATS (
        ch_bam
    )
    ch_versions = ch_versions.mix(BAM_STATS.out.versions)

    // MODULE: Convert mapped reads to FASTQ
    MAPPED_FASTQ (
        BAM_STATS.out.read_list.join(ch_ref_list, by: [0,1])
    )
    ch_versions = ch_versions.mix(MAPPED_FASTQ.out.versions)

    /* 
    =============================================================================================================================
        CONDENSE DUPLICATE ASSEMBLIES
    =============================================================================================================================
    */
    CONDENSE(
        ch_consensus
            .groupTuple(by: 0)
            .join(BAM_STATS.out.coverage.groupTuple(by: 0), by: 0)
            .map{ meta, ref_id1, assemblies, ref_id2, stats -> [ meta, assemblies, stats ] }

    )
    ch_versions = ch_versions.mix(CONDENSE.out.versions)
    CONDENSE
        .out
        .assembly
        .transpose()
        .map{ meta, assembly -> [ meta, getRefName(assembly).replaceAll("${meta.id}_", ''), assembly ] }
        .set{ ch_consensus }

    // Publish alignments of condensed assemblies
    EXPORT_ALN (
        ch_consensus
          .map{ meta, ref_id, assembly -> [ meta, ref_id ] }
          .join( ch_ref_list.map{ meta, ref_id, ref_path, reads -> [ meta, ref_id, ref_path ] }, by: [ 0, 1 ] )
          .join( ch_bam, by: [ 0, 1 ] )
    )


    /* 
    =============================================================================================================================
        ASSEMBLY QC METRICS
    =============================================================================================================================
    */
    // MODULE: Run Nextclade
    NEXTCLADE (
        ch_consensus
            .combine(ch_ref_list.map{ meta, ref_id, ref_path, reads -> [ ref_path, ref_id ] }.unique(), by: 1)
    )
    ch_versions = ch_versions.mix(NEXTCLADE.out.versions)
    
    // Load Nextclade results
    NEXTCLADE
        .out
        .json
        .splitJson()
        .filter{ meta, refId, samLen, refLen, item -> item.key == 'results' }
        .map{ meta, refId, samLen, refLen, ncResults -> [ meta: meta, refId: refId ] + gatherMetrics( samLen, refLen, ncResults ) }
        .set{ ch_nc }
    ch_nc.map{ [ it.meta, it.refId, it.rowCsv ] }.set{ ch_nc_stats }

    emit:
    bam_stats = BAM_STATS.out.stats // channel: [ val(meta), val(ref), path(stats)) ]
    nextclade = ch_nc_stats         // channel: [ val(meta), val(ref), path(stats)) ]
    consensus = ch_consensus        // channel: [ val(meta), val(ref_id), path(assembly) ]
    versions  = ch_versions
}

// Remove file extension for references - allows reference names to have periods in them (GCA_7839185.1.fa.gz)
def getRefName(filename){
    def refname = filename
    def tokens  = filename.getName().tokenize('\\.')
    if(tokens.size() > 1){
        refname = filename.toString().endsWith('.gz') ? tokens[0..-3].join('.') : tokens[0..-2].join('.')
    }
    
    return refname
}

// Function for gathering assembly metrics
def gatherMetrics( samLen, refLen, data ){
    def metrics            = data.value[0] ? data.value[0] : null
    def coverage           = metrics ? metrics.coverage : null
    def status             = metrics ? metrics.qc.overallStatus : null
    def totalSubstitutions = metrics ? metrics.totalSubstitutions : null
    def totalDeletions     = metrics ? metrics.totalDeletions : null
    def totalNonACGTNs     = metrics ? metrics.totalNonACGTNs : null
    def totalInsertions    = metrics ? metrics.totalInsertions : null
    def totalMissing       = metrics ? metrics.totalMissing : null
    def terminiMissing     = 0

    if(metrics){
        metrics
            .missing
            .findAll{ it.range.begin.toInteger() <= 0 || it.range.end.toInteger() >= refLen.toInteger() }
            .each{ terminiMissing += it.range.end.toInteger() - it.range.begin.toInteger() }
    }
    def rowCsv = [ samLen, refLen, totalSubstitutions, totalDeletions, totalInsertions, totalMissing, terminiMissing, totalNonACGTNs, coverage, status ].join(',')
    return [ rowCsv: rowCsv ]
}
