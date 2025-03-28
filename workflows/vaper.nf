/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PRINT PARAMS SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { paramsSummaryLog; paramsSummaryMap } from 'plugin/nf-validation'

def logo = NfcoreTemplate.logo(workflow, params.monochrome_logs)
def citation = '\n' + WorkflowMain.citation(workflow) + '\n'
def summary_params = paramsSummaryMap(workflow)

// Print parameter summary log to screen
log.info logo + paramsSummaryLog(workflow) + citation

WorkflowVaper.initialise(params, log)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONFIG FILES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

ch_multiqc_config          = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
ch_multiqc_custom_config   = params.multiqc_config ? Channel.fromPath( params.multiqc_config, checkIfExists: true ) : Channel.empty()
ch_multiqc_logo            = params.multiqc_logo   ? Channel.fromPath( params.multiqc_logo, checkIfExists: true ) : Channel.empty()
ch_multiqc_custom_methods_description = params.multiqc_methods_description ? file(params.multiqc_methods_description, checkIfExists: true) : file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Local modules
//
include { SUMMARYLINE          } from '../modules/local/create-summaryline'
include { COMBINE_SUMMARYLINES } from '../modules/local/combine-summary'

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
include { PREPARE  } from '../subworkflows/local/prepare'
include { CLASSIFY } from '../subworkflows/local/classify'
include { ASSEMBLE } from '../subworkflows/local/assemble'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//
include { MULTIQC                     } from '../modules/nf-core/multiqc/main'
include { CUSTOM_DUMPSOFTWAREVERSIONS } from '../modules/nf-core/custom/dumpsoftwareversions/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Info required for completion email and summary
def multiqc_report = []

workflow VAPER {

    ch_versions = Channel.empty()

    /* 
    =============================================================================================================================
        PREPARE INPUT
    =============================================================================================================================
    */
    // SUBWORKFLOW: Read in samplesheet, validate and stage input files
    PREPARE (
        file(params.input),
        file(params.refs)
    )
    ch_versions = ch_versions.mix(PREPARE.out.versions)

    PREPARE.out.reads.set{ ch_reads }
    PREPARE.out.reads_qc.set{ ch_reads_qc }
    PREPARE.out.refs_man.set{ ch_refs_man }
    PREPARE.out.refs.set{ ch_refs }

    /* 
    =============================================================================================================================
        CLASSIFY VIRUSES
    =============================================================================================================================
    */ 

    //
    // SUBWORKFLOW: Classify viruses
    //
    CLASSIFY (
        ch_reads,
        ch_refs,
        ch_refs_man
    )
    ch_versions = ch_versions.mix(CLASSIFY.out.versions)

    /* 
    =============================================================================================================================
        ASSEMBLE VIRUSES
    =============================================================================================================================
    */
    // SUBWORKFLOW: Create consensus assemblies
    ASSEMBLE (
        CLASSIFY.out.ref_list.combine(ch_reads, by: 0)
    )
    ch_versions = ch_versions.mix(ASSEMBLE.out.versions)

    /* 
    =============================================================================================================================
        SUMMARIZE RESULTS
    =============================================================================================================================
    */

    // Combine outputs of samples that had references
    ASSEMBLE
        .out
        .bam_stats
        .join(ASSEMBLE.out.nextclade, by: [0,1])
        .set{ ch_assembly_list }
    // Create channel of samples with no reference
    if(! CLASSIFY.out.ref_list ){ println "No reference, no assembly!" }
    ch_reads
        .map{ meta, reads -> [ meta ] }
        .join( CLASSIFY.out.ref_list.map{ meta, ref_id, ref -> [ meta, ref_id ] }.ifEmpty( [ null, null ] ), by: 0, remainder: true )
        .filter{ meta, ref_id -> ref_id == null}
        .map{ meta, ref_id -> [ meta, "No_Reference", [], [] ] }
        .set{ ch_no_assembly_list }

    // Combine the reference and non-reference channels & add read QC & Sourmash results
    ch_assembly_list
        .concat(ch_no_assembly_list)
        .combine(ch_reads_qc.map{ sample, fastqc, fastp_raw, fastp_clean -> [ sample, fastp_raw, fastp_clean ] }, by: 0)
        .combine(CLASSIFY.out.sm_summary, by: 0)
        .set{ all_list }

    // MODULE: Create summaryline for each sample 
    SUMMARYLINE (
       all_list.combine( CLASSIFY.out.refsheet )
    )
    ch_versions = ch_versions.mix(SUMMARYLINE.out.versions)

    // MODULE: Combine summarylines
    SUMMARYLINE
        .out
        .summaryline
        .map{ meta, summaryline -> [ summaryline ] }
        .collect()
        .set{ all_summaries }
    
    COMBINE_SUMMARYLINES (
        all_summaries
    )
    ch_versions = ch_versions.mix(COMBINE_SUMMARYLINES.out.versions)

    CUSTOM_DUMPSOFTWAREVERSIONS (
        ch_versions.unique().collectFile(name: 'collated_versions.yml')
    )

    /* 
    =============================================================================================================================
        DEFAULTS
    =============================================================================================================================
    */

    // MODULE: MultiQC
    workflow_summary    = WorkflowVaper.paramsSummaryMultiqc(workflow, summary_params)
    ch_workflow_summary = Channel.value(workflow_summary)

    methods_description    = WorkflowVaper.methodsDescriptionText(workflow, ch_multiqc_custom_methods_description, params)
    ch_methods_description = Channel.value(methods_description)

    ch_multiqc_files = Channel.empty()
    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect())
    ch_multiqc_files = ch_multiqc_files.mix(ch_reads_qc.collect{it[1]}.ifEmpty([]))

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList()
    )
    multiqc_report = MULTIQC.out.report.toList()
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COMPLETION EMAIL AND SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow.onComplete {
    if (params.email || params.email_on_fail) {
        NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, multiqc_report)
    }
    NfcoreTemplate.dump_parameters(workflow, params)
    NfcoreTemplate.summary(workflow, params, log)
    if (params.hook_url) {
        NfcoreTemplate.IM_notification(workflow, params, summary_params, projectDir, log)
    }

    // animalFacts()
}

// Random fact bot
def animalFacts(){
    try {    
        def facts = file("https://github.com/ekohrt/animal-fun-facts-dataset/raw/refs/heads/main/animal-fun-facts-dataset.csv").splitCsv(header: true)
        randomFact = facts[new Random().nextInt(facts.size())]
        factMessage = "[${randomFact.animal_name}] ${randomFact.text.replaceAll('"','')} (source: ${randomFact.source})"
    } catch (Exception e) {
        factMessage = "No facts for you!"
    }
    println "Animal Fact: ${factMessage}"
}
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
