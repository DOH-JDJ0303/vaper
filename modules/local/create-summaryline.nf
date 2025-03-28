process SUMMARYLINE {
    tag "$meta.id"
    label 'process_low'

    input:
    tuple val(meta), val(ref_id), path(bam_stats), val(nextclade), path(fastp_raw, stageAs: 'raw.json'), path(fastp_clean, stageAs: 'clean.json'), path(sm_summary), path(refsheet)

    output:
    tuple val(meta), path("*.summaryline.csv"), emit: summaryline
    path "versions.yml",                        emit: versions


    when:
    task.ext.when == null || task.ext.when

    prefix = task.ext.prefix ?: "${meta.id}_${ref_id}"
    script:
    """
    # Format reports to tables
    ## Fastp
    fastp2tbl.sh ${fastp_raw} ${fastp_clean} > ${prefix}.fastp2tbl.csv
    ## Mapping stats
    if [ -f "${bam_stats}" ]
    then
        samtoolstats2tbl.sh ${bam_stats} > ${prefix}.samtoolstats2tbl.csv
    fi
    ## Nextclade
    if [ "${nextclade}" != "[]" ]
    then
        echo "${nextclade}" > nextclade.csv
    fi 
    
    # Extract refsheet
    ${refsheet.name.endsWith('.gz') ? 'zcat' : 'cat'} ${refsheet} > refsheet

    # Create summaryline
    summaryline.R ${prefix}.fastp2tbl.csv "${sm_summary}" ${prefix}.samtoolstats2tbl.csv nextclade.csv "${meta.id}" "${ref_id}" refsheet
    # rename using prefix and reference
    mv summaryline.csv "${prefix}.summaryline.csv"
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        summaryline.R: \$(summaryline.R version)
    END_VERSIONS
    """
}
