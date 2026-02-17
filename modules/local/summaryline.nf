process SUMMARYLINE {
    tag "$meta.id"
    label 'process_low'

    input:
    tuple val(meta), path(fastp), path(metagenome), path(ref_info), path(bam_stats), path(assembly_stats)

    output:
    tuple val(meta), path("*.json.gz"), emit: json
    path "versions.yml",                emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    tool = 'vaper_summary.py'
    """
    ${tool} \\
        --max_depth ${params.cons_max_depth} \\
        --qc_depth_min ${params.qc_depth} \\
        --qc_gf_min ${params.qc_genfrac} \\
        --sample ${meta.id} \\
        ${fastp ? '--fastp ' + fastp : ''} \\
        ${metagenome ? '--metagenome ' + metagenome : ''} \\
        ${ref_info ? '--ref_info ' + ref_info.join(' ') : ''} \\
        ${bam_stats ? '--bam_stats ' + bam_stats.join(' ') : ''} \\
        ${assembly_stats ? '--assembly_stats ' + assembly_stats.join(' ') : ''}

    # version info
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ${tool}: "\$(${tool} --version 2>&1 | tr -d '\\r')"
    END_VERSIONS
    """
}
