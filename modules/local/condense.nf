process CONDENSE {
    tag "${prefix}"
    label "process_high"
    stageInMode 'copy'

    input:
    tuple val(meta), path(assemblies, stageAs: 'input/*'), path(read_stats, stageAs: 'input/*')

    output:
    tuple val(meta), path("*.fa.gz"), emit: assembly
    tuple val(meta), path("*.csv"),   emit: summary, optional: true
    path "versions.yml",              emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = "${meta.id}"
    tool = 'vaper_condense.py'
    """
    ${tool} \\
        --fasta ${assemblies} \\
        --stats ${read_stats} \\
        --prefix ${prefix}

    # version info
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ${tool}: "\$(${tool} --version 2>&1 | tr -d '\\r')"
    END_VERSIONS
    """
}