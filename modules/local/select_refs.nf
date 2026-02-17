process SELECT_REFS {
    tag "$meta.id"
    label 'process_low'

    input:
    tuple val(meta), path(assembly), val(include), val(exclude)
    path refs
    
    output:
    tuple val(meta), path("subset.fa.gz"),             emit: fa,    optional: true
    tuple val(meta), path("subset.jsonl.gz"),          emit: jsonl, optional: true
    tuple val(meta), path("genome-fraction.csv"), emit: gf, optional: true

    when:
    task.ext.when == null || task.ext.when

    script:

    prefix = task.ext.prefix ?: "${meta.id}"
    tool = "vaper_refs.py"
    """
    ${tool} \\
        --refs "${refs}" \\
        --query "${assembly}" \\
        --genfrac ${params.ref_genfrac} \\
        --dist ${params.ref_dist} \\
        ${include ? "--include '" + include.join(',') + "'" : ''} \\
        ${exclude ? "--exclude '" + exclude.join(',') + "'" : ''}

    # version info
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ${tool}: "\$(${tool} --version 2>&1 | tr -d '\\r')"
    END_VERSIONS
    """
}
