process VALIDATE_REFS {
    label 'process_low'

    input:
    path refs, stageAs: 'input/'

    output:
    path "refs.valid.jsonl.gz", emit: jsonl
    path "versions.yml",        emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:

    tool = "vaper_refs.py"
    """
    # compress references
    pigz -p ${task.cpus} input/*.jsonl || true

    # combine references
    cat input/*.jsonl.gz > refs.valid.jsonl.gz

    # validate
    ${tool} \\
        --validate \\
        --refs refs.valid.jsonl.gz

    # version info
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ${tool}: "\$(${tool} --version 2>&1 | tr -d '\\r')"
    END_VERSIONS
    """
}
