process REF_CSV2JSONL {
    label 'process_low'

    input:
    path refsheet
    path assemblies

    output:
    path "*.jsonl.gz",   emit: jsonl
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:

    tool = "ref_csv2jsonl.py"
    """
    ${tool} ${refsheet}

    # version info
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ${tool}: "\$(${tool} --version 2>&1 | tr -d '\\r')"
    END_VERSIONS
    """
}
