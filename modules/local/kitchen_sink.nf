process KITCHEN_SINK {
    label 'process_low'
    container 'public.ecr.aws/o8h2f0o1/ncbi-datasets:16.15.0'

    input:
    path summaries

    output:
    path "contigs/*",         emit: fa, optional: true
    path "metadata.jsonl.gz", emit: jsonl, optional: true
    path "manifest.csv",      emit: csv, optional: true
    path "versions.yml",      emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    args   = task.ext.args ?: ''
    tool = "vaper_sink.py"
    """
    ${tool} \\
        ${summaries}

    # version info
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ${tool}: "\$(${tool} --version 2>&1 | tr -d '\\r')"
    END_VERSIONS
    """
}


