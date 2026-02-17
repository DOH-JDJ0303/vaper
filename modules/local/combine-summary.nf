process COMBINE_SUMMARYLINES {
    label 'process_low'

    input:
    path summaries
    path params_json

    output:
    path "VAPER-summary.json.gz", emit: json
    path "VAPER-summary.csv",     emit: csv
    path "versions.yml",          emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    tool = 'combine_summary.py'
    """
    
    ${tool} \\
        --params ${params_json} \\
        --results ${summaries}

    # version info
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ${tool}: "\$(${tool} --version 2>&1 | tr -d '\\r')"
    END_VERSIONS
    """
}
