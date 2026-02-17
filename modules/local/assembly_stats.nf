process ASSEMBLY_STATS {
    tag "${prefix}"
    maxRetries 1
    
    input:
    tuple val(meta), val(ref_id), path(consensus), path(ref)

    output:
    tuple val(meta), val(ref_id), path("*.json"), emit: json
    path "versions.yml",                          emit: versions
    
    script:
    prefix = "${meta.id}_${ref_id}"
    tool = "vaper_stats.py"
    """
    ${tool} \\
        --ref_name ${ref_id} \\
        --query_name ${meta.id} \\
        --ref ${ref} \\
        --query ${consensus} \\
        ${task.attempt > 1 ? '--memsave' : ''}

        
    # version info
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ${tool}: "\$(${tool} --version 2>&1 | tr -d '\\r')"
    END_VERSIONS
    """
}
