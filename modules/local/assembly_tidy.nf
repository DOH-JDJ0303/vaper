process ASSEMBLY_TIDY {
    tag "${prefix}"
    
    input:
    tuple val(meta), val(ref_id), path(consensus)

    output:
    tuple val(meta), val(ref_id), path("*.tidy.fa.gz"), emit: fasta
    path "versions.yml",                                emit: versions
    
    script:
    prefix = "${meta.id}_${ref_id}"
    tool = "vaper_tidy.py"
    """
    ${tool} \\
        ${consensus} \\
        ${params.cons_prune_termini ? '--no_terminal_n' : '' } \\
        ${params.cons_no_mixed_sites ? '--no_mixed_sites' : '' }
        
        
    # version info
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ${tool}: "\$(${tool} --version 2>&1 | tr -d '\\r')"
    END_VERSIONS
    """
}
