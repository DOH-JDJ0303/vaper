process PAF_TO_SAM {
    tag "$meta.id-$ref_id"
    label 'process_medium'

    conda "bioconda::minimap2=2.24 bioconda::samtools=1.17"
    container "${ workflow.containerEngine == 'singularity' ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-66534bcbb7031a148b13e2ad42583020b9cd25c4:1679e915ddb9d6b4abda91880c4b48857d471bd8-0' :
        'biocontainers/mulled-v2-66534bcbb7031a148b13e2ad42583020b9cd25c4:1679e915ddb9d6b4abda91880c4b48857d471bd8-0' }"

    input:
    tuple val(meta), val(ref_id), path(paf), path(reference), path(reads)

    output:
    tuple val(meta), val(ref_id), path("*.bam"), emit: bam
    tuple val(meta), val(ref_id), path("*.bai"), emit: bai
    path "versions.yml"                         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = "${meta.id}_${ref_id}"
    def reads_input = reads instanceof List ? reads.join(' ') : reads
    """
    # Re-run minimap2 with -a flag to get SAM output
    # PAF cannot be directly converted to SAM/BAM as it lacks sequence information
    minimap2 \\
        -a \\
        -x sr \\
        -t ${task.cpus} \\
        $args \\
        $reference \\
        $reads_input \\
        | samtools view -b -F 4 \\
        | samtools sort -@ ${task.cpus} -o ${prefix}.bam

    samtools index ${prefix}.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        minimap2: \$(minimap2 --version 2>&1)
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """

    stub:
    def prefix = "${meta.id}_${ref_id}"
    """
    touch ${prefix}.bam
    touch ${prefix}.bai
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        minimap2: \$(minimap2 --version 2>&1)
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}
