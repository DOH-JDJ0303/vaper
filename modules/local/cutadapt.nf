process CUTADAPT {
    tag "$meta.id"
    label 'process_medium'

    conda "cutadapt"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/cutadapt:5.0--py39hbcbf7aa_0':
        'biocontainers/cutadapt:5.0--py310h1fe012e_0' }"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.fastq.gz"), emit: reads
    path "versions.yml",                 emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    cutadapt \\
        $args \\
        -o ${prefix}-cutadapt_1.fastq.gz \\
        -p ${prefix}-cutadapt_2.fastq.gz \\
        -j $task.cpus \\
        ${reads[0]} \\
        ${reads[1]}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cutadapt: \$(cutadapt --version)
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    touch ${prefix}.merged.fq.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        NGmerge: \$(echo \$(NGmerge --version 2>&1) | sed 's/^.*NGmerge, version //; s/ Copyright.*// ; s/: //g' ))
    END_VERSIONS
    """
}
