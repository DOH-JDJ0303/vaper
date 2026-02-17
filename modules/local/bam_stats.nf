process BAM_STATS {
    tag "${prefix}"
    label 'process_high'

    conda "bioconda::bwa"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-fe8faa35dbf6dc65a0f7f5d4ea12e31a79f73e40:219b6c272b25e7e642ae3ff0bf0c5c81a5135ab4-0' :
        'biocontainers/mulled-v2-fe8faa35dbf6dc65a0f7f5d4ea12e31a79f73e40:219b6c272b25e7e642ae3ff0bf0c5c81a5135ab4-0' }"

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("${prefix}.coverage.txt"), emit: coverage
    tuple val(meta), path("*.stats.txt"),            emit: stats
    tuple val(meta), path("*.read-list.txt"),        emit: read_list
    path "versions.yml",                             emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = "${meta.id}"
    """
    # setup for pipe
    set -euxo pipefail

    samtools index ${bam}

    # gather read stats
    samtools coverage ${bam} > ${prefix}.coverage.txt

    for CONTIG in \$(samtools idxstats ${bam} | cut -f 1 | grep -v '\\*'); do
        samtools stats ${bam} "\$CONTIG" > "${prefix}-\${CONTIG}.stats.txt"
        samtools view ${bam} "\$CONTIG" | cut -f 1 > "${prefix}-\${CONTIG}.read-list.txt"
    done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools version 2>&1) | cut -f 2 -d ' ')
    END_VERSIONS
    """
}
