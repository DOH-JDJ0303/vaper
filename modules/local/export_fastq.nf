process EXPORT_FASTQ {
    tag "${meta.id}"
    label 'process_single'

    conda "bioconda::seqtk=1.4"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/seqtk:1.4--he4a0461_1' :
        'biocontainers/seqtk:1.4--he4a0461_1' }"

    input:
    tuple val(meta), path(read_list), path(reads)

    output:
    path "*.gz"         , emit: reads
    path "versions.yml" , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args   ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    for FILE in ${read_list}; do
        BN=\$(basename \${FILE} | sed 's/.read-list.txt//g')
        # forward reads
        seqtk \\
            subseq \\
            $args \\
            ${reads[0]} \\
            "\${FILE}"  | \\
            gzip --no-name > "\${BN}_R1.fastq.gz"

        # reverse reads
        seqtk \\
            subseq \\
            $args \\
            ${reads[1]} \\
            "\${FILE}"  | \\
            gzip --no-name > "\${BN}_R2.fastq.gz"
    done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        seqtk: \$(echo \$(seqtk 2>&1) | sed 's/^.*Version: //; s/ .*\$//')
    END_VERSIONS
    """
}
