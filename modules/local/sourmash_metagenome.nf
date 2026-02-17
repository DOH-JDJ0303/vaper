process SOURMASH_METAGENOME {
    tag "$meta.id"
    label 'process_single'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/sourmash:4.9.4--hdfd78af_0':
        'biocontainers/sourmash:4.9.4--hdfd78af_0' }"

    input:
    tuple val(meta), path(gather_results)
    path(taxonomy)

    output:
    tuple val(meta), path('*.csv'),       emit: result
    tuple val(meta), path('*.krona.tsv'), emit: krona
    path "versions.yml"           ,       emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    sourmash tax metagenome \\
        $args \\
        -o ${prefix} \\
        -g ${gather_results} \\
        -t ${taxonomy} \\
        -r species

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sourmash: \$(echo \$(sourmash --version 2>&1) | sed 's/^sourmash //' )
    END_VERSIONS
    """
}
