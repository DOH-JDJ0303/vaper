process FORMAT_REFS {
    label 'process_low'

    conda "bioconda::seqtk"
    container 'docker.io/staphb/seqtk:1.3'

    input:
    path refs_tar

    output:
    tuple path("refs.fa.gz"), path("refs-comp.txt.gz"), path("*.tar.gz", includeInputs: true), path("refsheet.csv.gz"), emit: refs
    path "versions.yml",  emit: versions

    when:
    task.ext.when == null || task.ext.when

    version = "1.0"
    script: // This script is bundled with the pipeline, in nf-core/waphlviral/bin/
    """
    # Extract references
    tar xvzf ${refs_tar}

    # Copy refsheet
    cp */*.csv refsheet.csv && gzip refsheet.csv

    # Combine sequences into single file with renamed fasta headers
    for ref in \$(ls */references/)
    do
        seqtk seq -C */references/\${ref} > tmp.fa
        seqtk rename tmp.fa \${ref##*/} >> refs.fa
    done
    rm tmp.fa
    gzip refs.fa

    # Get contig lengths
    seqtk comp refs.fa.gz | gzip > refs-comp.txt.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        seqtk: \$(seqtk 2>&1 | grep "Version" | cut -f 2 -d ' ')
    END_VERSIONS
    """
}
