process TAR2REFS {
    tag "${refs_tar}"
    label 'process_low'
    stageInMode 'copy'

    input:
    tuple val(ref_ids), path(refs_tar)

    output:
    path "select/*",    emit: refs

    when:
    task.ext.when == null || task.ext.when

    script: // This script is bundled with the pipeline, in nf-core/waphlviral/bin/
    """
    # Save selected references to file
    echo -e "${ref_ids.join('\n')}" > select.csv

    # Extract references
    tar -xzf ${refs_tar}

    # Move over reference files
    mkdir select
    for REF in \$(cat select.csv)
    do
        cp */references/\${REF} select/
    done
    """
}
