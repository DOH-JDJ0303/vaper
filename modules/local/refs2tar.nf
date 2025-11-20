process REFS2TAR {
    tag "${refsheet}"
    label 'process_low'
    stageInMode 'copy'

    input:
    tuple path(refsheet), path(assemblies)

    output:
    path "refset.tar.gz",  emit: tar

    when:
    task.ext.when == null || task.ext.when

    script: // This script is bundled with the pipeline, in nf-core/waphlviral/bin/
    """
    mkdir refset || true
    mv * refset || true
    tar -cvzf refset.tar.gz refset/
    """
}
