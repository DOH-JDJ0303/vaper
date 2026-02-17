process SRA_HUMAN_SCRUBBER {
    tag "${meta.id}"
    
    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.scrubbed.fastq.gz"), emit: reads
    path "versions.yml",                          emit: versions
    
    script:
    """
    # interleave paired reads or just decompress
    ${meta.single_end ? 'seqtk seq ${reads} > reads.fq' : "seqtk mergepe ${reads} > reads.fq" }

    # remove human reads
    scrub.sh \\
        -i reads.fq \\
        -o scrubbed.fq \\
        -r -u 'human_reads' ${meta.single_end ? '-s' : ''}

    # split interleaved reads or rename scrubbed.fq
    ${meta.single_end ? 'mv scrubbed.fq ${meta.id}.scrubbed.fastq' : "seqtk seq -1 scrubbed.fq > ${meta.id}_R1.scrubbed.fastq" }
    ${meta.single_end ? '' : "seqtk seq -2 scrubbed.fq > ${meta.id}_R2.scrubbed.fastq" }

    # compress reads
    pigz -p ${task.cpus} *.scrubbed.fastq
        
    # version info
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sra-human-scrubber: 2.2.1
        seqtk: \$(echo \$(seqtk 2>&1) | sed 's/^.*Version: //; s/ .*\$//')

    END_VERSIONS
    """
}
