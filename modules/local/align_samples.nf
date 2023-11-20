process ALIGN_SAMPLES {
    container "docker.io/biocontainers/bowtie2:v2.4.1_cv1"
    label 'process_high'

    input:
        tuple val(meta), path(reads), path(summary)
        path(index_path)
    output:
        tuple val(meta), path("*_human.sam"), path("${meta.id}.sam"), path("${meta.id}_summary.csv"), emit: sams

    script:
    """
    #!/bin/bash

    bowtie2 -p ${task.cpus} -x ${index_path}/GRCh37/GRCh37 -1 ${reads[0]} -2 ${reads[1]} --local --no-unal -p 10 -S ${meta.id}_human.sam 
    bowtie2 -p ${task.cpus} -x ${index_path}/${meta.reference}/${meta.reference} -1 ${reads[0]} -2 ${reads[1]} --local --no-unal -p 6 -S ${meta.id}.sam
    """
}
