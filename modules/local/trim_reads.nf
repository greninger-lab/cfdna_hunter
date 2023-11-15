process TRIM_READS { 

    container "docker.io/biocontainers/fastp:v0.20.1_cv1"
    label 'process_medium'

    input:
        tuple val(meta), path(reads)
    output:
        tuple val(meta), path("*_37bp.fastq.gz"), path("*.csv"), emit: trimmed_reads

    
    script:
    """
    #!/bin/bash

    fastp -w ${task.cpus} -i ${reads[0]} -I ${reads[1]} -o ${meta.id}_R1_37bp.fastq.gz -O ${meta.id}_R2_37bp.fastq.gz -b 37 -B 37 --detect_adapter_for_pe > fastp_out.txt 2>&1

    trimmed_reads=\$(grep "reads passed filter:" fastp_out.txt | cut -f2 -d":" | cut -f2 -d" ")

    echo "sample_name,raw_reads,trimmed_reads,human_reads,dedup_human_reads,virus_reads,dedup_virus_reads" > ${meta.id}_summary.csv
    num_r1_untrimmed=\$(gunzip -c ${reads[0]} | wc -l)
    num_r2_untrimmed=\$(gunzip -c ${reads[1]} | wc -l)
    num_untrimmed=\$((\$((num_r1_untrimmed + num_r2_untrimmed))/4))
    printf "${meta.id},\$num_untrimmed,\$trimmed_reads," >> ${meta.id}_summary.csv

    """
}