process PROCESS_BAMS {
    container "docker.io/jefffurlong/samtools_picard"
    label 'process_medium'

    input:
       tuple val(meta), path(human_sam), path(reference_sam), path(summary)
    output:
        path("*.bam") 
        path("*final_summary.csv")
        path("*dedup_inserts.txt")

    script:
    """
    #!/bin/bash

    samtools view -@ 6 -Sb ${reference_sam} | samtools sort -o ${meta.id}_sorted.bam 
    samtools view -@ 6 -Sb ${human_sam} | samtools sort -o ${meta.id}_human_sorted.bam

    java -jar /usr/bin/picard.jar MarkDuplicates INPUT=${meta.id}_sorted.bam OUTPUT=${meta.id}_deduplicated.bam VERBOSITY=ERROR REMOVE_DUPLICATES=TRUE VALIDATION_STRINGENCY=SILENT ASSUME_SORTED=TRUE METRICS_FILE=${meta.id}_deduplicated_metrics.txt
    java -jar /usr/bin/picard.jar MarkDuplicates INPUT=${meta.id}_human_sorted.bam OUTPUT=${meta.id}_human_deduplicated.bam VERBOSITY=ERROR REMOVE_DUPLICATES=TRUE VALIDATION_STRINGENCY=SILENT ASSUME_SORTED=TRUE METRICS_FILE=${meta.id}_human_deduplicated_metrics.txt

    human_bam_reads=\$(samtools flagstat ${meta.id}_human_sorted.bam | grep "mapped (" | awk '{print \$1}')
    human_dedup_bam_reads=\$(samtools flagstat ${meta.id}_human_deduplicated.bam | grep "mapped (" | awk '{print \$1}')
    bam_reads=\$(samtools flagstat ${meta.id}_sorted.bam | grep "mapped (" | awk '{print \$1}')
    dedup_bam_reads=\$(samtools flagstat ${meta.id}_deduplicated.bam | grep "mapped (" | awk '{print \$1}')

    samtools view -f66 ${meta.id}_deduplicated.bam | cut -f9 | awk '{print sqrt(\$0^2)}' > ${meta.id}_dedup_inserts.txt

    cp ${summary} ${meta.id}_final_summary.csv
    printf "\$human_bam_reads,\$human_dedup_bam_reads,\$bam_reads,\$dedup_bam_reads" >> ${meta.id}_final_summary.csv

    """
} 