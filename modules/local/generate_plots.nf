process GENERATE_PLOTS {
    container "docker.io/jefffurlong/r-ggplot"
    label 'process_high'

    input:
        tuple val(meta), path(inserts)
        path(density)
    output:
        path("*.png")

    script:
    """
    #!/bin/bash
    Rscript ${density} ${meta.sample_name} ${meta.reference} ${meta.id}_dedup_inserts.txt ${meta.id}_human_dedup_inserts.txt
    """
}
