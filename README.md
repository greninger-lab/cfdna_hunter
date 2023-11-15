# cfdna_hunter
Nextflow pipeline for searching for non human genomic material in cfDNA. 

## Usage
Install [`Nextflow`](https://www.nextflow.io/docs/latest/getstarted.html#installation) (`>=21.10.3`)

Install [`Docker`](https://docs.docker.com/engine/installation/)


### Command line:
    nextflow run greninger-lab/cfdna_hunter \
        --input PATH_TO_SAMPLE_CSV \                      # required
        --outdir PATH_TO_OUTPUT_FOLDER \                  # required
        --index_path PATH_TO_BOWTIE2_INDEX_ROOT_FOLDER \  # required (path to folders containing bowtie2 indexes e.g. /my_indexes/NC_001798/NC_001798.1.bt2 etc)
        -profile docker \                                 # required
        -with-tower \                                     # optional (use if you want to use Nextflow Tower)
        -c nextflow_aws.config \                          # optional (AWS account config info) 
        -r main                                           # required (use the github main branch)

#### Input sample csv format (note: the pairing of "sample" and "reference" must be unique):
---------
   sample,fastq_1,fastq_2,reference
   SRR1024623,/Users/jfurlong/dev/cfDNA/SRR1024623_R1.fastq.gz,/Users/jfurlong/dev/cfDNA/SRR1024623_R2.fastq.gz,FR751470
---------

