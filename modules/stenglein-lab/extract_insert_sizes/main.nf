/*
  This process extracts insert sizes from samtools stats output
 */
process EXTRACT_INSERT_SIZES {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.18--h50ea8bc_1' :
        'quay.io/biocontainers/samtools:1.18--h50ea8bc_1' }"

    input:
    tuple val(meta), path(bam), val(refseq_id)

    when:
    task.ext.when == null || task.ext.when

    output:
    tuple val(meta), path("*.insert_sizes.txt"), emit: insert_sizes

    shell:
    def new_name = bam.name.replaceAll(/.bam$/, ".insert_sizes.txt")

    """
    # run samtools stats and extract insert sizes
    # and prepend with sample ID and refseq ID
    samtools stats $bam | grep "^IS" | cut -f 2- | awk '{print "${meta.id}" "\t" "$refseq_id.id" "\t" \$0}' > ${new_name}
    """
}
