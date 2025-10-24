process PREPEND_TSV_WITH_ID {
    tag "$tsv"
    label 'process_low'

    // we just need a base linux environment for this module
    // which is an assumption of the nextflow pipeline

    input:
    tuple val(meta), path(tsv)                                                

    output:
    tuple val(meta), path ("*prepended*", includeInputs: false) , emit: prepended
    path ("*prepended*"), includeInputs: false                  , emit: prepended_file
    path "versions.yml"                                         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:

    def extension = tsv.extension
    def new_name  = tsv.name.replaceAll(/\Q${extension}\E$/, "prepended.${extension}")

    // prepend tsv output with a column containing sample ID (from meta.id)
    // ignore lines beginning with # (comment lines)
    """
    grep -v -e "^#" ${tsv} | awk '{print "${meta.id}" "\t" \$0}' > ${new_name}
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
    END_VERSIONS
    """
}
